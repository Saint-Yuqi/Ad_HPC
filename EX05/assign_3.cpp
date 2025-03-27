// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <new>
#include <complex>
#include <stdlib.h>
#include "blitz/array.h"
#include "fftw3.h"
#include "tipsy.h"
#include "aweights.h"
using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// A separate version is created for each different "Order".
// This allows the compiler to optimize the process for each of the four orders
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R,Array<float,1> &M) {
    auto nGrid = grid.rows(); // total number of particles
    // C++ Lambda to apply the periodic wrap of the grid index
    auto wrap = [nGrid](int i) {
        if (i<0) i+=nGrid;
        else if (i>=nGrid) i-=nGrid;
        return i;
    };
    #pragma omp parallel for
    for(int pn=0; pn<R.rows(); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),Hy((y+0.5f)*nGrid),Hz((z+0.5f)*nGrid);
        for(auto i=0; i<Order; ++i) {
            for(auto j=0; j<Order; ++j) {
                for(auto k=0; k<Order; ++k) {
                    #pragma omp atomic
                    grid(wrap(Hx.i+i),wrap(Hy.i+j),wrap(Hz.i+k)) += m * Hx.H[i]*Hy.H[j]*Hz.H[k];
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    std::locale::global(std::locale("")); // e.g., LC_ALL=en_GB.UTF-8 or de_CH.UTF-8
    std::cerr.imbue(std::locale()); // e.g., LC_ALL=en_GB.UTF-8
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size] [order]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100;
    if (argc>2) nGrid = atoi(argv[2]);

    int iOrder = 1;
    if (argc>3) iOrder = atoi(argv[3]);

    // Read the position and masses from the file
    auto t0 = hrc::now();
    std::ifstream io(argv[1],std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h),sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << std::fixed << N << " particles" << std::endl;
    std::cerr << "Loading " << N << " particles" << std::endl;
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // Load the particles
    tipsy::dark d;
    for(int i=0; i<N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d),sizeof(d))) {
            std::cerr << "error reading particle" << std::endl;
            return errno;
        }
        r(i,0) = d.pos[0];
        r(i,1) = d.pos[1];
        r(i,2) = d.pos[2];
        m(i) = d.mass;
    }

    duration dt = hrc::now() - t0;
    std::cerr << "File reading took " << std::setw(9) << dt.count() << " seconds." << std::endl;

    // Create Mass Assignment Grid
    t0 = hrc::now();

    // 1: Calculate the number of floats to allocate.
    // 1: Use size_t to avoid integer overflow for large grid sizes
    // 1: Allocate the storage ourselves
    // 2: Allocate n x n x n+2 to account for FFT padding
    // 2: Create the array with the new shape
    // 3: Create a subarray from the "raw" array
    auto k_nz = nGrid/2 + 1;
    auto n_floats = size_t(1) * nGrid * nGrid * 2*k_nz; // Careful. For odd nGrid 2*k_nz != nGrid + 2
    float *data = new (std:: align_val_t (64)) float [n_floats]; // 512-bit alignment
    Array<float,3> raw_grid(data,shape(nGrid,nGrid,2*k_nz),deleteDataWhenDone);
    Array<float,3> grid(raw_grid(Range(0,nGrid-1),Range(0,nGrid-1),Range(0,nGrid-1)));

    // 4: Create a complex view of the data to match the inplace memory order
    // 4: memory policy is neverDeleteData -- only one of the Array should delete []
    Array <std::complex<float>,3> kgrid(reinterpret_cast<std::complex<float>*>(data),shape(nGrid,nGrid,k_nz),neverDeleteData);

    // Assign the mass to the grid
    grid = 0;

    // This creates four different versions of "assign mass", one for each order
    std::cerr << "Assigning mass to the grid using order " << iOrder <<std::endl;
    switch(iOrder) {
    case 1:
        assign_mass<1>(grid,r,m);
        break;
    case 2:
        assign_mass<2>(grid,r,m);
        break;
    case 3:
        assign_mass<3>(grid,r,m);
        break;
    case 4:
        assign_mass<4>(grid,r,m);
        break;
    default:
        std::cerr << "Invalid order " << iOrder << " (must be 1, 2, 3 or 4)" << std::endl;
    }
    dt = hrc::now() - t0;
    std::cerr << "Mass assignment took " << std::setw(9) << dt.count() << " seconds." << std::endl;
    // Compute mean density
    /*The first exercise of HW5*/
    float total_mass = sum(grid);  // Blitz handles sum over whole array
    float mean_density = total_mass; // Since box volume is 1

    // Convert mass to density contrast δ(r)
    Array<float,3> delta(shape(nGrid, nGrid, nGrid));
    delta = grid / mean_density - 1.0f;
 
    // Calculate projected density
    t0 = hrc::now();
    Array<float,2> projected(nGrid,nGrid);
    thirdIndex ii;
    projected = max(grid,ii);
    dt = hrc::now() - t0;
    std::cerr << "Density projection took " << std::setw(9) << dt.count() << " seconds." << std::endl;

    // Write out the 2D map
    std::ofstream of("density.dat",std::ios::binary);
    of.write(reinterpret_cast<char*>(projected.data()),projected.size()*sizeof(float));
    std::cerr << "Finished writing density.dat successfully.\n";

    // Calculate the FFT
    auto plan = fftwf_plan_dft_r2c_3d(nGrid, nGrid, nGrid,
            delta.dataFirst(),  // ✅ use δ(r)
            reinterpret_cast<fftwf_complex*>(kgrid.dataFirst()),
            FFTW_ESTIMATE);

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    std::cerr << "FFT completed with no crash.\n";
    
    //Exercise 3
    // Let's define how many bins we want
    int nbins = 80;

    // 1) Determine kmax
    //    In a real-to-complex FFT with size nGrid, the largest frequency in each dimension
    //    is nGrid/2 in magnitude. So the largest possible |k| is sqrt((nGrid/2)^2 * 3).
    float kmax = std::sqrt( float(nGrid/2) * float(nGrid/2) * 3.0f );
    // If your box size is 1, then this wave number is dimensionless (like your earlier approach).
    // If you also do 2π factors in real codes, you would multiply by that factor, etc.
    // For now, we assume k = sqrt(kx^2 + ky^2 + kz^2) exactly as in your previous code.

    // 2) The bin width = kmax / nbins (assuming kmin=0).
    float delta_bin = kmax / float(nbins);

    // Prepare accumulation arrays
    std::vector<float> fPowerVar(nbins, 0.0f);
    std::vector<int>   nPowerVar(nbins, 0);

    // Loop over k-space grid
    //int k_nz = nGrid/2 + 1;
    for (int ix = 0; ix < nGrid; ix++) {
        // Convert ix to kx in [-nGrid/2, +nGrid/2]
        int kx = (ix <= nGrid/2) ? ix : (ix - nGrid);

        for (int iy = 0; iy < nGrid; iy++) {
            int ky = (iy <= nGrid/2) ? iy : (iy - nGrid);

            for (int iz = 0; iz < k_nz; iz++) {
                int kz = iz;  // 0..nGrid/2

                // Complex δ(k)
                std::complex<float> val = kgrid(ix, iy, iz);

                // Power = |δ(k)|^2
                float p_k = std::norm(val);

                // Magnitude of wave vector
                float k_abs = std::sqrt(float(kx*kx + ky*ky + kz*kz));

                // Convert k_abs to a bin index iBin
                if (k_abs > 0.0f) {
                    int iBin = int(std::floor(k_abs / delta_bin));
                    if (iBin < nbins) {
                        fPowerVar[iBin] += p_k;
                        nPowerVar[iBin]++;
                    }
                }
            }
        }
    }

    std::ofstream outfile("power_variable_3.dat");
    outfile << "# binCenter  P(k)\n";

    for (int i = 0; i < nbins; i++) {
        float k_center = (i + 0.5f) * delta_bin;
        if (nPowerVar[i] > 0) {
            float avgPower = fPowerVar[i] / float(nPowerVar[i]);
            outfile << k_center << " " << avgPower << "\n";
        } else {
            outfile << k_center << " 0.0\n";
        }
    }
    outfile.close();
}

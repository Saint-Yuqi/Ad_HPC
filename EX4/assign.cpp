// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <new>    // for std::align_val_t
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <stdlib.h>
#include "blitz/array.h"
#include "tipsy.h"
#include "aweights.h"
#include <fftw3.h>  // Single-precision FFTW: fftwf_*
using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// A separate version is created for each different "Order".
// This allows the compiler to optimize the process for each of the four orders
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R,Array<float,1> &M) {
    auto nGrid = grid.rows(); // number of points along the first dimension
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
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),
                                      Hy((y+0.5f)*nGrid),
                                      Hz((z+0.5f)*nGrid);

        for(auto i=0; i<Order; ++i) {
            for(auto j=0; j<Order; ++j) {
                for(auto k=0; k<Order; ++k) {
                    #pragma omp atomic
                    grid(wrap(Hx.i+i), wrap(Hy.i+j), wrap(Hz.i+k))
                        += m * Hx.H[i] * Hy.H[j] * Hz.H[k];
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    std::locale::global(std::locale("")); // e.g., LC_ALL=en_GB.UTF-8 or de_CH.UTF-8
    std::cerr.imbue(std::locale());       // e.g., LC_ALL=en_GB.UTF-8
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size] [order]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100;
    std::size_t totalSize = std::size_t(nGrid) * std::size_t(nGrid) * std::size_t(nGrid + 2);

    // Allocate single-precision float memory
    float* data = new (std::align_val_t(64)) float[ totalSize ];

    // Create a Blitz array of shape (nGrid, nGrid, nGrid+2) for the real domain
    blitz::TinyVector<int,3> shapeReal(nGrid, nGrid, nGrid + 2);
    blitz::Array<float,3> padded(data, shapeReal, blitz::deleteDataWhenDone);

    // Also create a Blitz array of shape (nGrid, nGrid, nGrid/2+1) for the complex domain
    // (typical for r2c transforms). We cast the pointer to std::complex<float>*:
    blitz::TinyVector<int,3> shapeCplx(nGrid, nGrid, (nGrid/2 + 1));
    blitz::Array<std::complex<float>,3> kdata(
    reinterpret_cast<std::complex<float>*>(data), 
    shapeCplx, 
    blitz::neverDeleteData
    );

    // 1) Create a 3D plan for real-to-complex in-place transform
    fftwf_plan plan = fftwf_plan_dft_r2c_3d(
        nGrid, nGrid, nGrid,
        reinterpret_cast<float*>(padded.data()), 
        reinterpret_cast<fftwf_complex*>(kdata.data()),
        FFTW_ESTIMATE
    );

    // 2) Execute the plan
    fftwf_execute(plan);

    // Now 'kdata' contains the frequency-domain data in place

    // 3) Destroy the plan
    fftwf_destroy_plan(plan);

    if (argc>2) nGrid = atoi(argv[2]);

    int iOrder = 1;
    if (argc>3) iOrder = atoi(argv[3]);

    // Read the position and masses from the file
    auto t0 = hrc::now();
    std::ifstream io(argv[1], std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h), sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << std::fixed << N << " particles" << std::endl;
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // Load the particles
    tipsy::dark d;
    for(int i=0; i<N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
            std::cerr << "error reading particle" << std::endl;
            return errno;
        }
        r(i,0) = d.pos[0];
        r(i,1) = d.pos[1];
        r(i,2) = d.pos[2];
        m(i) = d.mass;
    }

    duration dt = hrc::now() - t0;
    std::cerr << "File reading took " << std::setw(9)
              << dt.count() << " seconds." << std::endl;

    // Create Mass Assignment Grid
    t0 = hrc::now();

    // 1) Allocate enough space for N x N x (N+2)
    //    (We add +2 in the last dimension for in-place FFT padding.)
    std::size_t totalSize = std::size_t(nGrid) * std::size_t(nGrid) * std::size_t(nGrid + 2);

    // 2) Allocate with 64-byte alignment (C++17)
    float* memGrid = new (std::align_val_t(64)) float[ totalSize ];

    // 3) Blitz++ shape is now (nGrid, nGrid, nGrid+2)
    blitz::TinyVector<int,3> shape(nGrid, nGrid, nGrid + 2);

    // 4) Construct the array with the external pointer, telling Blitz++ to delete[] it internally
    blitz::Array<float,3> grid(memGrid, shape, blitz::deleteDataWhenDone);

    // Initialize the entire array to zero
    grid = 0.0f;

    //---------------------------------------------------------
    // Create a subarray view of shape N x N x N
    // subGrid references only the "true" portion of the data.
    //---------------------------------------------------------
    blitz::Array<float,3> subGrid = 
        grid(Range::all(), Range::all(), Range(0, nGrid-1));

    // Assign the mass to the grid
    std::cerr << "Assigning mass to the grid using order " << iOrder << std::endl;
    switch(iOrder) {
        case 1:
            assign_mass<1>(subGrid, r, m);
            break;
        case 2:
            assign_mass<2>(subGrid, r, m);
            break;
        case 3:
            assign_mass<3>(subGrid, r, m);
            break;
        case 4:
            assign_mass<4>(subGrid, r, m);
            break;
        default:
            std::cerr << "Invalid order " << iOrder 
                    << " (must be 1, 2, 3 or 4)" << std::endl;
        }
    dt = hrc::now() - t0;
    std::cerr << "Mass assignment took "
              << std::setw(9) << dt.count() << " seconds." << std::endl;

    // Calculate projected density
    t0 = hrc::now();
    Array<float,2> projected(nGrid, nGrid);
    thirdIndex ii;
    projected = max(grid, ii);
    dt = hrc::now() - t0;
    std::cerr << "Density projection took "
              << std::setw(9) << dt.count() << " seconds." << std::endl;

    // Write out the 2D map
    std::ofstream of("density.dat", std::ios::binary);
    of.write(reinterpret_cast<char*>(projected.data()),
             projected.size() * sizeof(float));

    // The allocated memory (memGrid) will be freed automatically
    // by Blitz++ (deleteDataWhenDone policy) when grid goes out of scope.
    return 0;
}

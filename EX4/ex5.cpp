#include <iostream>
#include <new>                 // for std::align_val_t
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <cstdlib>             // for atoi
#include <complex>             // for std::complex<float>
#include <fftw3.h>            // single-precision FFTW
#include "blitz/array.h"
#include "tipsy.h"
#include "aweights.h"

using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// --------------------------------------------------------
// Template for mass assignment
// --------------------------------------------------------
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R, Array<float,1> &M)
{
    // dimension along the first axis
    int nGrid = grid.extent(0);

    // Lambda for periodic wrapping
    auto wrap = [nGrid](int i) {
        if (i < 0)      i += nGrid;
        else if (i>=nGrid) i -= nGrid;
        return i;
    };

    #pragma omp parallel for
    for(int pn = 0; pn < R.rows(); ++pn) {
        float x = R(pn, 0);
        float y = R(pn, 1);
        float z = R(pn, 2);
        float m = M(pn);

        // Compute assignment weights
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),
                                      Hy((y+0.5f)*nGrid),
                                      Hz((z+0.5f)*nGrid);

        // Accumulate contributions
        for(auto i=0; i < Order; ++i) {
            for(auto j=0; j < Order; ++j) {
                for(auto k=0; k < Order; ++k) {
                    #pragma omp atomic
                    grid(wrap(Hx.i + i), 
                         wrap(Hy.i + j), 
                         wrap(Hz.i + k))
                       += m * Hx.H[i] * Hy.H[j] * Hz.H[k];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{

    std::locale::global(std::locale("")); 
    std::cerr.imbue(std::locale());

    if (argc <= 1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size] [order]\n";
        return 1;
    }

    // Default grid size = 100
    int nGrid = 100;
    if (argc > 2) {
        nGrid = std::atoi(argv[2]);
    }

    // Default assignment order = 1
    int iOrder = 1;
    if (argc > 3) {
        iOrder = std::atoi(argv[3]);
    }

    auto t0 = hrc::now();

    std::ifstream io(argv[1], std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return 1;
    }

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h), sizeof(h))) {
        std::cerr << "Error reading header\n";
        return 1;
    }

    std::uint64_t N = h.nDark;  // number of dark matter particles
    std::cerr << "Loading " << N << " particles\n";

    // Allocate arrays for positions (Nx3) and masses (Nx1)
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // Load the particles
    tipsy::dark d;
    for(int i=0; i < N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
            std::cerr << "Error reading particle " << i << std::endl;
            return 1;
        }
        r(i,0) = d.pos[0];
        r(i,1) = d.pos[1];
        r(i,2) = d.pos[2];
        m(i)   = d.mass;
    }

    duration dt = hrc::now() - t0;
    std::cerr << "File reading took " << dt.count() << " seconds.\n";

    t0 = hrc::now();

    std::size_t totalSize = std::size_t(nGrid) 
                          * std::size_t(nGrid) 
                          * std::size_t(nGrid + 2);

    // Allocate float data, 64-byte aligned
    float* data = new (std::align_val_t(64)) float[ totalSize ];

    // Real Blitz++ array
    TinyVector<int,3> shapeReal(nGrid, nGrid, nGrid + 2);
    Array<float,3> padded(data, shapeReal, blitz::deleteDataWhenDone);

    // We will interpret the SAME memory as complex after the transform
    TinyVector<int,3> shapeCplx(nGrid, nGrid, (nGrid/2 + 1));
    Array<std::complex<float>, 3> kdata(
        reinterpret_cast<std::complex<float>*>(data),
        shapeCplx,
        blitz::neverDeleteData
    );

    // Zero the entire real array
    padded = 0.0f;

    // Create a subarray (nGrid, nGrid, nGrid) to skip the last 2 columns
    Array<float,3> subGrid = padded(Range::all(), Range::all(), Range(0, nGrid-1));

   
    std::cerr << "Assigning mass (Order=" << iOrder << ") ...\n";
    assign_mass<1>(subGrid, r, m); 
    // or: assign_mass<iOrder>(subGrid, r, m);

    dt = hrc::now() - t0;
    std::cerr << "Mass assignment took " << dt.count() << " seconds.\n";

  
    t0 = hrc::now();
    Array<float,2> projected(nGrid, nGrid);
    thirdIndex ii;
    projected = max(subGrid, ii);  
    // or max(padded, ii). They share the first nGrid columns.

    dt = hrc::now() - t0;
    std::cerr << "Projection took " << dt.count() << " seconds.\n";

    // Write out the 2D map (so you can plot and verify)
    std::ofstream of("density.dat", std::ios::binary);
    of.write(reinterpret_cast<char*>(projected.data()),
             projected.size() * sizeof(float));
    of.close();

   
    fftwf_plan plan = fftwf_plan_dft_r2c_3d(
        nGrid, nGrid, nGrid,
        reinterpret_cast<float*>(padded.data()),
        reinterpret_cast<fftwf_complex*>(kdata.data()),
        FFTW_ESTIMATE
    );

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);


    // Blitz++ frees 'data' automatically when 'padded' goes out of scope.
    std::cerr << "3D in-place FFT completed.\n";
    return 0;
}

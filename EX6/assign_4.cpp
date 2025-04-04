// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <new>
#include <complex>
#include <stdlib.h>
#include <vector>
#include <omp.h>         // For OpenMP
#include <limits>        // For std::numeric_limits
#include "blitz/array.h"
#include "fftw3.h"
#include "tipsy.h"
#include "aweights.h"

using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// A separate version is created for each different "Order".
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R, Array<float,1> &M) {
    auto nGrid = grid.rows();
    auto wrap = [nGrid](int i) {
        if (i<0) i += nGrid;
        else if (i>=nGrid) i -= nGrid;
        return i;
    };

    #pragma omp parallel for
    for(int pn=0; pn<R.rows(); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);
        AssignmentWeights<Order,float> 
            Hx((x+0.5f)*nGrid),
            Hy((y+0.5f)*nGrid),
            Hz((z+0.5f)*nGrid);
        for(int i=0; i<Order; ++i) {
            for(int j=0; j<Order; ++j) {
                for(int k=0; k<Order; ++k) {
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

int main(int argc, char *argv[]) {
    std::locale::global(std::locale(""));
    std::cerr.imbue(std::locale());

    if (argc<=2) {
        std::cerr << "Usage: " << argv[0] 
                  << " tipsyfile.std grid-size[/[L]bin-count] [order]"
                  << std::endl;
        return 1;
    }

    int nGrid = atoi(argv[2]);
    const auto iNyquist = nGrid / 2;
    int nBins = iNyquist;
    bool bLog = false;
    auto p = strchr(argv[2], '/');
    if (p) {
        *p++ = '\0';
        if (*p == 'L') {
            bLog = true;
            ++p;
        }
        nBins = atoi(p);
    }

    int iOrder = 1;
    if (argc > 3) iOrder = atoi(argv[3]);

    // 1) Read just the header (single-thread), then close
    auto t0 = hrc::now();
    tipsy::header h;
    {
        std::ifstream headerFile(argv[1], std::ifstream::binary);
        if (!headerFile) {
            std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
            return errno;
        }
        if (!headerFile.read(reinterpret_cast<char*>(&h), sizeof(h))) {
            std::cerr << "error reading header" << std::endl;
            return errno;
        }
        headerFile.close();
    }

    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles (parallel read)" << std::endl;

    // We assume N <= INT_MAX so we can cast to int safely.
    // If N might exceed INT_MAX, you need a different approach or 64-bit indexing library.
    if (N > std::numeric_limits<int>::max()) {
        std::cerr << "Error: N = " << N 
                  << " exceeds INT_MAX, cannot cast to int for blitz." 
                  << std::endl;
        return 1;
    }

    // Allocate arrays for positions and masses
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // 2) Parallel read of dark particles
    auto t0read = hrc::now();
    bool errorFlag = false;  // Will be set true if a thread fails opening

    #pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int tid      = omp_get_thread_num();

        std::uint64_t chunkSize = N / nthreads;
        std::uint64_t remainder = N % nthreads;

        std::uint64_t start = tid * chunkSize + std::min<std::uint64_t>(tid, remainder);
        std::uint64_t end   = (tid+1)*chunkSize + std::min<std::uint64_t>(tid+1, remainder);

        if (start >= end) {
            // No particles for this thread
            // Must NOT return here -> just skip
        }
        else {
            // Open separate ifstream
            std::ifstream ioThread(argv[1], std::ios::binary);
            if (!ioThread) {
                // Indicate error, but do not return inside parallel block
                #pragma omp critical
                {
                    std::cerr << "Thread " << tid << ": error opening file\n";
                    errorFlag = true;
                }
            }
            else {
                // Optional custom buffer
                static thread_local std::vector<char> localBuf(8*1024*1024);
                ioThread.rdbuf()->pubsetbuf(localBuf.data(), localBuf.size());

                // Seek to correct position:
                // offset = sizeof(header) + start*sizeof(tipsy::dark)
                std::streamoff offset = static_cast<std::streamoff>(
                    sizeof(tipsy::header) 
                    + start * sizeof(tipsy::dark)
                );
                ioThread.seekg(offset, std::ios::beg);

                tipsy::dark d;
                for (std::uint64_t i = start; i < end; i++) {
                    if (!ioThread.read(reinterpret_cast<char*>(&d), sizeof(d))) {
                        #pragma omp critical
                        {
                            std::cerr << "Thread " << tid 
                                      << ": error reading particle " << i 
                                      << std::endl;
                            errorFlag = true;
                        }
                        break;  // exit loop, but not the function
                    }

                    // Cast i to int so blitz::Array can do (int, int)
                    int idx = static_cast<int>(i);
                    r(idx,0) = d.pos[0];
                    r(idx,1) = d.pos[1];
                    r(idx,2) = d.pos[2];
                    m(idx)   = d.mass;
                }
                ioThread.close();
            } 
        } // end else
    } // end omp parallel

    if (errorFlag) {
        // If any thread had an error
        std::cerr << "Parallel read failed in one or more threads." << std::endl;
        return 1;
    }

    duration dtRead = hrc::now() - t0read;
    duration dt = hrc::now() - t0;
    double secRead  = dtRead.count();
    std::cerr << "Parallel file reading took " 
              << secRead << " seconds." << std::endl;
    double bytesRead = (sizeof(tipsy::header) + N*sizeof(tipsy::dark));
    double readRateMBs = bytesRead / secRead / (1024.0 * 1024.0);
    std::cerr << "File reading rate " 
              << readRateMBs << " MB/s." << std::endl;

    // = The rest of your code: mass assignment, FFT, etc. =
    auto k_nz     = nGrid / 2 + 1;
    auto n_floats = size_t(nGrid) * nGrid * (2*k_nz);
    float *data = new (std::align_val_t(64)) float[n_floats];
    Array<float,3> raw_grid(
        data,
        shape(nGrid,nGrid,2*k_nz),
        deleteDataWhenDone
    );
    Array<float,3> grid(
        raw_grid(Range(0,nGrid-1),
                 Range(0,nGrid-1),
                 Range(0,nGrid-1))
    );

    Array<std::complex<float>,3> kgrid(
        reinterpret_cast<std::complex<float>*>(data),
        shape(nGrid, nGrid, k_nz),
        neverDeleteData
    );

    // Assign the mass to the grid
    auto t0Mass = hrc::now();
    grid = 0;
    std::cerr << "Assigning mass to the grid using order " 
              << iOrder << std::endl;
    switch(iOrder) {
    case 1:
        assign_mass<1>(grid, r, m);
        break;
    case 2:
        assign_mass<2>(grid, r, m);
        break;
    case 3:
        assign_mass<3>(grid, r, m);
        break;
    case 4:
        assign_mass<4>(grid, r, m);
        break;
    default:
        std::cerr << "Invalid order " << iOrder 
                  << " (must be 1, 2, 3, or 4)" << std::endl;
    }
    duration dtMass = hrc::now() - t0Mass;
    std::cerr << "Mass assignment took " 
              << dtMass.count() << " seconds." << std::endl;
    std::cerr << "Total mass assigned is " 
              << blitz::sum(grid) << std::endl;

     // Calculate projected density
     t0 = hrc::now();
     Array<float,2> projected(nGrid,nGrid);
     thirdIndex ii;
     projected = max(grid,ii);
     dt = hrc::now() - t0;
     std::cerr << "Density projection took " << std::setw(9) << dt.count() << " seconds." << std::endl;
 
     // Write out the 2D map
     std::ofstream of("density_parallel.dat",std::ios::binary);
     of.write(reinterpret_cast<char*>(projected.data()),projected.size()*sizeof(float));
    // Calculate delta
    float total_mass = blitz::sum(grid);
    float diRhoBar   = (1.0f * nGrid * nGrid * nGrid) / total_mass;
    grid = grid * diRhoBar - 1.0f;
    grid /= (nGrid*nGrid*nGrid); // FFT normalization

    // Calculate FFT
    auto t0FFT = hrc::now();
    auto plan = fftwf_plan_dft_r2c_3d(
        nGrid, nGrid, nGrid,
        grid.dataFirst(),
        reinterpret_cast<fftwf_complex*>(kgrid.dataFirst()),
        FFTW_ESTIMATE
    );
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    duration dtFFT = hrc::now() - t0FFT;
    std::cerr << "FFT took " 
              << dtFFT.count() << " seconds." << std::endl;

    Array<double,1> ak(nBins), pk(nBins), nk(nBins);
    ak = 0; pk = 0; nk = 0;

    std::cerr << "Using " << nBins 
              << (bLog ? " logarithmic" : " linear") 
              << " bins." << std::endl;

    AssignmentWindow W(nGrid, iOrder);
    for (auto ii = kgrid.begin(); ii != kgrid.end(); ++ii) {
        auto pos = ii.position();
        auto bin = [iNyquist, nGrid](int k) {
            return (k <= iNyquist) ? k : k - nGrid;
        };
        auto kx = bin(pos[0]);
        auto ky = bin(pos[1]);
        auto kz = pos[2];
        float kval = std::sqrt(kx*kx + ky*ky + kz*kz);

        // Correction for mass assignment
        *ii *= W[std::abs(kx)] * W[std::abs(ky)] * W[kz];

        int idx;
        if (bLog) {
            idx = int(std::log(kval) / std::log(iNyquist) * nBins);
        } else {
            idx = int(kval / iNyquist * nBins);
        }
        if (idx > 0 && idx < nBins) {
            ak(idx) += kval;
            pk(idx) += std::norm(*ii);
            nk(idx) += 1;
        }
    }
    ak = ak / nk;
    pk = pk / nk;

    for(int i=1; i<nBins; ++i) {
        if (nk(i)) {
            printf("%.10g %.10g\n", ak(i), pk(i));
        }
    }

    return 0;
}

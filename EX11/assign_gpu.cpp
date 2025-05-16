// assign_gpu.cpp
//
// GPU version that:
//  1) Reads TIPSY data (stub).
//  2) Assigns mass onto a 3D grid (CPU).
//  3) Projects to 2D by taking max along z.
//  4) Writes the projected 2D array to "density_gpu.dat" (binary).
//  5) Performs a 2D FFT on the GPU (cuFFT).
//  6) Validates with a CPU inverse transform.
//
// Usage example:
//   ./assign_gpu tipsyfile.std 128 2
//
// Compile example:
//   nvcc -o assign_gpu assign_gpu.cpp -I. -lblitz -lfftw3f -lcufft -lcudart \
//        -O3 -std=c++17
// (Adjust paths/libraries for your system.)

#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdlib>
#include <new>

// Blitz
#include "blitz/array.h"

// FFTW (single-precision)
#include "fftw3.h"

// CUDA runtime + cuFFT
#include <cuda_runtime.h>
#include <cufft.h>

// Local stubs
#include "tipsy.h"
#include "aweights.h"

using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

//======================
// CUDA helper macros
//======================
#define CUDA_CHECK(err) do {                                  \
    cudaError_t e = (err);                                    \
    if (e != cudaSuccess) {                                   \
        std::cerr << "CUDA Error: " << cudaGetErrorString(e)  \
                  << " at line " << __LINE__ << std::endl;    \
        exit(EXIT_FAILURE);                                   \
    }                                                         \
} while(0)

#define CUFFT_CHECK(err) do {                                 \
    cufftResult r = (err);                                    \
    if (r != CUFFT_SUCCESS) {                                 \
        std::cerr << "CUFFT Error: " << r                     \
                  << " at line " << __LINE__ << std::endl;    \
        exit(EXIT_FAILURE);                                   \
    }                                                         \
} while(0)

//=========================
// CPU 2D inverse transform
//=========================
bool validate_2D(Array<float,2> &original, Array<std::complex<float>,2> &kdata)
{
    int nx = original.extent(0);
    int ny = original.extent(1);

    Array<float,2> recovered(nx, ny);

    fftwf_plan plan = fftwf_plan_dft_c2r_2d(nx, ny,
        reinterpret_cast<fftwf_complex*>(kdata.data()),
        recovered.data(), FFTW_ESTIMATE);

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Normalize
    recovered /= float(nx*ny);

    // Compare
    float tol = 1e-5f;
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            if (std::fabs(recovered(i,j) - original(i,j)) > tol){
                return false;
            }
        }
    }
    return true;
}

//=========================
// Assign mass function
//=========================
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R, Array<float,1> &M)
{
    int nGrid = grid.extent(0);
    auto wrap = [nGrid](int i){
        if (i < 0) i += nGrid;
        else if (i >= nGrid) i -= nGrid;
        return i;
    };

    #pragma omp parallel for
    for(int p=0; p < R.extent(0); p++){
        float x = R(p,0);
        float y = R(p,1);
        float z = R(p,2);
        float m = M(p);

        // "Center-of-cell" shift
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),
                                       Hy((y+0.5f)*nGrid),
                                       Hz((z+0.5f)*nGrid);

        for(int i=0; i<Order; i++){
            for(int j=0; j<Order; j++){
                for(int k=0; k<Order; k++){
                    #pragma omp atomic
                    grid(wrap(Hx.i + i), wrap(Hy.i + j), wrap(Hz.i + k))
                        += m * Hx.H[i] * Hy.H[j] * Hz.H[k];
                }
            }
        }
    }
}

//=========================
// Minimal TIPSY load
//=========================
bool load_tipsy_data(const char* fname, Array<float,2> &r, Array<float,1> &m, std::uint64_t &N)
{
    // Stub: fill random positions
    // Real code would read a TIPSY file "fname" and fill (r, m, N)
    N = 10;
    r.resize(N, 3);
    m.resize(N);

    for(int i=0; i<N; i++){
        r(i,0) = float(rand())/RAND_MAX;
        r(i,1) = float(rand())/RAND_MAX;
        r(i,2) = float(rand())/RAND_MAX;
        m(i)   = float(rand())/RAND_MAX;
    }
    std::cerr << "Stub: loaded " << N << " random particles from " << fname << "\n";
    return true;
}

//=========================
// Main
//=========================
int main(int argc, char* argv[])
{
    if (argc < 3){
        std::cerr << "Usage: " << argv[0] << " tipsyfile grid-size [order]\n";
        return 1;
    }

    int nGrid = std::atoi(argv[2]);
    int iOrder = 1;
    if (argc > 3) {
        iOrder = std::atoi(argv[3]);
    }

    // 1) Load tipsy data
    std::uint64_t N = 0;
    Array<float,2> r; // shape: Nx3
    Array<float,1> m; // shape: Nx1
    load_tipsy_data(argv[1], r, m, N);

    // 2) Create 3D grid (with padding in last dimension for possible R2C)
    int k_nz = nGrid / 2 + 1;
    size_t n_floats = (size_t)nGrid * nGrid * (2 * k_nz);
    float* data = new float[n_floats];

    // raw_grid is shape(nGrid, nGrid, 2*k_nz)
    Array<float,3> raw_grid(data, shape(nGrid, nGrid, 2*k_nz), neverDeleteData);

    // grid is the subregion ignoring last padding: shape(nGrid, nGrid, nGrid)
    Array<float,3> grid = raw_grid(Range(0, nGrid - 1),
                                   Range(0, nGrid - 1),
                                   Range(0, nGrid - 1));

    raw_grid = 0.0f;

    // 3) Assign mass
    std::cerr << "Assigning mass with order=" << iOrder << "...\n";
    switch (iOrder) {
        case 1: assign_mass<1>(grid, r, m); break;
        case 2: assign_mass<2>(grid, r, m); break;
        case 3: assign_mass<3>(grid, r, m); break;
        case 4: assign_mass<4>(grid, r, m); break;
        default: assign_mass<1>(grid, r, m); break;
    }

    // Output total mass in 3D
    float total_mass_3d = sum(grid);
    std::cerr << "Total mass assigned (3D) = " << total_mass_3d << "\n";

    // 4) Project to 2D via max over z
    Array<float,2> projected(nGrid, nGrid);
    projected = max(grid, thirdIndex());

    // Output sum of projected array
    float sum_proj = sum(projected);
    std::cerr << "Sum of 2D projection = " << sum_proj << "\n";

    // Write out the 2D map if desired
    {
        std::ofstream of("density_gpu.dat", std::ios::binary);
        of.write(reinterpret_cast<char*>(projected.data()),
                 projected.size() * sizeof(float));
        of.close();
        std::cerr << "Wrote 2D projection to density_gpu.dat (binary)\n";
    }

    // 5) GPU 2D FFT (in-place R2C)
    std::cerr << "Performing 2D GPU FFT...\n";
    int ny_padded = 2 * (nGrid / 2 + 1);
    Array<float,2> raw_data_2d(nGrid, ny_padded);
    raw_data_2d = 0.0f;

    // The "real" subregion (size nGrid x nGrid)
    Array<float,2> data_2d = raw_data_2d(Range(0, nGrid - 1),
                                        Range(0, nGrid - 1));
    data_2d = projected; // copy the 2D projection

    // Keep a copy for validation
    Array<float,2> data_2d_copy(data_2d.copy());

    // Allocate on GPU
    size_t size_bytes = raw_data_2d.size() * sizeof(float);
    float* d_rdata2d = nullptr;
    CUDA_CHECK(cudaMalloc(&d_rdata2d, size_bytes));
    CUDA_CHECK(cudaMemcpy(d_rdata2d, raw_data_2d.data(),
                          size_bytes, cudaMemcpyHostToDevice));

    // Create cuFFT plan for 2D R2C (in-place)
    cufftHandle plan;
    CUFFT_CHECK(cufftPlan2d(&plan, nGrid, nGrid, CUFFT_R2C));
    CUFFT_CHECK(cufftExecR2C(plan,
        reinterpret_cast<cufftReal*>(d_rdata2d),
        reinterpret_cast<cufftComplex*>(d_rdata2d)));
    CUFFT_CHECK(cufftDestroy(plan));

    // Copy data back to host
    CUDA_CHECK(cudaMemcpy(raw_data_2d.data(), d_rdata2d,
                          size_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_rdata2d));

    // Interpret result as complex
    Array<std::complex<float>,2> kdata(
        reinterpret_cast<std::complex<float>*>(raw_data_2d.data()),
        shape(nGrid, nGrid / 2 + 1),
        neverDeleteData
    );

    // 6) Validate with a CPU inverse transform
    bool ok = validate_2D(data_2d_copy, kdata);
    std::cerr << "GPU 2D FFT validation: " << (ok ? "MATCH" : "MISMATCH") << "\n";

    // Cleanup
    delete[] data;
    std::cerr << "Done.\n";
    return 0;
}

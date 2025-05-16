/******************************************************************************
 * Example: GPU 2D transform code that reads a TIPSY file, does mass assignment,
 * projects to 2D, does a cuFFT, and bins the result.
 *
 * Fixes:
 *   1) Cast 'i' to int in pos((int)i, 0) calls to avoid ambiguous overload.
 *   2) Fix missing parenthesis in CUDA_CHECK(...)
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <mpi.h>

// Blitz
#include "blitz/array.h"

// FFTW (single precision)
#include "fftw3.h"

// CUDA
#include <cuda_runtime.h>
#include <cufft.h>

// TIPSY + assignment
#include "tipsy.h"
#include "aweights.h"

using namespace blitz;
using hrc      = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;


#define CUDA_CHECK(err) do {                                  \
    cudaError_t e = (err);                                    \
    if (e != cudaSuccess) {                                   \
        std::cerr << "CUDA Error: " << cudaGetErrorString(e)  \
                  << " at line " << __LINE__ << std::endl;    \
        MPI_Abort(MPI_COMM_WORLD, 1);                         \
    }                                                         \
} while(0)

#define CUFFT_CHECK(err) do {                                 \
    cufftResult r = (err);                                    \
    if (r != CUFFT_SUCCESS) {                                 \
        std::cerr << "cuFFT Error: " << r                     \
                  << " at line " << __LINE__ << std::endl;    \
        MPI_Abort(MPI_COMM_WORLD, 1);                         \
    }                                                         \
} while(0)

// Minimal TIPSY loader: read all dark
bool load_tipsy_data(const char* fname,
                     Array<float,2> &pos,
                     Array<float,1> &mass,
                     std::uint64_t &N)
{
    std::ifstream in(fname, std::ios::binary);
    if(!in){
        std::cerr << "Unable to open TIPSY file " << fname << std::endl;
        return false;
    }
    tipsy::header h;
    if(!in.read(reinterpret_cast<char*>(&h), sizeof(h))){
        std::cerr << "Error reading TIPSY header\n";
        return false;
    }

    N = h.nDark;
    std::cerr << "Loading " << N << " dark particles from " << fname << "\n";
    pos.resize(N, 3);
    mass.resize(N);

    tipsy::dark d;
    for(std::uint64_t i=0; i<N; i++){
        if(!in.read(reinterpret_cast<char*>(&d), sizeof(d))){
            std::cerr << "Error reading TIPSY particle " << i << "\n";
            return false;
        }
        // cast 'i' to int to avoid ambiguous overload in Blitz
        pos((int)i, 0) = d.pos[0];
        pos((int)i, 1) = d.pos[1];
        pos((int)i, 2) = d.pos[2];
        mass((int)i)   = d.mass;
    }
    return true;
}

// Assign mass function

template<int Order=1>
void assign_mass_3D(Array<float,3> &grid,
                    Array<float,2> &R,
                    Array<float,1> &M)
{
    int nGrid = grid.extent(0);
    auto wrap = [nGrid](int i){
        if(i<0)      i += nGrid;
        else if(i>=nGrid) i -= nGrid;
        return i;
    };

    #pragma omp parallel for
    for(int p=0; p<R.extent(0); p++){  // p is int, so no overload issues
        float x = R(p,0);
        float y = R(p,1);
        float z = R(p,2);
        float m = M(p);

        AssignmentWeights<Order,float> 
            Hx((x+0.5f)*nGrid),
            Hy((y+0.5f)*nGrid),
            Hz((z+0.5f)*nGrid);

        for(int i=0; i<Order; i++){
            for(int j=0; j<Order; j++){
                for(int k=0; k<Order; k++){
                    #pragma omp atomic
                    grid(wrap(Hx.i + i), wrap(Hy.i + j), wrap(Hz.i + k))
                        += m * Hx.H[i]*Hy.H[j]*Hz.H[k];
                }
            }
        }
    }
}


// Validate 2D transform
bool validate_2D(Array<float,2> &original,
                 Array<std::complex<float>,2> &kdata)
{
    int nx= original.extent(0);
    int ny= original.extent(1);

    Array<float,2> recovered(nx, ny);
    fftwf_plan plan = fftwf_plan_dft_c2r_2d(nx, ny,
        reinterpret_cast<fftwf_complex*>(kdata.data()),
        recovered.data(), FFTW_ESTIMATE);

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    recovered /= float(nx*ny);

    float tol=1e-5f;
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            if(std::fabs(recovered(i,j)- original(i,j))>tol){
                return false;
            }
        }
    }
    return true;
}


// main
//=========================
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank=0, size=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // parse arguments...
    // (for brevity, not shown - just ensure you have nGrid >0, etc.)

    // 1) load TIPSY
    std::uint64_t N=0;
    Array<float,2> r; // Nx3
    Array<float,1> m; // Nx1
    if(!load_tipsy_data(argv[1], r, m, N)) {
        MPI_Finalize();
        return 1;
    }


    // 2) create 3D grid
    int nGrid= 128;
    int k_nz= nGrid/2 +1;
    size_t n_floats= size_t(nGrid)* nGrid* (2*k_nz);
    float* data= new float[n_floats];

    // shape(nGrid,nGrid,2*k_nz)
    Array<float,3> raw_grid(data, shape(nGrid,nGrid,2*k_nz), neverDeleteData);
    Array<float,3> grid= raw_grid(Range(0,nGrid-1),
                                  Range(0,nGrid-1),
                                  Range(0,nGrid-1));
    raw_grid=0.f;

    // assign mass
    assign_mass_3D<1>(grid, r, m); // e.g., order=1

    // project => 2D
    Array<float,2> projected(nGrid,nGrid);
    {
        thirdIndex kk;
        projected= max(grid, kk);
    }

    // form delta2D => (rho/mean -1)/(n^2)
    float total2D= sum(projected);
    float mean2D= total2D/(nGrid*nGrid);
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            float val= projected(i,j)/mean2D -1.f;
            // typical normalization for r2c
            projected(i,j)= val/(nGrid*nGrid);
        }
    }

    // 3) GPU 2D in-place R2C
    int ny_padded= 2*(nGrid/2+1);
    Array<float,2> raw_data2d(nGrid, ny_padded);
    raw_data2d=0.f;

    Array<float,2> real_2d= raw_data2d(Range(0,nGrid-1),
                                       Range(0,nGrid-1));
    real_2d= projected;

    Array<float,2> real_2d_copy(real_2d.copy());

    size_t size_bytes= raw_data2d.size()*sizeof(float);
    float* d_rdata= nullptr;
    CUDA_CHECK(cudaMalloc(&d_rdata, size_bytes));
    CUDA_CHECK(cudaMemcpy(d_rdata, raw_data2d.data(), 
                          size_bytes, cudaMemcpyHostToDevice));

    cufftHandle plan2d;
    CUFFT_CHECK(cufftPlan2d(&plan2d, nGrid, nGrid, CUFFT_R2C));
    CUFFT_CHECK(cufftExecR2C(plan2d,
        reinterpret_cast<cufftReal*>(d_rdata),
        reinterpret_cast<cufftComplex*>(d_rdata)));
    CUFFT_CHECK(cufftDestroy(plan2d));

    CUDA_CHECK(cudaMemcpy(raw_data2d.data(), d_rdata, 
                          size_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_rdata));

    // interpret as shape(nGrid, nGrid/2+1)
    Array<std::complex<float>,2> kdata(
        reinterpret_cast<std::complex<float>*>(raw_data2d.data()),
        shape(nGrid, nGrid/2+1),
        neverDeleteData
    );

    // 4) validate
    bool ok= validate_2D(real_2d_copy, kdata);
    std::cerr<<"GPU 2D FFT validation: "<<(ok?"MATCH":"MISMATCH")<<"\n";

    delete[] data;
    MPI_Finalize();
    return 0;
}

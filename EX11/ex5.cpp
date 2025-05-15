#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <complex>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <iomanip>

// For GPU
#include <cuda_runtime.h>
#include <cufft.h>

// For CPU FFT
#include <fftw3.h>

// For arrays and mass assignment
#include "blitz/array.h"
#include "aweights.h"     // Must contain AssignmentWeights template
#include "tipsy.h"        // Must contain tipsy::header, tipsy::dark structs

using namespace blitz;
using std::complex;

static void fill_array(Array<float,2> &data, float value=0.0f) {
    data = value;
}

//--------------------------------------------------------------------
// Simple function: read x,y,z from tipsy file, store x,y in a 2D array
// Also store mass in 1D array. We do no parallel/MPI here:
// just read the entire file in rank=0 style. 
//--------------------------------------------------------------------
void read_tipsy_2d(const char* filename,
                   std::vector<float> &vx,
                   std::vector<float> &vy,
                   std::vector<float> &vz, // we won't use z, but let's read it
                   std::vector<float> &vm)
{
    std::ifstream io(filename, std::ifstream::binary);
    if(!io) {
        std::cerr << "Cannot open tipsy file " << filename << std::endl;
        exit(1);
    }
    tipsy::header h;
    if(!io.read(reinterpret_cast<char*>(&h), sizeof(h))) {
        std::cerr << "Error reading tipsy header\n";
        exit(1);
    }

    std::uint64_t N = h.nDark; 
    vx.resize(N); 
    vy.resize(N);
    vz.resize(N);
    vm.resize(N);

    // read dark matter particles
    tipsy::dark d;
    for(uint64_t i=0; i<N; i++) {
        if(!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
            std::cerr << "Error reading dark matter particle\n";
            exit(1);
        }
        vx[i] = d.pos[0];
        vy[i] = d.pos[1];
        vz[i] = d.pos[2];
        vm[i] = d.mass;
    }
}

//--------------------------------------------------------------------
// Example mass assignment in 2D: we ignore z, using x,y only
// We deposit onto an nGrid x nGrid array
//--------------------------------------------------------------------
template<int Order>
void assign_mass_2d(Array<float,2> &grid,
                    const std::vector<float> &vx,
                    const std::vector<float> &vy,
                    const std::vector<float> &vm)
{
    int nGrid = grid.rows(); // We assume square grid
    auto wrap = [nGrid](int i) {
        if(i<0) i += nGrid;
        else if(i>=nGrid) i -= nGrid;
        return i;
    };

    #pragma omp parallel for
    for(size_t p=0; p<vx.size(); p++) {
        float x = vx[p];
        float y = vy[p];
        float m = vm[p];

        // We shift by +0.5 so that we deposit 
        // in a cell-based approach (like the 3D code).
        AssignmentWeights<Order,float> Hx((x + 0.5f)*nGrid);
        AssignmentWeights<Order,float> Hy((y + 0.5f)*nGrid);

        for(int i=0; i<Order; i++) {
            for(int j=0; j<Order; j++) {
                #pragma omp atomic
                grid(wrap(Hx.i + i), wrap(Hy.i + j)) += m * Hx.H[i] * Hy.H[j];
            }
        }
    }
}


//--------------------------------------------------------------------
// For debugging, let's do a CPU inverse transform with FFTW
// and compare to the original real array, to confirm 
// forward transform is correct
//--------------------------------------------------------------------
bool validate(Array<float,2> &orig, Array<std::complex<float>,2> &kdata)
{
    // We'll do c2r with FFTW in 2D
    int nx = orig.rows();
    int ny = orig.cols();
    Array<float,2> rtemp(nx, ny);

    fftwf_plan plan = fftwf_plan_dft_c2r_2d(nx, ny,
        reinterpret_cast<fftwf_complex*>(kdata.data()),
        rtemp.data(), FFTW_ESTIMATE);

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Because forward transform is unscaled, the inverse 
    // returns sum-of data * (nx * ny). So we divide:
    rtemp /= float(nx*ny);

    // Check difference
    float tol = 1e-5f;
    // Use Blitz: all(abs(orig - rtemp) < 1e-5)
    Array<float,2> diff(orig - rtemp);
    float maxdiff = max(abs(diff));
    // For a large NxN, some floating differences might exceed 1e-5 occasionally.
    // But let's do it anyway:
    bool matched = all(abs(diff) < tol);
    std::cerr << "max difference = " << maxdiff << std::endl;
    return matched;
}


//====================================================================
// MAIN
//====================================================================
int main(int argc, char** argv)
{
    if(argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " tipsyfile.std  nGrid  [order=1..4]\n";
        return 1;
    }
    const char* tipsyfile = argv[1];
    int nGrid = atoi(argv[2]);
    int iOrder = (argc > 3 ? atoi(argv[3]) : 1);
    if(iOrder <1 || iOrder>4) {
        std::cerr << "Invalid order " << iOrder << " (must be 1..4)\n";
        return 1;
    }

    // 1) Read tipsy data into std::vectors
    std::vector<float> vx, vy, vz, vm;
    read_tipsy_2d(tipsyfile, vx, vy, vz, vm);
    std::cerr << "Loaded " << vx.size() << " particles.\n";

    // 2) Construct a 2D real grid of shape (nGrid x nGrid) 
    //    to deposit mass
    Array<float,2> real_grid(nGrid, nGrid);
    fill_array(real_grid, 0.0f);

    // 3) Assign mass with the specified order (1..4).
    std::cerr << "Assigning mass with order " << iOrder << std::endl;
    switch(iOrder) {
        case 1: assign_mass_2d<1>(real_grid, vx, vy, vm); break;
        case 2: assign_mass_2d<2>(real_grid, vx, vy, vm); break;
        case 3: assign_mass_2d<3>(real_grid, vx, vy, vm); break;
        case 4: assign_mass_2d<4>(real_grid, vx, vy, vm); break;
    }

    // 4) Convert to overdensity if you want. 
    //    Typically, sum_of_grid = total mass. Then we do delta = rho/rhoBar - 1
    //    For demonstration, let's just keep it as-is 
    //    (or do something akin to your 3D code).
    float total_mass=0.f;
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            total_mass += real_grid(i,j);
        }
    }
    std::cerr << "Total mass on grid = " << total_mass << std::endl;
    // Example normalization:
    float inv_mean = (nGrid*nGrid)/ (std::max(total_mass,1e-20f));
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            real_grid(i,j) = real_grid(i,j)*inv_mean - 1.f;
        }
    }

    // 5) Prepare for in-place R2C transform on GPU
    //    We need a padded array: shape (n, 2*(n/2 +1)) in float.
    int n_padded = 2*(nGrid/2 + 1);
    Array<float,2> raw_data2d(nGrid, n_padded);
    // "view" the left portion as (nGrid x nGrid) real data
    Array<float,2> data2d = raw_data2d(Range(0,nGrid-1), Range(0,nGrid-1));

    // Copy from "real_grid" into data2d
    data2d = real_grid;

    // Save a copy for validation
    Array<float,2> original_copy(data2d.copy());

    // 6) Allocate GPU memory for the entire padded array
    size_t size_in_bytes = sizeof(float)*raw_data2d.size();
    float* d_data;
    cudaMalloc(&d_data, size_in_bytes);

    // Copy host -> device
    cudaMemcpy(d_data, raw_data2d.data(), size_in_bytes, cudaMemcpyHostToDevice);

    // 7) Create cuFFT plan for an in-place R2C 2D transform
    cufftHandle plan;
    cufftPlan2d(&plan, nGrid, nGrid, CUFFT_R2C);

    // 8) Execute in-place: input= d_data, output= d_data
    cufftExecR2C(plan,
                 reinterpret_cast<cufftReal*>(d_data),
                 reinterpret_cast<cufftComplex*>(d_data));
    cufftDestroy(plan);

    // 9) Copy device -> host
    cudaMemcpy(raw_data2d.data(), d_data, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaFree(d_data);

    // 10) Now interpret raw_data2d's left portion as complex
    // Because R2C in-place puts complex data in shape (nGrid, nGrid/2+1)
    // We can use:
    Array<std::complex<float>,2> kdata2d(
        reinterpret_cast<std::complex<float>*>(data2d.data()),
        shape(nGrid, nGrid/2 + 1),
        neverDeleteData
    );

    // 11) Validate by performing CPU inverse FFT with fftwf_plan_dft_c2r_2d
    bool good = validate(original_copy, kdata2d);
    std::cout << "GPU in-place FFT " << (good?"match":"MISMATCH") << std::endl;

    return 0;
}

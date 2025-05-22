#include <blitz/array.h>
#include "cufft.h"
#include <cuda_runtime.h>

cufftHandle gpu_make_plan_2D(int nGrid);
void gpu_fft_2D_R2C(blitz::Array<float,2> &grid,void *slab,cufftHandle plan, cudaStream_t stream);
void *gpu_allocate_slab(size_t nGrid);

cufftHandle gpu_make_plan_1D(int nGrid);
void gpu_fft_1D_C2C(blitz::Array<std::complex<float>,2> &grid, void *slab,
                    cufftHandle plan, cudaStream_t stream);
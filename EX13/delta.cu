#include "delta.h"
#include <stdio.h>


//EX3
__global__
void delta_kernel(float *slab,
                  int    nGrid,         
                  int    embed,          
                  float  inv_rhoBar,     
                  float  inv_nGrid3)    
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;  
    const int j = blockIdx.y * blockDim.y + threadIdx.y;  
    if (i >= nGrid || j >= nGrid) return;
    const int stride_k = embed * nGrid;        
    int idx = i + j * embed;                   
    for (int k = 0; k < nGrid; ++k, idx += stride_k)
    {
        float v = slab[idx];                  
        v = v * inv_rhoBar - 1.0f;            
        slab[idx] = v * inv_nGrid3;            
    }
}
//EX2
// __global__
// void delta_kernel(float *data,
//                   int    nGrid,          // cube side length
//                   float  inv_rhoBar,     // 1 / ρ̄
//                   float  inv_nGrid3)     // 1 / (nGrid³)
// {
//     /* 2-D thread coordinates in the x-y plane */
//     int ix = blockIdx.x * blockDim.x + threadIdx.x;   // 0 … nGrid-1
//     int iy = blockIdx.y * blockDim.y + threadIdx.y;   // 0 … nGrid-1
//     if (ix >= nGrid || iy >= nGrid) return;           // outside the slab

//     /* Each thread loops over z for its (x,y) column */
//     for (int iz = 0; iz < nGrid; ++iz) {
//         int idx = ix * nGrid * nGrid + iy * nGrid + iz;   // flatten (x,y,z)
//         float val = data[idx];
//         val  = val * inv_rhoBar - 1.0f;
//         data[idx] = val * inv_nGrid3;
//     }
// }


// void compute_delta(float       *d_slab,
//                    int          nGrid,
//                    float        rhoBar,
//                    cudaStream_t stream)
// {
//     const float inv_rho    = 1.0f / rhoBar;
//     const float inv_nGrid3 = 1.0f / float(nGrid * nGrid * nGrid);

//     /* 32×32 = 1024-thread blocks (fits CUDA limit) */
//     constexpr int blockXY = 32;
//     dim3 dimBlock(blockXY, blockXY, 1);

//     /* number of blocks per dimension, rounded up */
//     int nBlocks = (nGrid + blockXY - 1) / blockXY;
//     dim3 dimGrid(nBlocks, nBlocks, 1);

//     delta_kernel<<<dimGrid, dimBlock, 0, stream>>>(d_slab,
//                                                    nGrid,
//                                                    inv_rho,
//                                                    inv_nGrid3);
// }



//EX1
// __global__
// void delta_kernel(float *data, int n, float inv_rhoBar, float inv_nGrid3) {
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx < n) {
//         float val = data[idx];
//         val = val * inv_rhoBar - 1.0f;
//         data[idx] = val * inv_nGrid3;
//     }
// }

// void compute_delta(float *d_slab, int nGrid, float rhoBar, cudaStream_t stream) {
//     int n = nGrid * nGrid * nGrid;
//     float inv_rho = 1.0f / rhoBar;
//     float inv_n = 1.0f / float(n);

//     int threads = 256;
//     int blocks = (n + threads - 1) / threads;
//     delta_kernel<<<blocks, threads, 0, stream>>>(d_slab, n, inv_rho, inv_n);
// }
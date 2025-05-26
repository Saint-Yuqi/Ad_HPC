#include "delta.h"
#include <stdio.h>

__global__
void delta_kernel(float *data, int n, float inv_rhoBar, float inv_nGrid3) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        float val = data[idx];
        val = val * inv_rhoBar - 1.0f;
        data[idx] = val * inv_nGrid3;
    }
}

void compute_delta(float *d_slab, int nGrid, float rhoBar, cudaStream_t stream) {
    int n = nGrid * nGrid * nGrid;
    float inv_rho = 1.0f / rhoBar;
    float inv_n = 1.0f / float(n);

    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    delta_kernel<<<blocks, threads, 0, stream>>>(d_slab, n, inv_rho, inv_n);
}

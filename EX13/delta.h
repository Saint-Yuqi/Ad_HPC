#ifndef DELTA_H
#define DELTA_H
#include <cuda_runtime.h>

void compute_delta(float *d_slab, int nGrid, float rhoBar, cudaStream_t stream);

#endif // DELTA_H
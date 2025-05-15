#include <cuda_runtime.h>
#include <cufft.h>
#include <cstdio>
#include <cstdlib>

int main() {
    int N = 1024;  // example size
    // total real size = N * 2*(N/2+1)
    size_t realSize = sizeof(float) * N * 2*(N/2 + 1);

    // Allocate device memory for the padded real array:
    float* d_data;
    cudaMalloc(&d_data, realSize);

    // (You would copy your host data to d_data here, if needed)
    //   cudaMemcpy(d_data, hostData, realSize, cudaMemcpyHostToDevice);

    // Plan many parameters:
    int n[2]       = { N, N };
    int inembed[2] = { N, 2*(N/2 + 1) };
    int onembed[2] = { N, (N/2 + 1) };

    int istride = 1;
    int ostride = 1;
    int idist   = N * 2*(N/2 + 1);
    int odist   = N * (N/2 + 1);
    int batch   = 1;  // single 2D FFT

    cufftHandle plan;
    cufftResult r = cufftPlanMany(&plan,
        2, n,               // rank=2, transform size
        inembed, istride, idist,
        onembed, ostride, odist,
        CUFFT_R2C, batch);
    if (r != CUFFT_SUCCESS) {
        printf("cufftPlanMany failed, error=%d\n", r);
        return 1;
    }

    // Execute the in-place transform
    r = cufftExecR2C(plan,
        reinterpret_cast<cufftReal*>(d_data),
        reinterpret_cast<cufftComplex*>(d_data));
    if (r != CUFFT_SUCCESS) {
        printf("cufftExecR2C failed, error=%d\n", r);
        return 1;
    }

    cufftDestroy(plan);

    // Now d_data has the complex output in the first N*(N/2+1) complexes
    // (plus some leftover in the padded region).
    // You could copy it back or do more processing on GPU.

    cudaFree(d_data);
    return 0;
}

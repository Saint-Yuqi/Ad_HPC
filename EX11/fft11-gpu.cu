#include <blitz/array.h>
#include <fftw3.h>
#include <complex>
#include <cmath>
#include <iostream>
#include <cuda.h>
#include <cufft.h>
using namespace blitz;
using std::complex;

void fill_array(Array<float,2> &data) {
    // Set the grid to the sum of two sine functions
    for (int i=0; i < data.rows(); i++) {
        for (int j=0; j < data.cols(); j++) {
            float x = (float)i / 25.0; // Period of 1/4 of the box in x
            float y = (float)j / 10.0; // Period of 1/10 of the box in y
            data(i,j) = sin(2.0 * M_PI * x) + sin(2.0 * M_PI * y);
        }
    }
}

// Verify the FFT (kdata) of data by performing a reverse transform and comparing
bool validate(Array<float,2> &data,Array<std::complex<float>, 2> kdata) {
    Array<float,2> rdata(data.extent());
    fftwf_plan plan = fftwf_plan_dft_c2r_2d(data.rows(), data.cols(),
        reinterpret_cast<fftwf_complex*>(kdata.data()), rdata.data(), FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    rdata /= data.size(); // Normalize for the FFT
    return all(abs(data - rdata) < 1e-5);
}

int main() {
    int n = 10000;

    // Out of place: rdata1 -> kdata1 without padding
    Array<float,2> rdata1(n,n);
    Array<std::complex<float>, 2> kdata1(n, n/2 + 1);
    fftwf_plan plan1  = fftwf_plan_dft_r2c_2d(n, n,
        rdata1.data(), reinterpret_cast<fftwf_complex*>(kdata1.data()), FFTW_ESTIMATE);
    fill_array(rdata1);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    std::cout << ">>> Out of place FFT " << (validate(rdata1,kdata1)?"match":"MISMATCH") << endl;

    // in-place: rdata2 -> kdata2 (overlays rdata2)
    Array<float,2> raw_data2(n,2*(n/2+1));
    Array<float,2> rdata2 = raw_data2(Range(0,n-1),Range(0,n-1));
    fftwf_plan plan2  = fftwf_plan_dft_r2c_2d(n, n,
        rdata2.data(), reinterpret_cast<fftwf_complex*>(rdata2.data()), FFTW_ESTIMATE);
    fill_array(rdata2);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);
    Array<std::complex<float>, 2> kdata2(reinterpret_cast<std::complex<float>*>(rdata2.data()),
        shape(n, n/2 + 1),neverDeleteData);
    // rdata2 is destroyed now, but it would have the same data as rdata1 so compare with rdata1
    std::cout << ">>> In-place FFT " << (validate(rdata1,kdata2)?"match":"MISMATCH") << endl;
    
    //Ex3
    Array<float,2> raw_data3(n, 2*(n/2 + 1));
    Array<float,2> data3 = raw_data3(Range(0,n-1), Range(0,n-1));
    fill_array(data3);

    // Create a second raw_data4/data4 array which remains empty initially
    Array<float,2> raw_data4(n, 2*(n/2 + 1));
    Array<float,2> data4 = raw_data4(Range(0,n-1), Range(0,n-1));

    size_t size_in_bytes = raw_data3.size() * sizeof(float);
    void *device_data;
    cudaMalloc(&device_data, size_in_bytes);

    // Transfer raw_data3 to the device (HostToDevice)
    cudaMemcpy(device_data, raw_data3.data(), size_in_bytes, cudaMemcpyHostToDevice);

    // Transfer data back into raw_data4 (DeviceToHost)
    cudaMemcpy(raw_data4.data(), device_data, size_in_bytes, cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(device_data);
    bool gpu_copy_ok = all(abs(data3 - data4) < 1e-5);
    std::cout << ">>> GPU copy test " << (gpu_copy_ok ? "match" : "MISMATCH") << std::endl;
    
    //EX4
    Array<float,2> raw_data5(n, 2 * (n/2 + 1));
    Array<float,2> data5 = raw_data5(Range(0,n-1),Range(0,n-1));
    fill_array(data5);

    Array<float, 2> data5_copy(data5.copy());

    float* rdata5;
    size_t size_bytes = raw_data5.size() * sizeof(float);

    cudaMalloc(&rdata5, size_bytes);
    cudaMemcpy(rdata5, raw_data5.data(), size_bytes, cudaMemcpyHostToDevice);

    cufftHandle plan5;
    cufftPlan2d(&plan5, n, n, CUFFT_R2C);

    cufftExecR2C(plan5,
        reinterpret_cast<cufftReal*>(rdata5),
        reinterpret_cast<cufftComplex*>(rdata5));
    cufftDestroy(plan5);

    cudaMemcpy(raw_data5.data(), rdata5, size_bytes, cudaMemcpyDeviceToHost);
    cudaFree(rdata5);

    Array<std::complex<float>, 2> kdata5(
        reinterpret_cast<std::complex<float>*>(raw_data5.data()),
        shape(n, n/2 + 1),
        neverDeleteData
    );

    std::cout << ">>> GPU in-place FFT "
          << (validate(data5_copy, kdata5) ? "match" : "MISMATCH")
          << std::endl;
 
    return 0;

}

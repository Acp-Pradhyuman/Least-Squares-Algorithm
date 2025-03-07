#include <iostream>
#include <omp.h>
#include <cmath>
#include <cuda_runtime.h>
using namespace std;

const double LsEpsilon = 1.0e-12;

// CUDA kernel for calculating sums for least squares
__global__ void calcLeastSquaresKernel(const double* x, const double* y, int n, 
    double* sumX, double* sumY, double* sumXY, double* sumXX) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        atomicAdd(sumX, x[idx]);
        atomicAdd(sumY, y[idx]);
        atomicAdd(sumXX, x[idx] * x[idx]);
        atomicAdd(sumXY, x[idx] * y[idx]);
    }
}

// y = mx + b (Least Squares Calculation)
bool calcLeastSquaresCUDA(const double* x, const double* y, int n, 
    double* m, double* b) {
    if (n <= 0) {
        return false;
    }

    double *d_x, *d_y, *d_sumX, *d_sumY, *d_sumXY, *d_sumXX;
    double h_sumX = 0, h_sumY = 0, h_sumXY = 0, h_sumXX = 0;
    
    cudaMalloc(&d_x, n * sizeof(double));
    cudaMalloc(&d_y, n * sizeof(double));
    cudaMalloc(&d_sumX, sizeof(double));
    cudaMalloc(&d_sumY, sizeof(double));
    cudaMalloc(&d_sumXY, sizeof(double));
    cudaMalloc(&d_sumXX, sizeof(double));

    cudaMemcpy(d_x, x, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sumX, &h_sumX, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sumY, &h_sumY, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sumXY, &h_sumXY, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sumXX, &h_sumXX, sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int gridSize = (n + blockSize - 1) / blockSize;

    // Start time before launching the kernel
    // double start_time = omp_get_wtime();
    // End time
    double start_time = static_cast<double>(clock()) / CLOCKS_PER_SEC;

    // Launch kernel to calculate sums in parallel on the GPU
    calcLeastSquaresKernel<<<gridSize, blockSize>>>(d_x, d_y, n, 
        d_sumX, d_sumY, d_sumXY, d_sumXX);

    // Wait for GPU to finish (ensure synchronization before measuring time)
    cudaDeviceSynchronize();

    // End time after synchronization
    // double end_time = omp_get_wtime();
    // End time
    double end_time = static_cast<double>(clock()) / CLOCKS_PER_SEC;

    double time_taken = end_time - start_time;
    printf("Time: %f seconds\n", time_taken);

    // Copy the results back to host
    cudaMemcpy(&h_sumX, d_sumX, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_sumY, d_sumY, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_sumXY, d_sumXY, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_sumXX, d_sumXX, sizeof(double), cudaMemcpyDeviceToHost);

    double denom = n * h_sumXX - h_sumX * h_sumX;
    if (LsEpsilon >= fabs(denom)) {
        return false;
    }

    *m = (n * h_sumXY - h_sumX * h_sumY) / denom;
    *b = (h_sumXX * h_sumY - h_sumX * h_sumXY) / denom;

    // Free device memory
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_sumX);
    cudaFree(d_sumY);
    cudaFree(d_sumXY);
    cudaFree(d_sumXX);

    return true;
}

const int n = 5000 * 5000;
double x[n];
double y[n];

int main() {
    const double slope = 1.0;
    const double y_int = 0.5;

    // set to 8 threads (since 8 threads gives better parallelization)
    int num_threads = 8;
    omp_set_num_threads(num_threads);
        
    // Initialize the input vectors
    // All points lie on the line y = 1*x + 0.5
    #pragma omp parallel for
    for (int i = 0; i < n; i += 4) {
        x[i] = i;
        y[i] = slope * x[i] + y_int;
        x[i+1] = (i+1);
        y[i+1] = slope * x[i+1] + y_int;
        x[i+2] = (i+2);
        y[i+2] = slope * x[i+2] + y_int;
        x[i+3] = (i+3);
        y[i+3] = slope * x[i+3] + y_int;        
    }

    double m1 = 0, b1 = 0;
    bool rv1;

    // Use CUDA to calculate least squares
    rv1 = calcLeastSquaresCUDA(x, y, n, &m1, &b1);

    printf("slope = %6.4lf, m1 = %6.4lf\n", slope, m1);
    cout << "y_int = " << y_int << ", b1 = " << b1 << endl;
    cout << "rv1 = " << rv1 << endl;

    return 0;
}
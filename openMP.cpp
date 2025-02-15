#include <iostream>
#include <omp.h>
#include <cmath>
using namespace std;

const double LsEpsilon = 1.0e-12;

//y = mx + b
bool calcLeastSquaresCPP(const double* x, const double* y, 
    int n, double* m, double* b);

const int n = 5000 * 5000;
// alignas(64) double x[n];
// alignas(64) double y[n];
double x[n];
double y[n];

int main()
{
    const double slope = 1.0;
    const double y_int = 0.5;
    double start_time, end_time;
    double time_with_1_thread;

    // Loop over different thread counts
    for (int num_threads = 1; num_threads <= 64; num_threads *= 2) {
        omp_set_num_threads(num_threads);
        // Initialize the input vectors
        // All points lie on the line y = 1*x + 0.5
        start_time = omp_get_wtime();

        #pragma omp parallel for
        for (int i = 0; i < n; i+=4)
        {
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

        rv1 = calcLeastSquaresCPP(x, y, n, &m1, &b1);

        end_time = omp_get_wtime();

        // Capture the time with 1 thread for calculating speedup
        if (num_threads == 1) {
            time_with_1_thread = end_time - start_time;
        }

        double time_taken = end_time - start_time;
        double speedup = time_with_1_thread / time_taken;

        printf("slope = %6.4lf, m1 = %6.4lf\n", slope, m1);
        cout << "y_int = " << y_int << ", b1 = " << b1 << endl;
        cout << "rv1 = " << rv1 << endl;

        printf("Threads: %d, Time: %f seconds, Speedup: %f\n", 
            num_threads, time_taken, speedup);
    }
}

bool calcLeastSquaresCPP(const double* x, const double* y, int n, 
    double* m, double* b) {
    if (n <= 0) {
        return false;
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    #pragma omp parallel for reduction(+:sumX, sumY, sumXY, sumXX)
    for(int i=0; i<n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXX += x[i] * x[i];
        sumXY += x[i] * y[i];
    }

    double denom = n * sumXX - sumX * sumX;
    if (LsEpsilon >= fabs(denom)) return false;

    *m = (n * sumXY - sumX * sumY) / denom;
    *b = (sumXX*sumY - sumX*sumXY) / denom;

    return true;
}
#include <iostream>
#include <cmath>
#include <ctime>  // For time measurement

using namespace std;

const double LsEpsilon = 1.0e-12;

// y = mx + b (Least Squares Calculation)
bool calcLeastSquaresCPP(const double* x, const double* y, int n, 
    double* m, double* b) {
    if (n <= 0) {
        return false;
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

    double start_time, end_time;

    // Start time
    start_time = static_cast<double>(clock()) / CLOCKS_PER_SEC;
    
    for(int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXX += x[i] * x[i];
        sumXY += x[i] * y[i];
    }

    // End time
    end_time = static_cast<double>(clock()) / CLOCKS_PER_SEC;

    // Calculate elapsed time
    double time_taken = end_time - start_time;

    // Print execution time
    printf("Time: %f seconds\n", time_taken);

    double denom = n * sumXX - sumX * sumX;
    if (LsEpsilon >= fabs(denom)) return false;

    *m = (n * sumXY - sumX * sumY) / denom;
    *b = (sumXX * sumY - sumX * sumXY) / denom;

    return true;
}

const int n = 10000 * 10000; // Size of data
double x[n];
double y[n];

int main() {
    const double slope = 1.0;
    const double y_int = 0.5;

    // Initialize the input vectors (serial initialization)
    for (int i = 0; i < n; i++)
    {
        x[i] = i;
        y[i] = slope * x[i] + y_int;		
    }

    double m1 = 0, b1 = 0;
    bool rv1;

    // Call the least squares calculation function
    rv1 = calcLeastSquaresCPP(x, y, n, &m1, &b1);

    // Output results
    printf("slope = %6.4lf, m1 = %6.4lf\n", slope, m1);
    cout << "y_int = " << y_int << ", b1 = " << b1 << endl;
    cout << "rv1 = " << rv1 << endl;

    return 0;
}
// leastSquaresAlgo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
using namespace std;

extern "C" double LsEpsilon;

//y = mx + b
extern "C" bool calcLeastSquaresASM(const double* x, const double* y, int n, double* m, 
    double *b);
bool calcLeastSquaresCPP(const double* x, const double* y, int n, double* m, double* b);

int main()
{
    //std::cout << "Hello World!\n";
    const int n = 6;
    double x[n] = {0, 2, 4, 6, 8, 10};
    double y[n] = {51.23, 34.6, 12.3, 56.8, 90.1, 111.9};

    double m1 = 0, m2 = 0, b1 = 0, b2 = 0;
    bool rv1, rv2;

    rv1 = calcLeastSquaresASM(x, y, n, &m1, &b1);
    rv2 = calcLeastSquaresCPP(x, y, n, &m2, &b2);

    for(int i=0; i<n; i++) {
        printf("x[%d] = %12.4lf, y[%d] = %12.4lf\n", i, x[i], i, y[i]);
    }

    printf("m1 = %12.4lf, m2 = %12.4lf\n", m1, m2);
    cout << "b1 = " << b1 << ", b2 = " << b2 << endl;
    cout << "rv1 = " << rv1 << ", rv2 = " << rv2 << endl;
    // std::cout << "m1 = " << m1 << ", m2 = " << m2 << std::endl;
    // cout << "b1 = " << b1 << ", b2 = " << b2 << endl;
    // cout << "rv1 = " << rv1 << ", rv2 = " << rv2 << endl;
}

bool calcLeastSquaresCPP(const double* x, const double* y, int n, double* m, double* b) {
    if (n <= 0) {
        return false;
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
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

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

// leastSquaresAlgo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono> 
using namespace std;

extern "C" double LsEpsilon;

//y = mx + b
extern "C" bool calcLeastSquaresASM(const double* x, const double* y, int n, double* m, 
    double *b);
bool calcLeastSquaresCPP(const double* x, const double* y, int n, double* m, double* b);
bool calcLeastSquaresAvxFma(double x[], double y[], int n, double* m, double* b);

extern "C" void avxAddPackedFp(double x[], double y[], int size, double* sumX, 
    double* sumY);
extern "C" void fmaPackedFp(double x[], double y[], int size, double* sumXX, 
    double* sumXY);

const int n = 512 * 512;
__declspec(align(64)) double x[n];
__declspec(align(64)) double y[n];

int main()
{
    // //std::cout << "Hello World!\n";
    // const int n = 6;
    // double x[n] = {0, 2, 4, 6, 8, 10};
    // double y[n] = {51.23, 34.6, 12.3, 56.8, 90.1, 111.9};

    /*const int n = 256 * 128;
    __declspec(align(64)) double x[n];
    __declspec(align(64)) double y[n];*/

    // Initialize the input vectors
	// All points lie on the line y = 1*x + 0.5
	const double slope = 1.0;
	const double y_int = 0.5;

	for (int i = 0; i < n; i += 4)
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

    double m1 = 0, m2 = 0, m3 = 0, b1 = 0, b2 = 0, b3 = 0;
    bool rv1, rv2, rv3;

    chrono::high_resolution_clock::time_point start = 
        chrono::high_resolution_clock::now();
    rv1 = calcLeastSquaresASM(x, y, n, &m1, &b1);
    chrono::high_resolution_clock::time_point end =
        chrono::high_resolution_clock::now();
    auto time1 = end - start;

    start = chrono::high_resolution_clock::now();
    rv2 = calcLeastSquaresCPP(x, y, n, &m2, &b2);
    end = chrono::high_resolution_clock::now();
    auto time2 = end - start;

    start = chrono::high_resolution_clock::now();
    rv3 = calcLeastSquaresAvxFma(x, y, n, &m3, &b3);
    end = chrono::high_resolution_clock::now();
    auto time3 = end - start;

    // for(int i=0; i<n; i++) {
    //     printf("x[%d] = %12.4lf, y[%d] = %12.4lf\n", i, x[i], i, y[i]);
    // }

    long long t1 = time1.count();
    long long t2 = time2.count();
    long long t3 = time3.count();

    printf("slope = %6.4lf, m1 = %6.4lf, m2 = %6.4lf, m3 = %6.4lf\n", 
        slope, m1, m2, m3);
    cout << "y_int = " << y_int << ", b1 = " << b1 << ", b2 = " << b2 << 
        ", b3 = " << b3 << endl;
    cout << "rv1 = " << rv1 << ", rv2 = " << rv2 << ", rv3 = " << rv3 << endl;
    // std::cout << "m1 = " << m1 << ", m2 = " << m2 << std::endl;
    // cout << "b1 = " << b1 << ", b2 = " << b2 << endl;
    // cout << "rv1 = " << rv1 << ", rv2 = " << rv2 << endl;

    /*printf(
    "\nTime t1 = %lld nsec. = %lld usec., t2 = %lld nsec. = %lld usec., t3 = %lld nsec. = %lld usec.\n", 
        t1, t1 / 1000, t2, t2 / 1000, t3, t3 / 1000);*/
    printf("\nTime t1 = %lld usec, t2 = %lld usec, t3 = %lld usec\n",
        t1 / 1000, t2 / 1000, t3 / 1000);
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

    /*cout << endl << "calcLeastSquaresCPP: " << "sumX = " << sumX << ", sumY = "
        << sumY << ", sumXX = " << sumXX << ", sumXY = " << sumXY << endl;*/

    double denom = n * sumXX - sumX * sumX;
    if (LsEpsilon >= fabs(denom)) return false;

    *m = (n * sumXY - sumX * sumY) / denom;
    *b = (sumXX*sumY - sumX*sumXY) / denom;

    return true;
}

bool calcLeastSquaresAvxFma(double x[], double y[], int n, double* m,
    double* b) {
    if (n <= 0) {
        return false;
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

    // chrono::high_resolution_clock::time_point start =
    //     chrono::high_resolution_clock::now();
    avxAddPackedFp(x, y, n, &sumX, &sumY);
    // chrono::high_resolution_clock::time_point end =
    //     chrono::high_resolution_clock::now();
    // auto timeSum = end - start;

    // start = chrono::high_resolution_clock::now();
    fmaPackedFp(x, y, n, &sumXX, &sumXY);
    // end = chrono::high_resolution_clock::now();
    // auto timeProd = end - start;

    // long long tSum = timeSum.count();
    // long long tProd = timeProd.count();

    // printf("\nTime tSum = %lld usec, tProd = %lld usec\n\n",
    //     tSum / 1000, tProd / 1000);

    double denom = n * sumXX - sumX * sumX;
    if (LsEpsilon >= fabs(denom)) return false;

    /*cout << endl << "calcLeastSquaresAvxFma: " << "sumX = " << sumX << ", sumY = " 
        << sumY << ", sumXX = " << sumXX << ", sumXY = " << sumXY << endl;*/

    *m = (n * sumXY - sumX * sumY) / denom;
    *b = (sumXX * sumY - sumX * sumXY) / denom;

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

#include <iostream>
#include <cmath>
#include <mpi.h>
using namespace std;

const double LsEpsilon = 1.0e-12;

// y = mx + b (Least Squares Calculation)
bool calcLeastSquaresMPI(const double* x, const double* y, int n, 
    double* m, double* b, int rank, int size) {
    
    int chunk_size = n / size;
    int remainder = n % size;
    int local_n = chunk_size + (rank < remainder ? 1 : 0);
    int offset = rank * chunk_size + min(rank, remainder);

    double local_sumX = 0.0, local_sumY = 0.0;
    double local_sumXY = 0.0, local_sumXX = 0.0;

    for (int i = 0; i < local_n; ++i) {
        double xi = x[offset + i];
        double yi = y[offset + i];
        local_sumX += xi;
        local_sumY += yi;
        local_sumXX += xi * xi;
        local_sumXY += xi * yi;
    }

    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0;

    MPI_Reduce(&local_sumX, &sumX, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sumY, &sumY, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sumXY, &sumXY, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sumXX, &sumXX, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double denom = n * sumXX - sumX * sumX;
        if (LsEpsilon >= fabs(denom)) {
            return false;
        }

        *m = (n * sumXY - sumX * sumY) / denom;
        *b = (sumXX * sumY - sumX * sumXY) / denom;

        return true;
    }

    return false; // other ranks don't compute final result
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int n = 5000 * 5000;
    const double slope = 1.0;
    const double y_int = 0.5;

    double* x = nullptr;
    double* y = nullptr;

    if (rank == 0) {
        x = new double[n];
        y = new double[n];

        #pragma omp parallel for
        for (int i = 0; i < n; i += 4) {
            x[i] = i;
            y[i] = slope * x[i] + y_int;
            x[i+1] = i + 1;
            y[i+1] = slope * x[i+1] + y_int;
            x[i+2] = i + 2;
            y[i+2] = slope * x[i+2] + y_int;
            x[i+3] = i + 3;
            y[i+3] = slope * x[i+3] + y_int;
        }
    }

    // Broadcast input arrays to all processes
    if (rank != 0) {
        x = new double[n];
        y = new double[n];
    }

    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double m1 = 0.0, b1 = 0.0;
    double start_time = MPI_Wtime();
    bool success = calcLeastSquaresMPI(x, y, n, &m1, &b1, rank, size);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        cout << "slope = " << slope << ", m1 = " << m1 << endl;
        cout << "y_int = " << y_int << ", b1 = " << b1 << endl;
        cout << "Success = " << success << endl;
        printf("Total time in calcLeastSquaresMPI: %.6f seconds\n", 
            end_time - start_time);
    }

    delete[] x;
    delete[] y;

    MPI_Finalize();
    return 0;
}
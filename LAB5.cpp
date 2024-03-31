#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

const int N = 10;

double RandomDouble(double min, double max) {
    return min + (max - min) * (rand() / static_cast<double>(RAND_MAX));
}

double* FindY(double* A, double* X, int n) {
    double* Y = (double*)malloc(n * sizeof(double));;
    for (int i = 0; i < n; i++) {
        Y[i] = 0;
        for (int j = 0; j < n; j++) {
            Y[i] += A[i * n + j] * X[j];
        }
    }
    return Y;
}

int main(int argc, char* argv[])
{	//теги сообщений 1 - от мастера, 2 - от рабочих
    
    double* X = new double[N];
    
    int ProcAmout, rank;
    MPI_Status status;
    double beginTime, endTime;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcAmout);

    int blockSize = N / ProcAmout; //Считаем, что количество строк делится на количество процессов без остатка

    double* matrixA = NULL;
    double* rightHandB = NULL;
    double* X2 = NULL;
    if (rank == 0) { //Мастер нормализует и отправляет
    beginTime = MPI_Wtime();
    for (int k = 0; k < blockSize; k++)
    {
        column = rank * blockSize + k;
        for (int i = column + 1; i < N; i++) {
            subA[k * N + i] = subA[k * N + i] / subA[k * N + column];
        }
        subB[k] /= subA[k * N + column];
        subA[k * N + column] = 1;
        b = subB[k];
        memcpy(row, &subA[k * N], sizeof(double) * N);
        MPI_Bcast(row, N, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        MPI_Bcast(&b, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        for (int j = k + 1; j < blockSize; j++) {
            for (int i = column + 1; i < N; i++) {
                subA[j * N + i] = subA[j * N + i] - subA[j * N + k] * row[i];

            }
            subB[j] -= subA[j * N + k] * b;
            subA[j * N + column] = 0;

        }

    }

}
  else { //Остальные процессы ждут строку от мастера
    for (int k = 0; k < rank * blockSize; k++) {//Каждый последующий процесс обрабатывает больше строк
        MPI_Bcast(row, N, MPI_DOUBLE, k / blockSize, MPI_COMM_WORLD);
        MPI_Bcast(&b, 1, MPI_DOUBLE, k / blockSize, MPI_COMM_WORLD);
        for (int j = 0; j < blockSize; j++) {
            for (int i = k + 1; i < N; i++) {
                subA[j * N + i] = subA[j * N + i] - subA[j * N + k] * row[i];
            }
            subB[j] -= subA[j * N + k] * b;
            subA[j * N + k] = 0;
        }
    }
  
    MPI_Finalize();
    
}

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
    if (rank == 0) { //Мастер генерирует 
        matrixA = (double*)malloc(N * N * sizeof(double*));
        for (int i = 0; i < N * N; i++) { 
            matrixA[i] = RandomDouble(-100, 100);
        }
        X2 = (double*)malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) { //Генерация вектора X2, с которым будет сравниваться вектор X, получаемый из метода Гаусса 
            X2[i] = RandomDouble(-100, 100);
        }
        rightHandB = (double*)malloc(N * sizeof(double));
        rightHandB = FindY(matrixA, X2, N);
    }
    double* subA = (double*)malloc(blockSize * N * sizeof(double*)); //Блок строк у каждого процесса
    double* subB = (double*)malloc(blockSize * sizeof(double));

    MPI_Scatter(matrixA, blockSize * N, MPI_DOUBLE, subA, blockSize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(rightHandB, blockSize, MPI_DOUBLE, subB, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double* row = (double*)malloc(N * sizeof(double));
    double b;
    int column;

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
                    subA[j * N + i] = subA[j * N + i] - subA[j * N + column] * row[i];

                }
                subB[j] -= subA[j * N + column] * b;
                subA[j * N + column] = 0;
            }

        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(subA, blockSize * N, MPI_DOUBLE, matrixA, blockSize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(subB, blockSize, MPI_DOUBLE, rightHandB, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int k = N - 1; k >= 0; k--) {
            X[k] = rightHandB[k];
            for (int i = N - 1; i > k; i--) {
                X[k] = X[k] - matrixA[k * N + i] * X[i];
            }
            X[k] = X[k] / matrixA[k * N + k];
        }
        endTime = MPI_Wtime();
        std::cout << "Time = " << endTime - beginTime << std::endl;
    }
    
    MPI_Finalize();
    if (rank == 0) {
        int c = 0;
        for (int i = 0; i < N; i++) {
            if (std::abs(X[i] - X2[i]) > 0.00001) {
                c++;
            }
        }
        std::cout << "Number of non - matching X = " << c << std::endl;
    }

}

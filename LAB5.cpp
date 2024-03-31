#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

const int N = 10;

double RandomDouble(double min, double max) {
    return min + (max - min) * (rand() / static_cast<double>(RAND_MAX));
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
        
    }
    
    MPI_Finalize();
    
}

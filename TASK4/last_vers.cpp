#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <malloc.h>
#include <math.h>
#include <map>
#include <set>
#include <mpi.h>
using namespace std;

//TASK4

double* slau_method_mpi(double* x, int N,  int number_iter, int rank, int size) {
    double criteria = 1;
    double tmp = 0;
    double tau = 0;
    double num0 = 0;
    double num1 = 0;
    double sum0 = 0;
    double sum1 = 0;
    double norma_b = 0;
    double sum[2] = { 0 };
    double epsilon = 1e-6;
    int iter = 0;
    double* y = new(nothrow) double[N];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* b = new(nothrow) double[N];
    if (b == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* A = new(nothrow) double[N * N];
    if (A == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    for (int i = 0; i < N; i++) {
        y[i] = 0;
        b[i] = N + 1;
        for (long long int j = 0; j < N; j++) {
            if (i == j) {
                A[i * N + j] = 2;
            }
            else {
                A[i * N + j] = 1;
            }
        }
    }

    int* sendcounts = new(nothrow) int[size];
    if (sendcounts == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 1");
        exit(-1);
    }
    int* displa = new(nothrow) int[size];
    if (displa == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 1");
        exit(-1);
    }

    int* sendcounts_matrix = new(nothrow) int[size];
    if (sendcounts == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 1");
        exit(-1);
    }
    int* displa_matrix = new(nothrow) int[size];
    if (displa == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 1");
        exit(-1);
    }
    int nmin = N / size;
    int nextra = N % size;
    int k = 0;
    int k_matrix = 0;
    
    for (int i = 0; i < size; i++)
    {
        sendcounts[i] = nmin;
        sendcounts_matrix[i] = nmin * N;
        if (i < nextra) {
            sendcounts[i] = nmin + 1;
            sendcounts_matrix[i] = (nmin + 1) * N;
        }

        displa[i] = k;
        k = k + sendcounts[i];
        displa_matrix[i] = k_matrix;
        k_matrix = k_matrix + sendcounts_matrix[i];
    }

    int size_A_local = sendcounts_matrix[rank];
    double* A_local = new(nothrow) double[size_A_local];
    if (A_local == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        delete[] A;
        delete[] b;
        delete[] x;
        exit(-1);
    }

    for (int i = 0; i < size_A_local; i++) {
        A_local[i] = 0;
    }


    /////////////////////////

    MPI_Scatterv(A, sendcounts_matrix, displa_matrix, MPI_DOUBLE, A_local, sendcounts_matrix[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
   // MPI_Scatterv(y, sendcounts, displa, MPI_DOUBLE, y_local, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #pragma omp parallel //%
    {//%

        #pragma omp  for reduction(+:iter)
        for (int i = 0; i < rank; i++) {
            iter += sendcounts[rank];
        }

        #pragma omp  for reduction(+:norma_b)//%
        for (int i = 0; i < N; i++) {
            norma_b += pow(b[i], 2);
        }

        #pragma omp single//%
        {//%
            num1 = sqrt(norma_b);
        }//%

        for (int j = 0; j < number_iter; j++) {

            //Yn = AXn - b
            #pragma omp  for private(tmp) reduction(+:num0)//%
            for (int i = 0; i < sendcounts[rank]; i++) {
                tmp = 0;
                for (int k = 0; k < N; k++) {
                    tmp += A_local[i * N + k] * x[k];

                }

                y[i + iter] = tmp - b[i];
                num0 += pow(tmp - b[i], 2);

            }

           // MPI_Allgatherv(y_local, sendcounts[rank], MPI_DOUBLE, y, sendcounts, displa, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y, sendcounts, displa, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &num0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            //CRITERIA 
            #pragma omp single//%
            {//%
                criteria = sqrt(num0) / num1;
                num0 = 0;
            }//%

            if (criteria < epsilon) {
                break;
            }

            //TAU
            #pragma omp  for private(tmp) reduction(+:sum0,sum1)//%
            for (int i = 0; i < sendcounts[rank]; i++) {
                tmp = 0;
                for (int k = 0; k < N; k++) {
                    tmp += A_local[i * N + k] * y[k];

                }

                sum0 += y[i] * tmp;
                sum1 += tmp * tmp;
            }
            #pragma omp single
            {
                sum[0] = sum0;
                sum[1] = sum1;
                sum0 = 0;
                sum1 = 0;
            }

            MPI_Allreduce(MPI_IN_PLACE, &sum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            #pragma omp single//%
            {//%
                
                if (sum[1] != 0) {
                    tau = sum[0] / sum[1];
                    sum[0] = 0;
                    sum[1] = 0;
                }
                else {
                    fprintf(stderr, "cannot calculate criteria");
                    exit(-1);
                }
            }//%

            //Xn+1 = Xn - taun*yn 
            #pragma omp  for private(tmp)//%
            for (int i = 0; i < N; i++) {
                tmp = x[i] - tau * y[i];
                x[i] = tmp;

            }

        }
    }//%

    delete[] y;
    delete[] sendcounts;
    delete[] sendcounts_matrix;
    delete[] displa;
    delete[] displa_matrix;
    
    delete[] A_local;
    delete[] b;
    delete[] A;

    return x;

}

int test_slau_simple(double* x, int N) {
    double tmp = 0;
    double eps = 1e-7;
    for (int i = 0; i < N; i++) {
        tmp += x[i];
    }
    if (fabs(tmp - N) > eps) {
        return -1;
    }
    return 0;
}

int main(int argc, char* argv[])
{

    if (argc < 4) {

        fprintf(stderr, "Enter the arguments in the format: program.exe N number_ter num_thread");
        return 1;
    }

    int i = 0;
    long long int N = 0;
    long long int number_iter = 0;
    int num_thread = 0;
    i = sscanf(argv[1], "%lld", &N);
    if (i != 1) {
        fprintf(stderr, "The size must be an long long integer");
        printf("%d", i);
        return 1;
    }
    
    i = sscanf(argv[2], "%lld", &number_iter);
    if (i != 1) {
        fprintf(stderr, "The size must be an long long integer");
        printf("%d", i);
        return 1;
    }

    i = sscanf(argv[3], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The size must be an integer");
        printf("%d", i);
        return 1;
    }

    double* x = new(nothrow) double[N];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        
        exit(-1);
    }

    if (num_thread > 8 || num_thread < 1) {
        num_thread = omp_get_max_threads();
        omp_set_num_threads(num_thread);
    }
    else {
        omp_set_num_threads(num_thread);
    }

    //////////////////// MPI INIT ////////////////////

    double start_mpi = 0;
    double end_mpi = 0;
    int ret_mpi = 0;
    int rank, size, provided = 0;
   
    MPI_Init_thread(&argc, &argv,
        MPI_THREAD_SERIALIZED, &provided);
    if (provided != MPI_THREAD_SERIALIZED) {
        fprintf(stderr, "MPI_THREAD_SERIALIZED not available\n");
        return EXIT_FAILURE;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (long long int i = 0; i < N; i++) {
        x[i] = 0;
    }

    start_mpi = MPI_Wtime();
    x = slau_method_mpi(x, N,  number_iter, rank, size);
    end_mpi = MPI_Wtime();

   
    if (rank == 0)
    {
        ret_mpi = test_slau_simple(x, N);

        if ( ret_mpi != 0) {
            printf("%lld\t%lf\t%d\t%d\n", N, end_mpi - start_mpi, 1, num_thread);
            delete[] x;
            MPI_Finalize();
            return 1;
        }
        else {
            printf("%lld\t%lf\t%d\t%d\n", N,  end_mpi - start_mpi, 0, num_thread);
            delete[] x;
            MPI_Finalize();
            return 0;
        }

    }
    MPI_Finalize();

}

#include <iostream>
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
#include <stdlib.h> 
#include <mpi.h>
using namespace std;

void scatter_vector(
    double local_x[],
    int part_size,
    double x[],
    int rank,
    MPI_Comm comm
) {
 
    if (rank == 0) {   
        MPI_Scatter(x, part_size, MPI_DOUBLE, local_x, part_size, MPI_DOUBLE, 0, comm);
    }
    else {
        MPI_Scatter(x, part_size, MPI_DOUBLE, local_x, part_size, MPI_DOUBLE, 0, comm);
    }
}

void gather_vector(
    double local_res[],
    int part_size,
    double vector[],
    int my_rank,
    MPI_Comm comm
) {
   
    if (my_rank == 0) {
        MPI_Gather(local_res, part_size, MPI_DOUBLE, vector, part_size, MPI_DOUBLE, 0, comm);
    }
    else {
        MPI_Gather(local_res, part_size, MPI_DOUBLE, vector, part_size, MPI_DOUBLE, 0, comm);
    }
}

void vector_sum(
    double local_x[],   
    double local_y[],   
    double local_res[],
    int part_size,
    double a

) {
    for (int local_i = 0; local_i < part_size; local_i++) {
        local_res[local_i] = a*local_x[local_i] + local_y[local_i];
    }
}

int test_blas_daxpy(double* y, double* y_etalon, long long int vector_size) {
    double eps = 1e-07;
    double a = 0;
    double b = 0;
    for (long long int i = 0; i < vector_size; i++) {
        if (fabs(y[i] - y_etalon[i]) > eps) {
            return(-1);
        }
    }
    return 0;
}

int main(int argc, char* argv[])
{
    /*if (argc < 3) {

        fprintf(stderr, "Enter the arguments in the format: program.exe vector_size a ");
        return 1;
    }



    long long int vector_size = 0;
    int i = 0;
    i = sscanf_s(argv[1], "%lld", &vector_size);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an long long integer");
        printf("%d", i);
        return 1;
    }
    i = 0;

    double a = 0;
    i = sscanf_s(argv[2], "%lf", &a);
    if (i != 1) {
        fprintf(stderr, "The a must be an double");
        printf("%d", i);
        return 1;
    }*/
   
    long long int vector_size = 1000; 
    double a = 4;

    double* x = new(nothrow) double[vector_size];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 1");
        return 1;
    }


    double* y = new(nothrow) double[vector_size];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 2");
        delete[] x;
        return 1;
    }

    double* y_etalon = new(nothrow) double[vector_size];
    if (y_etalon == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 3");
        delete[] x;
        delete[] y;
        return 1;
    }

    double* y_res = new(nothrow) double[vector_size];
    if (y_res == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        return 1;
    }

   
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = i * 2;

        y_etalon[i] = a * i * 2;

    }


    double start_mpi = 0;
    double end_mpi = 0;

    long long int partVectorSize = 0;

    

    int rank, size, provided;
    MPI_Init_thread(&argc, &argv,
        MPI_THREAD_SERIALIZED, &provided);
    if (provided != MPI_THREAD_SERIALIZED) {
        fprintf(stderr, "MPI_THREAD_SERIALIZED not available\n");
        return EXIT_FAILURE;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    partVectorSize = vector_size / size;

    double* x_local = new(nothrow) double[partVectorSize];
    if (x_local == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 3");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        delete[] y_res;
        MPI_Finalize();
        return 1;
    }

    double* y_local = new(nothrow) double[partVectorSize];
    if (y_local == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        delete[] y_res;
        delete[] x_local;
        MPI_Finalize();
        return 1;
    }
    double* y_local_res = new(nothrow) double[partVectorSize];
    if (y_local == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        delete[] y_res;
        delete[] x_local;
        delete[] y_local;
        MPI_Finalize();
        return 1;
    }


    start_mpi = MPI_Wtime();

    scatter_vector(x_local, partVectorSize, x, rank, MPI_COMM_WORLD );
   
    scatter_vector(y_local, partVectorSize, y, rank, MPI_COMM_WORLD);

    vector_sum(x_local, y_local, y_local_res, partVectorSize,a);

    gather_vector(y_local_res, partVectorSize, y_res, rank, MPI_COMM_WORLD);
 
    end_mpi = MPI_Wtime();


    int res_test_mpi = test_blas_daxpy(y_res, y_etalon, vector_size);

    if (res_test_mpi != 0) {
        printf("%lld\t%lf\t%d\n", vector_size, end_mpi - start_mpi, 1);
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        delete[] y_res;
        delete[] y_local_res;
        delete[] x_local;
        delete[] y_local;
        MPI_Finalize();
        return 1;
    }
    else {
        printf("%lld\t%lf\t%d\n", vector_size, end_mpi - start_mpi, 0);
    }

    delete[] x;
    delete[] y;
    delete[] y_etalon;
    delete[] y_res;
    delete[] y_local_res;
    delete[] x_local;
    delete[] y_local;
    MPI_Finalize();
    return 0;

}



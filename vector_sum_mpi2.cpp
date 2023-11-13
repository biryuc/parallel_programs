#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
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

//TASK2
double* blas_daxpy(double* x, double* y, long long int vector_size, double a) {

   
    for (long long int i = 0; i < vector_size; i++) {
        y[i] = a * x[i] + y[i];
       
    }

    return y;
}


//double* blas_daxpy_omp(double* x, double* y,double* y_res, long long int vector_size, double a) {
//
//    double tmp = 0;
//
//    #pragma omp parallel for private(tmp)
//    for (long long int i = 0; i < vector_size; i++) {
//            tmp = a * x[i] + y[i];
//            y_res[i] = tmp;
//    }
//
//    return y_res;
//}
double* blas_daxpy_omp(double* x, double* y, long long int vector_size, double a) {

    double tmp = 0;

    #pragma omp parallel for private(tmp)
    for (long long int i = 0; i < vector_size; i++) {
        tmp = a * x[i] + y[i];
        y[i] = tmp;
    }

    return y;
}

//////////////////////// MPI ///////////////////////////

double* blas_daxpy_mpi(double* x, double* y, double* y_res, long long int vector_size, double a,int rank,long long int partVectorSize,double* x_local, double* y_local,double* y_local_res) {

    MPI_Scatter(x, partVectorSize, MPI_DOUBLE, x_local, partVectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, partVectorSize, MPI_DOUBLE, y_local, partVectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double tmp = 0;
    #pragma omp parallel for private(tmp)
    for (int local_i = 0; local_i < partVectorSize; local_i++) {
        tmp = a * x_local[local_i] + y_local[local_i];
        y_local_res[local_i] = tmp;
    }

    
    MPI_Gather(y_local_res, partVectorSize, MPI_DOUBLE, y_res, partVectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    return y_res;
}

int test_blas_daxpy(double* y, double* y_etalon, long long int vector_size) {
    double eps = 1e-07;
    double a = 0;
    double b = 0;
    for (long long int i = 0; i < vector_size; i++) {
        if (fabs(y[i] - y_etalon[ i]) > eps) {
            return(-1);
        }
    }
    return 0;
}

int main(int argc, char* argv[])
{

    if (argc < 4) {

        fprintf(stderr, "Enter the arguments in the format: program.exe vector_size a num_thread");
        return 1;
    }



    long long int vector_size = 0;
    int i = 0;
    i = sscanf(argv[1], "%lld", &vector_size);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an long long integer");
        printf("%d", i);
        return 1;
    }
    i = 0;

    double a = 0;
    i = sscanf(argv[2], "%lf", &a);
    if (i != 1) {
        fprintf(stderr, "The a must be an double");
        printf("%d", i);
        return 1;
    }
    int num_thread = 0;
    i = sscanf(argv[3], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The num_thread must be an int");
        printf("%d", i);
        return 1;
    }
  
    
    int res_test_first = 0;
    int res_test_second = 0;
    int res_test_third = 0;
    int res_test_omp = 0;
    int res_test_mpi = 0;

    double start_omp = 0;
    double end_omp = 0;
  
    double start = 0;
    double end = 0;


    if (num_thread < 1) {
        num_thread = omp_get_max_threads();
        omp_set_num_threads(num_thread);
    }
    else {
        omp_set_num_threads(num_thread);
    }
   
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
   
    //double* y_res = new(nothrow) double[vector_size];
    //if (y_res == nullptr) {
    //    fprintf(stderr, "Memory cannot be allocated 4");
    //    delete[] x;
    //    delete[] y;
    //    delete[] y_etalon;
    //    return 1;
    //}

    double* y_mpi = new(nothrow) double[vector_size];
    if (y_mpi == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        return 1;
    }


    double start_mpi = 0;
    double end_mpi = 0;

    long long int partVectorSize = 0;



    int rank, size, provided = 0;
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
        //delete[] y_res;
        MPI_Finalize();
        return 1;
    }

    double* y_local = new(nothrow) double[partVectorSize];
    if (y_local == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        //delete[] y_res;
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
       // delete[] y_res;
        delete[] x_local;
        delete[] y_local;
        MPI_Finalize();
        return 1;
    }

    //TEST3
    //i*2
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = i*2;
        y_etalon[i] = a*i*2 ;

    }

    start = omp_get_wtime();
    y = blas_daxpy(x, y, vector_size, a);
    end= omp_get_wtime();
    
    res_test_third = test_blas_daxpy(y, y_etalon, vector_size);
   
    ////////////////// OMP ////////////////////////
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = i * 2;
       // y_res[i] = 0;
        y_etalon[i] = a * i * 2;

    }
    start_omp = omp_get_wtime();
    y = blas_daxpy_omp(x, y, vector_size, a);
    end_omp = omp_get_wtime();
  
    res_test_omp = test_blas_daxpy(y, y_etalon, vector_size);

    ////////////////// MPI ////////////////////////
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = i * 2;
        y_mpi[i] = 0;
        y_etalon[i] = a * i * 2;

       

    }
    for (int i = 0; i < partVectorSize; i++) {
        x_local[i] = 0;
        y_local[i] = 0;
        y_local_res[i] = 0;
    }

    start_mpi = MPI_Wtime();
    y_mpi = blas_daxpy_mpi(x, y, y_mpi, vector_size, a, rank, partVectorSize, x_local, y_local, y_local_res);
    end_mpi = MPI_Wtime();

    res_test_mpi = test_blas_daxpy(y_mpi, y_etalon, vector_size);
 
    if (rank == 0) {


        if (res_test_third != 0 || res_test_omp != 0 || res_test_mpi != 0) {
            printf("%lld\t%lf\t%lf\t%lf\t%d\t%d\n", vector_size, end - start, end_omp - start_omp, end_mpi - start_mpi, 1, num_thread);
            delete[] x;
            delete[] y;
            delete[] y_etalon;
            //delete[] y_res;
            delete[] y_local_res;
            delete[] x_local;
            delete[] y_local;
            MPI_Finalize();
            return 1;
        }
        else {
            printf("%lld\t%lf\t%lf\t%lf\t%d\t%d\n", vector_size, end - start, end_omp - start_omp, end_mpi - start_mpi, 0, num_thread);
        }
    }
    
    delete[] x;
    delete[] y;
    delete[] y_etalon;
    //delete[] y_res;
    delete[] y_local_res;
    delete[] x_local;
    delete[] y_local;
    MPI_Finalize();
    return 0;
}



// vector_sum.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

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
using namespace std;

//TASK2
double* blas_daxpy(double* x, double* y, long long int vector_size, double a) {

   
    for (long long int i = 0; i < vector_size; i++) {
        y[i] = a * x[i] + y[i];
       
    }

    return y;
}


double* blas_daxpy_omp(double* x, double* y,double* y_res, long long int vector_size, double a) {

    double tmp = 0;

    #pragma omp parallel for private(tmp)
    for (long long int i = 0; i < vector_size; i++) {
            tmp = a * x[i] + y[i];
            y_res[i] = tmp;
    }

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
    }
    int num_thread = 0;
    i = sscanf_s(argv[3], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The num_thread must be an int");
        printf("%d", i);
        return 1;
    }

    //int a = 4;
    //int vector_size = 100000000;
    //int num_thread = 8;

       
    int r = 0;

    int res_test_first = 0;
    int res_test_second = 0;
    int res_test_third = 0;
    int res_test_omp = 0;

    double start_omp = 0;
    double end_omp = 0;
    //int num_threads = 0;
    double start = 0;
    double end = 0;


    if (num_thread > 8 || num_thread < 1) {
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
   
    double* y_res = new(nothrow) double[vector_size];
    if (y_res == nullptr) {
        fprintf(stderr, "Memory cannot be allocated 4");
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        return 1;
    }
   
    //TEST1
    //zero vector
    //for (int i = 0; i < vector_size; i++) {
    //   
    //    y[i] = 0;
    //    x[i] = 0;
    //    y_etalon[i] = 0;

    //}

    //y = blas_daxpy(x, y, vector_size, a);
    //res_test_first = test_blas_daxpy(y, y_etalon, vector_size);

    ////TEST2
    ////1 vector
    //for (int i = 0; i < vector_size; i++) {

    //    y[i] = 0;
    //   
    //    x[i] = 1;
    //    y_etalon[i] = a;

    //}
   
    //y = blas_daxpy(x, y, vector_size, a);
    //res_test_second =  test_blas_daxpy(y, y_etalon, vector_size);

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
   
    //OMP
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = i * 2;
        y_res[i] = 0;
        y_etalon[i] = a * i * 2;

    }
    start_omp = omp_get_wtime();
    y_res = blas_daxpy_omp(x, y, y_res, vector_size, a);
    end_omp = omp_get_wtime();
  
    res_test_omp = test_blas_daxpy(y_res, y_etalon, vector_size);
 

    if (res_test_third != 0 || res_test_second != 0 || res_test_first != 0 || res_test_omp != 0) {
        printf("%lld\t%lf\t%lf\t%d\t%d\n", vector_size, end - start, end_omp - start_omp, 1,num_thread);
        delete[] x;
        delete[] y;
        delete[] y_etalon;
        delete[] y_res;

        return 1;
    }
    else {
        printf("%lld\t%lf\t%lf\t%d\t%d\n", vector_size, end - start, end_omp - start_omp, 0, num_thread);
    }
    
    delete[] x;
    delete[] y;
    delete[] y_etalon;
    delete[] y_res;

    return 0;
}



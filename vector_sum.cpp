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

using namespace std;

//TASK2
double* blas_daxpy(double* x, double* y, int vector_size, int a) {
    for (int i = 0; i < vector_size; i++) {
        y[i] = a * x[i] + y[i];
    }

    return y;
}

int test_blas_daxpy(double* y, double* y_etalon,int vector_size) {
    double eps = 1e-07;
    for (int i = 0; i < vector_size; i++) {
        if (fabs(y[i] - y_etalon[ i]) > eps) {
            return(-1);
        }
    }
    return 0;
}

int main()
{
       
    int r = 0;
    int a = 10;
    int res_test = 0;
    int vector_size = 1000000;
    double* x = new(nothrow) double[vector_size];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* y = new(nothrow) double[vector_size];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    double* y_etalon = new(nothrow) double[vector_size];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    //TEST1
    //zero vector
    for (int i = 0; i < vector_size; i++) {
       
        y[i] = 0;
        x[i] = 0;
        y_etalon[i] = 0;

    }

    y = blas_daxpy(x, y, vector_size, a);


    res_test = test_blas_daxpy(y, y_etalon, vector_size);
    if (res_test != 0) {
        fprintf(stderr, "Test 1 failed,The vectors are different\n");

    }
    else {
        fprintf(stderr, "Test 1 ok\n");
    }


    //TEST2
    //1 vector
    for (int i = 0; i < vector_size; i++) {

        y[i] = 0;
        x[i] = 1;
        y_etalon[i] = 10;

    }

    y = blas_daxpy(x, y, vector_size, a);

    res_test =  test_blas_daxpy(y, y_etalon, vector_size);
    if (res_test != 0) {
        fprintf(stderr, "Test 2 failed,The vectors are different\n");
    }
    else {
        fprintf(stderr, "Test 2 ok\n");
    }


   

    
    delete[] x;
    delete[] y;
}



// slau_method.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


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

//TASK4

double* A_vec_mult(int* A,double*res, int N, double* x) {
    double tmp = 0;
   
    for (int i = 0; i < N; i++) {
        tmp = 0;
        for (int k = 0; k < N; k++) {
            tmp += A[i * N + k] * x[k];
            res[i] = tmp;
            
        }
    }

    return res;
}
double* slau_method(int* A, int* b, double* x, double* y, int N, double epsilon) {
    double criteria = 1;
    double tau = 0;
    double sum0 = 0;
    double sum1 = 0;
    double num0 = 0;
    double num1 = 0;

    double* res = new(nothrow) double[N];
    if (res == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    while (criteria > epsilon) {
        //Y = AX - b
        res = A_vec_mult(A,res, N, x);
        for (int i = 0; i < N; i++) {
            
            y[i] = res[i] - b[i];

        }
        //TAU
        res = A_vec_mult(A,res, N, y);

        for (int i = 0; i < N; i++) {
            sum0 += y[i] * res[i];
            sum1 += res[i] * res[i];
        }
        if (sum1 != 0) {
            tau = sum0 / sum1;
            sum0 = 0;
            sum1 = 0;
        }


        //X^n
        for (int i = 0; i < N; i++) {
            x[i] = x[i] - tau * y[i];
        }

        //CRITERIA
         num0 = 0;
         num1 = 0;
        res = A_vec_mult(A,res, N, x);
        for (int i = 0; i < N; i++) {
            num0 += pow(res[i] - b[i], 2);
            num1 += pow(b[i], 2);
        }
        num0 = sqrt(num0);
        num1 = sqrt(num1);
        if (num1 != 0) {
            criteria = num0 / num1;
            num0 = 0;
            num1 = 0;
        }
    }

    delete[] res;
    return x;

}

int test_slau(double* x,double* etalon,int N) {
    double eps = 1e-7;
    for (int i = 0; i < N; i++) {
        if (fabs(x[i] - etalon[i]) > eps) {
            return -1;
        }
    }
    return 0;
}
int main()
{

    int N = 30;
    double epsilon = 1e-6;
    int* A = new(nothrow) int[N * N];
    if (A == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    int* b = new(nothrow) int[N];
    if (b == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* x = new(nothrow) double[N];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* etalon = new(nothrow) double[N];
    if (etalon == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* y = new(nothrow) double[N];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    //Инициализация
    for (int i = 0; i < N; i++) {
        b[i] = N + 1;
        x[i] = 0;
        y[i] = 0;
        etalon[i] = 1;
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i * N + j] = 2;
            }
            else {
                A[i * N + j] = 1;
            }
        }
    }

    x = slau_method(A, b, x, y, N, epsilon);

    int ret = test_slau(x, etalon, N);

    if (ret != 0) {
        fprintf(stderr, "Test failed");
    }
    else {
        fprintf(stderr, "Test ok");
    }


    //for (int i = 0; i < N; i++) {
    //    printf("%lf\n", x[i]);
    //}

    delete[] A;
    delete[] b;
    delete[] x;
    delete[] y;
    delete[] etalon;

}



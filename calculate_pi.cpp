// calculate_pi.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
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
//TASK1
double calculate_pi(int N) {
    double pi = 0;

    for (int i = 0; i < N; i++) {
        pi += 4 * (pow(-1, i) / (2 * i + 1));
    }
    return pi;
}


//TASK2
int* blas_daxpy(int* x,int* y,int vector_size,int a) {
    for (int i = 0; i < vector_size; i++) {
        y[i] = a * x[i] + y[i];
    }

    return y;
}

//TASK3
/* первая матрица, размером k столбцов на m строк*/
/* вторая матрица, размером n столбцов на k строк*/
//принцип локальности
void matrix_mult(float* matrix, float* matrix1, float* matrix2, int n, int k, int m) {

    float tmp = 0.0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i * n + j] = 0;
            for (int f = 0; f < k; f++) {
                tmp = matrix1[i * k + f] * matrix2[f * n + j];
                matrix[i * n + j] += tmp;
            }
        }
    }


}

//TASK4

float* A_vec_mult(int*A,int N,float* vector) {
    int tmp = 0;
    float* res_vec = new(nothrow) float[N];
    if (res_vec == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            tmp += A[i * N + k] * vector[k];
            res_vec[i] = tmp;
        }
    }
   

    return res_vec;
}
float* slau_method(int N,float epsilon) {
    float* tmp = new(nothrow) float[N];
    if (tmp == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    
    
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
    float* x = new(nothrow) float[N];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    float* y = new(nothrow) float[N];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    //Инициализация
    for (int i = 0; i < N; i++) {
        b[i] = N + 1;
        x[i] = 0;
        y[i] = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i * N + j] = 2;
            }
            else {
                A[i * N + j] = 1;
            }
        }
    }
    float criteria = 2*epsilon;
    float tau = 0;
    float sum0 = 0;
    float sum1 = 0;

    while (criteria > epsilon) {
        //Y = AX - b
        for (int i = 0; i < N; i++) {
            tmp = A_vec_mult(A, N, x);
            y[i] = tmp[i] - b[i];

        }
        //TAU
        tmp = A_vec_mult(A, N, y);

        for (int i = 0; i < N; i++) {
            sum0 += y[i] * tmp[i];
            sum1 += tmp[i] * tmp[i];
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
        float num0 = 0;
        float num1 = 0;
        tmp = A_vec_mult(A, N, x);
        for (int i = 0; i < N; i++) {
            num0 += pow(tmp[i] - b[i], 2);
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

   
    return x;

}

int main()
{
    int N = 3;
    float epsilon = 0.5;
    float* x = new(nothrow) float[N];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }


    x = slau_method(N, epsilon);

    for (int i = 0; i < N; i++) {
        printf("%f\n", x[i]);
    }

    return 0;
}



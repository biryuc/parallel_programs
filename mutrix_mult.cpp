// mutrix_mult.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
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

//TASK3
/* первая матрица, размером k столбцов на m строк*/
/* вторая матрица, размером n столбцов на k строк*/
//принцип локальности
void matrix_mult(double* matrix, double* matrix1, double* matrix2, int n, int k, int m) {

    for (int i = 0; i < m; i++) {
        for (int f = 0; f < k; f++) {
            for (int j = 0; j < n; j++) {
                    matrix[i * n + j] += matrix1[i * k + f] * matrix2[f * n + j];
            }
        }
    }

}

int test_matrix_mult(double* matrix, double* matrix_etalon,  int n,  int m) {
    double epsilon = 1e-7;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(matrix[i * n + j] - matrix_etalon[i * n + j]) >epsilon) {
                return -1;
            }
           
        }
    }
    return 0;
}

int main()
{   
    int n = 3;
    int m = 3;
    int k = 3;


    double* x = new(nothrow) double[m*k];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* y = new(nothrow) double[n*k];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    double* result = new(nothrow) double[m*n];
    if (result == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* etalon = new(nothrow) double[m * n];
    if (result == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }

    for (int i = 0; i < m ; i++) {
        for (int j = 0; j < n; j++) {
            x[i*m+j] = i * i;
            if (i == j) {
                y[i * m + j] = 1;

            }
            else {
                y[i * m + j] = 0;
            }
            result[i * m + j] = 0;
            etalon[i * m + j] = i * i;
        }

    }

    matrix_mult(result, x, y, n, k, m);

    int ret = test_matrix_mult(result, etalon, n, m);

    if (ret != 0) {
        fprintf(stderr, "Test failed");
    }
    else {
        fprintf(stderr, "Test ok");
    }

   /* for (int i = 0; i < m * n; i++) {
        printf("%lf\n", result[i]);
    }*/

}


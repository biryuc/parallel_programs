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
double* matrix_mult(double* matrix, double* matrix1, double* matrix2, long long int n, long long int k, long long int m) {
    for (long long int i = 0; i < m; i++) {
        for (long long int f = 0; f < k; f++) {
            for (long long int j = 0; j < n; j++) {
                matrix[i * n + j] = 0;
            }
        }
    }
    
    for (long long int i = 0; i < m; i++) {
        for (long long int f = 0; f < k; f++) {
            for (long long int j = 0; j < n; j++) {
                    
                    matrix[i * n + j] += matrix1[i * k + f] * matrix2[f * n + j];
            }
        }
    }

    return matrix;

}


double* matrix_mult_omp(double* matrix, double* matrix1, double* matrix2, long long int n, long long int k, long long int m) {

   
    
    for (long long int i = 0; i < m; i++) {
        for (long long int f = 0; f < k; f++) {
            for (long long int j = 0; j < n; j++) {
                matrix[i * n + j] = 0;
            }
        }
    }
   
        long long int i = 0, j = 0, f = 0;

        #pragma omp parallel for shared(matrix1, matrix2, matrix) private(i, j, f)
        for (i = 0; i < m; i++) {
            for (f = 0; f < k; f++) {
                for (j = 0; j < n; j++) {
                    matrix[i * n + j] += matrix1[i * k + f] * matrix2[f * n + j];
                }
            }
        }

    

    return matrix;

}

int test_matrix_mult(double* matrix, double* matrix_etalon, long long int n, long long  int m) {
    double epsilon = 1e-7;

    for (long long int i = 0; i < m; i++) {
        for (long long int j = 0; j < n; j++) {
            if (fabs(matrix[i * n + j] - matrix_etalon[i * n + j]) >epsilon) {
                return -1;
            }
           
        }
    }
    return 0;
}

int main(int argc, char * argv[])
{   

    if (argc < 4) {

        fprintf(stderr, "Enter the arguments in the format: program.exe n m k num_thread");
        return 1;
    }

    
    int i = 0;
    int n = 0;
    int m = 0;
    int k = 0;
    int num_thread = 0;

    i = sscanf_s(argv[1], "%d", &n);
    if (i != 1) {
        fprintf(stderr, "The size must be an long long integer");
        printf("%d", i);
        return 1;
    }
    i = sscanf_s(argv[2], "%d", &m);
    if (i != 1) {
        fprintf(stderr, "The size must be an long long integer");
        printf("%d", i);
        return 1;
    }
    i = sscanf_s(argv[3], "%d", &k);
    if (i != 1) {
        fprintf(stderr, "The size must be an long long integer");
        printf("%d", i);
        return 1;
    }

    i = sscanf_s(argv[4], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The size must be an integer");
        printf("%d", i);
        return 1;
    }
    //int n = 100;
    //int m = 100;
    //int k = 100;
    //int num_thread = 8;



    double* x = new(nothrow) double[m*k];
    if (x == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        exit(-1);
    }
    double* y = new(nothrow) double[n*k];
    if (y == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        delete[] x;
        exit(-1);
    }

    double* result = new(nothrow) double[m*n];
    if (result == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        delete[] x;
        delete[] y;
        exit(-1);
    }

    double* result_omp = new(nothrow) double[m * n];
    if (result_omp == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        delete[] x;
        delete[] y;
        delete[] result;
        exit(-1);
    }

    double* etalon = new(nothrow) double[m * n];
    if (etalon == nullptr) {
        fprintf(stderr, "Memory cannot be allocated");
        delete[] x;
        delete[] y;
        delete[] result;
        delete[] result_omp;
        exit(-1);
    }

    if (num_thread > 8 || num_thread < 1) {
        num_thread = omp_get_max_threads();
        omp_set_num_threads(num_thread);
    }
    else {
        omp_set_num_threads(num_thread);
    }

    double pi_omp = 0;
    double start_omp = 0;
    double end_omp = 0;
    double start = 0;
    double end = 0;
   

    for (int i = 0; i < m ; i++) {
        for (int j = 0; j < n; j++) {
            x[i*m+j] = i * i;
            if (i == j) {
                y[i * m + j] = 1;

            }
            else {
                y[i * m + j] = 0;
            }
            //result[i * m + j] = 0;
            //result_omp[i * m + j] = 0;
            etalon[i * m + j] = i * i;
        }

    }


    start = omp_get_wtime();
    matrix_mult(result, x, y, n, k, m);
    end = omp_get_wtime();

    int ret = test_matrix_mult(result, etalon, n, m);

    start_omp = omp_get_wtime();
    matrix_mult_omp(result_omp, x, y, n, k, m);
    end_omp = omp_get_wtime();

    int ret_omp = test_matrix_mult(result_omp, etalon, n, m);

    if (ret != 0 || ret_omp != 0) {
        printf("%lld\t%lf\t%lf\t%d\t%d\n", n*m, end - start, end_omp - start_omp, 1,num_thread);
        delete[] x;
        delete[] y;
        delete[] result;
        delete[] result_omp;
        delete[] etalon;
        return 1;
    }
    else {
        printf("%lld\t%lf\t%lf\t%d\t%d\n", n*m, end - start, end_omp - start_omp, 0,num_thread);
        delete[] x;
        delete[] y;
        delete[] result;
        delete[] result_omp;
        delete[] etalon;
        return 0;
    }

   /* for (int i = 0; i < m * n; i++) {
        printf("%lf\n", result[i]);
    }*/
   

}


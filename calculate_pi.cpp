#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <malloc.h>
#define _CRT_SECURE_NO_WARNINGS
#include <map>
#include <set>
#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;
//TASK1
double calculate_pi(long long int N) {
    double pi = 0;

    for (long long i = 0; i < N; i++) {
        pi += (pow(-1, i) / (2 * i + 1));
    }
    return 4*pi;
}


double calculate_pi_omp(long long int N) {
    double pi = 0;
   // schedule(guided,1 )
    #pragma omp parallel for default(none) shared(N) reduction(+:pi) 
        for (long long i = 0; i < N; i++) {
            pi += (pow(-1, i) / (2 * i + 1));

        }

    
   
    return 4 * pi;
}


int main(int argc , char* argv[])
{

    if (argc < 3) {
        
        fprintf(stderr, "Enter the arguments in the format: program.exe number_iteration num_thread");
        return 1;
    }

    long long int N = 0;
    int i = 0;
    int num_thread = 0;

    i = sscanf_s(argv[1], "%lld", &N);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an integer");
        printf("%d",i);
        return 1;
    }

    i = sscanf_s(argv[2], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an integer");
        printf("%d", i);
        return 1;
    }

   
    double pi = 0;
    double start = 0;
    double end = 0;

    double pi_omp = 0;
    double start_omp = 0;
    double end_omp = 0;
   // int num_threads = 0;

   
    start = omp_get_wtime();
    pi = calculate_pi(N);
    end = omp_get_wtime();

    if (num_thread > 8 || num_thread < 1) {
        num_thread = omp_get_max_threads();
        omp_set_num_threads(num_thread);
    }
    else {
        omp_set_num_threads(num_thread);
    }
    

    start_omp = omp_get_wtime();
    pi_omp = calculate_pi_omp(N);   
    end_omp = omp_get_wtime();


    int k = 0;
    if (fabs(pi - M_PI) > 0.01 || fabs(pi_omp - M_PI) > 0.01) {
       
        k = 1;
        printf("%lld\t%lf\t%lf\t%d\t%d\n", N, end - start, end_omp - start_omp, k, num_thread);
        return 1;
    }
    else {
        k = 0;
        printf("%lld\t%lf\t%lf\t%d\t%d\n", N, end - start, end_omp - start_omp, k, num_thread);
    }

    return 0;
}


//sscanf - возвращает количество сконвертированных значений

//omp_get_wtime()
//тестирование проводить в скрипте и записывать в файл , потом строить график в зависимости от количества потока, также в таблицу выводить 0 если правильно отработало 
//массив ptr_diff_t или size_t

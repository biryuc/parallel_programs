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
#include <mpi.h>


using namespace std;

double power(int n, int p) {
    double tmp = 1;
    for (int i = 0; i < p; i++) {
        tmp *= n;
    }
    return tmp;
}
//TASK1
double calculate_pi(long long int N) {
    double pi = 0;

    for (long long i = 0; i < N; i++) {
        pi += (power(-1, i) / (2 * i + 1));
    }
    return 4*pi;
}


double calculate_pi_omp(long long int N) {
    double pi = 0;
   // schedule(guided,1 )
    #pragma omp parallel for default(none) shared(N) reduction(+:pi) 
        for (long long i = 0; i < N; i++) {
            pi += (power(-1, i) / (2 * i + 1));

        }

    
   
    return 4 * pi;
}


//double calculate_pi_mpi(long long int N,int rank,int size) {
//    double pi = 0;
//    double tmp = 0;
//    
//   
//    for (long long i = rank; i < N; i+=size) {
//        tmp += (power(-1, i) / (2 * i + 1));
//
//    }
//
//    MPI_Reduce(&tmp, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//    return 4 * pi;
//}



int main(int argc , char* argv[])
{

   /* if (argc < 3) {
        
        fprintf(stderr, "Enter the arguments in the format: program.exe number_iteration num_thread");
        return 1;
    }*/
  

    long long int N = 0;
    int i = 0;
    int num_thread = 0;

    /*i = sscanf_s(argv[1], "%lld", &N);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an integer");
        printf("%d",i);
        MPI_Finalize();
        return 1;
    }

    i = sscanf_s(argv[2], "%d", &num_thread);
    if (i != 1) {
        fprintf(stderr, "The number of iterations must be an integer");
        printf("%d", i);
        MPI_Finalize();
        return 1;
    }*/

    N = 1000;
    num_thread = 8;
   
    double pi = 0;
    double start = 0;
    double end = 0;

    double pi_omp = 0;
    double start_omp = 0;
    double end_omp = 0;

    double start_mpi = 0;
    double end_mpi = 0;
    double pi_mpi = 0;

    /*int num_thread = 3;
    long long int N = 10;*/

    

   
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



    //////////////////// MPI ///////////////////////
    int size, rank, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
    if (provided != MPI_THREAD_SERIALIZED) {
        fprintf(stderr, "MPI_THREAD_SERIALIZED not available\n");
        return EXIT_FAILURE;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    start_mpi = MPI_Wtime();
    double tmp = 0;
    for (long long i = rank; i < N; i += size) {
        tmp += (power(-1, i) / (2 * i + 1));

    }
    MPI_Reduce(&tmp, &pi_mpi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end_mpi = MPI_Wtime();

    if (rank == 0)
    {
        pi_mpi = 4 *pi_mpi;
    }
    //////////////////// MPI ///////////////////////



    
    if (fabs(pi - M_PI) > 0.01 || fabs(pi_omp - M_PI) > 0.01 || fabs(pi_mpi - M_PI) > 0.01) {
       
        
        printf("%lld\t%lf\t%lf\t%lf\t%d\t%d\n", N, end - start, end_omp - start_omp,end_mpi - start_mpi, 1, num_thread);
        MPI_Finalize();
        return 1;
    }
    else {
       
        printf("%lld\t%lf\t%lf\t%lf\t%d\t%d\n", N, end - start, end_omp - start_omp, end_mpi - start_mpi, 0, num_thread);
    }


    MPI_Finalize();
    return 0;
}


//sscanf - возвращает количество сконвертированных значений

//omp_get_wtime()
//тестирование проводить в скрипте и записывать в файл , потом строить график в зависимости от количества потока, также в таблицу выводить 0 если правильно отработало 
//массив ptr_diff_t или size_t

//omp single

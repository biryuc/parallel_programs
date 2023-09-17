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

#include <map>
#include <set>

#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;
//TASK1
double calculate_pi(int N) {
    double pi = 0;

    for (int i = 0; i < N; i++) {
        pi += 4 * (pow(-1, i) / (2 * i + 1));
    }
    return pi;
}



int main()
{
    int N = 30000;
    double pi = 0;

    pi = calculate_pi(N);
    

    //TEST
    if (fabs(pi - M_PI) > 0.01) {
        fprintf(stderr, "Test failed");
    }
    else {
        fprintf(stderr, "TEST OK");
    }

    return 0;
}



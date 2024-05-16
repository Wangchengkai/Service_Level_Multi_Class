#pragma once
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>

#include "data.h"

using namespace std;

/* Definitions de constants */
/*==========================*/
//#define EPSILON 1.0e-5
//#define INFINI 1.0e10
#define period LengthOfOnePeriod
#define days 500000
#define qmax 100#define tao 0.5
#define N 80#define M 100
/* Definitions de fonctions */
/*==========================*/
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MOD(x) (x - 10000 * (x / 10000))
#define ALEA(x, m, M) (x = (int)(((double)rand() / (RAND_MAX)) * ((M) - (m) + 1)) + (m))
#define SQR(x) ((x) * (x))
#define Random ((double)rand() / (RAND_MAX))
#define EXP(mean) (-log(1.0 - ((double)rand() / (RAND_MAX))) * mean)

//group two
#define T NumberOfPeriodsCondsidered
#define mu Miu
//int p[T] = { 0 };
//double lamda[T] = {0};

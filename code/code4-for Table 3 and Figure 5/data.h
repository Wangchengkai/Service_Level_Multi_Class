#pragma once
#ifndef _DATA_
#define _DATA_

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
#include <algorithm>
#include<unordered_map>
#include<unordered_set>

using namespace std;



/* Definitions de constants */
/*==========================*/
#define NumberOfPriors 3

#define Qiangzhishoulian 200

#define INFINI 1.0e10
#define LengthOfOnePeriod 1.0
#define NumberOfPeriodsCondsidered 24
#define IsUsingMyselfCoding true
#define IsUsingEigen false
#define IsUsingSimu false
#define Fanzhen_1_OrYouhua_0 false//1 evolution  0 solve
#define NumberOfPiecesInAppro 40
#define ApproIndex 2


#define PunishForOneMoreServer 100.01//100.01

#define GTT 0.5
const double ThresholdOfServiceLevelList[NumberOfPeriodsCondsidered + 1] = {
	0.0,
	GTT,GTT,GTT,GTT,
	GTT,GTT,GTT,GTT,//5-8
	//GTT,GTT,GTT,GTT,
	//GTT,GTT,GTT,GTT,//13-16
	//GTT,GTT,GTT,GTT,
	//GTT,GTT,GTT,GTT //21-24
};

extern double Miu;
extern double MiuForEachPeriod[];


#define StartTimeOfNightShift 16
#define EndTimeOfNightShift 24
#define maxShiftTime 9
#define minShiftTime 1
#define TotalNumberOfServers 40

extern double torWaitingTime;
extern double torWaitingTimePL;


const double tor_SeveralCustomers[NumberOfPriors + 1] = { 0,  10.0 / 60.0,  0.5,  50.0 / 60.0 };

const double thresholdOfServiceLevel_SeveralCustomers[NumberOfPriors + 1] = { 0,0.95,0.8,0.7 };
const double thresholdOfServiceLevelTooHigh_SeveralCustomers[NumberOfPriors + 1] = { 0,0.99,0.90,0.85 };
//#define ThresholdOfServiceLevel 0.5
#define ThresholdOfServiceLevelTooHigh 0.9

#define MaxNumOfShifts 300
#define ExpendOperBound 1
#define PrecisionOfPoissonParameter 10000
#define PrecisionOfIntegralOfTInOnePeriod 100
#define MaxNumberOfTranstionTimesInUniform 300
#define MaxValueOfStatesInUniform 300
#define MaxValueOfStatesInIntegral MaxValueOfStatesInUniform
#define MaxNumberOfTranstionTimesInIntegral MaxNumberOfTranstionTimesInUniform
#define BreakThreshold 1.0e-9
#define EPSILON 1.0e-5
#define MaxNumOfK 500
#define MaxNumOfState  MaxValueOfStatesInUniform
#define MaxNumOfPeriod NumberOfPeriodsCondsidered
//#define EPSILON 1.0e-5
//#define LengthOfOnePeriod 1.0
//#define PrecisionOfIntegralOfTInOnePeriod 200
//*****************************************


#define QMAX1 30
#define QMAX2 150
#define QMAX3 300
#define SMAX 16
#define QMAXSUM QMAX1+QMAX2+QMAX3
#define QMAXSUMS QMAX1+QMAX2+SMAX



const string thiscase = " ";
extern double Lambda_initial[][NumberOfPeriodsCondsidered + 1];
extern double Lamda[];
extern double LamdaMP[];
extern double LamdaLP[];
extern double partA;
extern double partB;
extern double partC;

const double Lambda_initial_initial[NumberOfPeriodsCondsidered + 1] = { 0,

//9.0000 ,10.0000 ,11.0000 ,12.0000 ,13.0000 ,14.0000 ,15.0000 ,16.0000 ,17.0000 ,18.0000 ,19.0000 ,20.0000 ,21.0000 ,22.0000 ,23.0000 ,24.0000 ,1.0000 ,2.0000 ,3.0000 ,4.0000 ,5.0000 ,6.0000 ,7.0000 ,8.0000 ,
//23.0769 ,12.7692 ,16.3846 ,17.9231 ,14.7692 ,10.8462 ,12.6923 ,14.3077 ,14.3077 ,13.3846 ,11.7692 ,12.3077 ,21.1538 ,17.5385 ,13.9231 ,8.9231 ,6.6154 ,4.6154 ,3.3077 ,2.7692 ,2.6923 ,3.1538 ,3.4615 ,9.4615 ,
//20.9286 ,13.7143 ,16.1429 ,15.0714 ,14.1429 ,10.9286 ,11.0714 ,13.9286 ,12.5000 ,13.0714 ,12.2143 ,14.0000 ,17.8571 ,17.2857 ,13.4286 ,9.4286 ,6.4286 ,5.5714 ,4.0000 ,3.0000 ,2.8571 ,2.5714 ,3.6429 ,9.6429 ,
//22.3077 ,14.0000 ,13.9231 ,14.7692 ,13.0769 ,9.6154 ,12.2308 ,11.5385 ,12.0000 ,11.4615 ,11.6154 ,12.1538 ,19.0769 ,17.6923 ,13.1538 ,8.7692 ,6.8462 ,4.1538 ,3.4615 ,3.6154 ,2.6154 ,2.2308 ,4.8462 ,10.0769 ,
//22.8462 ,13.0000 ,16.2308 ,14.3846 ,11.3846 ,10.5385 ,11.9231 ,12.0769 ,13.0000 ,13.9231 ,10.2308 ,11.6923 ,18.2308 ,16.3077 ,14.4615 ,9.3846 ,6.5385 ,6.3077 ,3.4615 ,3.6154 ,2.5385 ,2.1538 ,2.3077 ,9.3846 ,
//21.9231 ,12.3846 ,14.0769 ,14.8462 ,13.3077 ,8.9231 ,12.8462 ,12.3846 ,13.0000 ,12.0000 ,12.1538 ,10.6923 ,19.2308 ,14.3077 ,13.5385 ,10.3846 ,6.1538 ,4.4615 ,2.3077 ,2.9231 ,3.7692 ,3.1538 ,3.6923 ,8.2308 ,
//18.3846 ,12.4615 ,11.3077 ,15.2308 ,12.5385 ,9.9231 ,8.7692 ,11.3846 ,12.2308 ,13.2308 ,12.0769 ,11.4615 ,17.6154 ,10.9231 ,12.0000 ,8.5385 ,6.7692 ,4.8462 ,3.2308 ,3.0769 ,3.3846 ,1.8462 ,4.7692 ,9.3846 ,
//23.0769 ,16.2308 ,19.1538 ,17.7692 ,12.4615 ,12.6154 ,18.5385 ,14.2308 ,13.0000 ,11.6154 ,9.6923 ,12.3077 ,16.9231 ,15.9231 ,11.8462 ,8.7692 ,5.7692 ,5.2308 ,2.8462 ,3.3077 ,2.8462 ,2.7692 ,4.3077 ,9.4615 ,
//
//23.0000 ,14.0000 ,19.2308 ,21.4615 ,16.3846 ,12.6923 ,14.6154 ,14.0000 ,11.8462 ,14.6154 ,10.5385 ,15.1538 ,19.9231 ,17.9231 ,12.8462 ,9.8462 ,6.3077 ,5.5385 ,3.5385 ,2.8462 ,2.7692 ,2.3077 ,4.6923 ,8.6154 ,
//22.2500 ,15.1667 ,16.5000 ,19.2500 ,15.0833 ,11.7500 ,13.6667 ,13.5833 ,12.2500 ,16.0833 ,13.7500 ,12.5000 ,21.6667 ,19.0000 ,13.1667 ,10.0833 ,6.6667 ,6.4167 ,4.3333 ,2.6667 ,2.3333 ,2.1667 ,5.0000 ,7.7500 ,
//22.9167 ,13.7500 ,16.8333 ,17.5833 ,14.1667 ,12.2500 ,11.1667 ,13.0000 ,13.2500 ,14.3333 ,13.8333 ,13.8333 ,16.5833 ,18.9167 ,13.0000 ,9.5000 ,8.0833 ,4.9167 ,4.0000 ,2.6667 ,1.7500 ,3.0000 ,7.0000 ,6.5833 ,
//23.3077 ,13.2308 ,17.9231 ,17.4615 ,14.2308 ,11.1538 ,12.4615 ,15.1538 ,12.0000 ,16.1538 ,11.6154 ,12.1538 ,18.0000 ,17.3846 ,12.6923 ,8.3077 ,5.9231 ,6.0000 ,2.9231 ,2.5385 ,2.0000 ,2.7692 ,4.0769 ,9.0000 ,
//26.7143 ,15.5000 ,17.2143 ,18.5714 ,14.5000 ,11.7143 ,11.0000 ,14.2857 ,10.9286 ,13.7143 ,12.0714 ,13.4286 ,15.8571 ,18.7857 ,14.1429 ,8.8571 ,7.3571 ,5.2857 ,5.0000 ,2.7857 ,2.5714 ,2.4286 ,4.9286 ,7.0000 ,
//23.7143 ,14.5714 ,20.2143 ,14.4286 ,15.0714 ,14.0714 ,14.4286 ,15.0000 ,14.0000 ,12.2143 ,12.7857 ,15.8571 ,17.2143 ,15.8571 ,14.3571 ,9.6429 ,7.5714 ,5.2143 ,4.9286 ,3.0000 ,2.7143 ,2.7143 ,4.7857 ,7.2857 ,
//28.0000 ,17.2308 ,21.9231 ,22.6154 ,17.3846 ,14.9231 ,17.0000 ,19.0769 ,14.6154 ,14.6154 ,10.0000 ,14.2308 ,17.4615 ,13.6923 ,12.6154 ,8.3846 ,5.4615 ,6.2308 ,3.9231 ,2.6923 ,2.5385 ,2.3846 ,5.6923 ,9.2308 ,
//
//
//11.7692 ,12.6154 ,17.1538 ,15.9231 ,13.8462 ,13.0769 ,10.3846 ,14.5385 ,10.2308 ,12.6154 ,13.7692 ,13.3846 ,20.8462 ,17.6154 ,15.2308 ,9.9231 ,6.8462 ,6.6154 ,5.0000 ,2.7692 ,2.7692 ,2.5385 ,3.6154 ,7.3077 ,
//10.7692 ,10.3846 ,13.8462 ,14.9231 ,11.9231 ,11.2308 ,10.3077 ,13.4615 ,12.0769 ,14.3846 ,11.5385 ,12.9231 ,20.9231 ,15.5385 ,12.6154 ,11.6154 ,7.6154 ,6.6923 ,4.5385 ,4.0000 ,4.6154 ,2.6923 ,5.2308 ,8.5385 ,
//11.8462 ,15.0769 ,17.3077 ,14.8462 ,12.6154 ,11.2308 ,17.3846 ,14.9231 ,13.3077 ,11.6154 ,11.3077 ,12.9231 ,15.6923 ,13.9231 ,11.5385 ,12.4615 ,7.0769 ,5.4615 ,3.4615 ,4.6923 ,2.6923 ,2.5385 ,5.2308 ,8.4615 ,

//24.48603866,13.83246395,15.7226143,16.58177355,14.86345505,9.966247315,14.3479595,13.83246395,14.51979135,13.40288432,13.57471617,11.94231359,17.47898128,14.98036207,13.12120282,11.59864989,6.87327401,4.983123658,2.577477754,3.264805155,4.209880331,3.52255293,4.123964406,9.193003989
//22.8462 ,13.0000 ,16.2308 ,14.3846 ,11.3846 ,10.5385 ,11.9231 ,12.0769 ,13.0000 ,13.9231 ,10.2308 ,11.6923 ,18.2308 ,16.3077 ,14.4615 ,9.3846 ,6.5385 ,6.3077 ,3.4615 ,3.6154 ,2.5385 ,2.1538 ,2.3077 ,9.3846 ,
//22.2500 ,15.1667 ,16.5000 ,19.2500 ,15.0833 ,11.7500 ,13.6667 ,13.5833 ,12.2500 ,16.0833 ,13.7500 ,12.5000 ,21.6667 ,19.0000 ,13.1667 ,10.0833 ,6.6667 ,6.4167 ,4.3333 ,2.6667 ,2.3333 ,2.1667 ,5.0000 ,7.7500 ,


//2.769926512,3.244771057,3.561334087,9.734313171,

//6.806105144,4.748445449,3.403052572,2.84906727,2.769926512,3.244771057,3.561334087,9.734313171,23.74222725,13.13736574,16.85698135,18.4397965,15.19502544,11.15884681,13.05822499,14.72018089,14.72018089,13.7704918,12.1085359,12.6625212,17.76370831,15.04409271,14.32447711,9.180327869

//12,12,12,12,12
//16,16,16,16,16
//12,12,12,0,0,0,0,0
//12,12

23.74222725,13.13736574,16.85698135,18.4397965,15.19502544,11.15884681,13.05822499,14.72018089,14.72018089,13.7704918,12.1085359,12.6625212,17.76370831,15.04409271,14.32447711,9.180327869,6.806105144,4.748445449,3.403052572,2.84906727,2.769926512,3.244771057,3.561334087,9.734313171
//22.24511931,14.57700651,17.15835141,16.01952278,15.03253796,11.61605206,11.76789588,14.80477223,13.28633406,13.89370933,12.98264642,14.88069414,15.98047722,15.37310195,13.27331887,10.02169197,6.8329718,5.921908894,4.251626898,3.188720174,3.036876356,2.73318872,3.872017354,10.2494577
//26.7143,15.5,17.2143,18.5714,14.5,11.7143,11,14.2857,10.9286,13.7143,12.0714,13.4286,15.8571,18.7857,14.1429,8.8571,7.3571,5.2857,5,2.7857,2.5714,2.4286,4.9286,7
//23,14,19.2308,21.4615,16.3846,12.6923,14.6154,14,11.8462,14.6154,10.5385,15.1538,19.9231,17.9231,12.8462,9.8462,    6.3077,5.5385,3.5385,2.8462,2.7692,2.3077,4.6923,8.6154


};

//important: initial state
const int insertOfClass1 = 0;
const int insertOfClass2 = 0;
const int insertOfClass3 = 0;


class OperData {
public:
	static void initialInitialData() {
		cout << "test djalhfkl\n";
		for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			Lambda_initial[1][indexOfPeriod] = Lambda_initial_initial[indexOfPeriod] * partA;
			Lambda_initial[2][indexOfPeriod] = Lambda_initial_initial[indexOfPeriod] * partB;
			Lambda_initial[3][indexOfPeriod] = Lambda_initial_initial[indexOfPeriod] * partC;
		}
	}

	static void initialData(int switchIndex) {
		switch (switchIndex)
		{
		case 1:
			std::cout << "A->B" << endl;
			torWaitingTime = tor_SeveralCustomers[1];
			//torWaitingTimePL = tor_SeveralCustomers[2];
			for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
				LamdaLP[indexOfPeriod] = Lambda_initial[2][indexOfPeriod];
			}
			break;
		case 2:
			std::cout << "A->B--C" << endl;
			//torWaitingTime = tor_SeveralCustomers[1];
			//torWaitingTimePL = tor_SeveralCustomers[2];
			for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
				LamdaMP[indexOfPeriod] = Lambda_initial[2][indexOfPeriod];
				LamdaLP[indexOfPeriod] = Lambda_initial[3][indexOfPeriod];
			}
			break;
		case 3:

		default:
			break;
		}
	}
};


const double dt = 1;




#endif 
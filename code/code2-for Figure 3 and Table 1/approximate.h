#pragma once
#ifndef _APPROXIMATE_
#define _APPROXIMATE_

#include "data.h"
#include "operFunctions.h"

class HistoryOfVectorOfUniformState;

class APPROXIMATE {
	Solve1 solveForApproximate;
	int offShift[NumberOfPeriodsCondsidered + 1] = { 0 };
public:
	double listOfTao;

	double readfileOfSkellam();

	double valueOfPai_Z_T(int z, double t);
	double valueOfA(double t, double tao);
	double valueOfAarr(double t, double tao);

	double integralPai(int z, double LB, double UB);
	double integralPaiMultA(int z, double LB, double UB, double tao);
	double integralPaiMultA_arr(int z, double LB, double UB, double tao);



	double expectedServiceRate_Abar(int l, int z, double LB, double UB, double tao);
	double expectedArrivalRate_Abar_arr(int l, int z, double LB, double UB, double tao);

	double skellam(double serviceRate, double arrivalRate, int value);
	double skellam_dealwith0(int indexOfArrivalRate, int indexOfServiceRate, double arrivalRate, double serviceRate, int value);

	double valueOfBeta2_Left(int l, int z, double tao);
	double valueOfBeta2_Right(int l, int z, double tao);
	double valueOfBeta2_new(int l, int z, double tao,double left,double right);
	double valueOfBeta2_Choose(int l, int z, double tao, double t);

	double valueOfBeta1_Left(int l, int z, double maxOfTao);
	double valueOfBeta1_Right(int l, int z, double maxOfTao);
	double valueOfBeta1_Choose(int l, int z, double maxOfTao, double t);

	double foreachZ(int l, int z, double maxOfTao);

	double serviceLevel(int l, Shifts shiftEnv, Solve1 solve, double maxOfTao, HistoryOfVectorOfUniformState p_iq);

	int start();

	void clearMap();

};






#endif
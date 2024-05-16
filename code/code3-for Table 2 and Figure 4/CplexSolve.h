#pragma once
//#include "data.h"
#include <ilcplex/ilocplex.h>
#include "OperFunctions.h"

class CplexSolve {
	//Solve1 solve;
	bool isShiftIncludePeriod[MaxNumOfShifts + 1][NumberOfPeriodsCondsidered + 1] = { 0 };
	int totalNumOfShifts;
	double totalLamda();
	//int solutionOnShifts[MaxNumOfShifts + 1] = { 0 };
public:
	void initial_AfterShiftEnvSetting();
	void initialTotalNumOfShifts();
	void initialIsShiftIncludePeriod();
	SolutionOnShift cplexModel(Solve1 solve1);
};




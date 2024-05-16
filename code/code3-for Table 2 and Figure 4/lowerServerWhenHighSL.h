#pragma once
#ifndef _LOWERSERVERWHENHIGHSL_
#define _LOWERSERVERWHENHIGHSL_
//#include "data.h"
#include "OperFunctions.h"

class LowerServerWhenHighSL {	
public:
	double computeSumOfServiceTime(SolutionOfServers solutionOfServers);

	void doIfOneServiceLeverOverHighThreshold(int indexOfPeriod, Solve1* solve, Shifts shiftsEnv);

private:
	LowBoundOfServers tryLowBoundMinusMinus(int indexOfPeriod,LowBoundOfServers lowBoundOfServers);

	//LB->cplex
	//cplex solution->st->SL

	bool isSaveTheChange(int indexOfCurrentPeriod, Solve1 tempSlove, double lastObjValue);

	void saveNewLowBound(Solve1* solve, Solve1 tempSolve);








};



























#endif // !1
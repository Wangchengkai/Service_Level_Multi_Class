#include "lowerServerWhenHighSL.h"
#include "CplexSolve.h"

SolutionOfServers convertFromShiftToServers(SolutionOnShift solutionOnShift, Shifts shiftEnv);

double LowerServerWhenHighSL::computeSumOfServiceTime(SolutionOfServers solutionOfServers) {
	double sum = 0;
	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
		sum += solutionOfServers.getSolutionOfServersForTime(indexOfPeriod) * LengthOfOnePeriod;
	return sum;
}

void LowerServerWhenHighSL::doIfOneServiceLeverOverHighThreshold(int indexOfPeriod, Solve1 *solve,Shifts shiftsEnv) {
	double lastObjValue = computeSumOfServiceTime((*solve).currentSolutionOfServers);
	Solve1 tempSolve = (*solve);
	if (indexOfPeriod >= StartTimeOfNightShift && indexOfPeriod <= EndTimeOfNightShift && (double(rand())/RAND_MAX < 0.2)) {
		for (int nightPeriod = StartTimeOfNightShift+1; nightPeriod <= EndTimeOfNightShift; nightPeriod++) {
			tempSolve.lowBoundOfServers = tryLowBoundMinusMinus(nightPeriod, tempSolve.lowBoundOfServers);
		}
	}
	else tempSolve.lowBoundOfServers = tryLowBoundMinusMinus(indexOfPeriod,tempSolve.lowBoundOfServers);
	
	CplexSolve tempCplexSolve;
	tempCplexSolve.initial_AfterShiftEnvSetting();	tempSolve.currentSolutionOnShift = tempCplexSolve.cplexModel(tempSolve);	tempSolve.currentSolutionOfServers = convertFromShiftToServers(tempSolve.currentSolutionOnShift, shiftsEnv);

	tempSolve.computeServiceLevel_SeveralType();
		
	tempSolve.periodUnits.setPeriodUnits((*solve).serviceLevelForSeveralCustomers);

	if (isSaveTheChange(indexOfPeriod, tempSolve, lastObjValue))
	{
		cout << "save" << endl;
		saveNewLowBound(solve, tempSolve);
	}
	else { cout << "not save" << endl;
}
}

LowBoundOfServers LowerServerWhenHighSL::tryLowBoundMinusMinus(int indexOfPeriod,LowBoundOfServers lowBoundOfServers) {
	lowBoundOfServers.setLowBoundOfServersForTime(indexOfPeriod,
		max(lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriod) - 1,0));
	return lowBoundOfServers;
}

//cplex solution->st->SL

bool LowerServerWhenHighSL::isSaveTheChange(int indexOfCurrentPeriod, Solve1 tempSlove, double lastObjValue) {
	bool isSLBothOverLowBound = true;

	for (int indexofperiod = max(1, indexOfCurrentPeriod - 1);
		indexofperiod <= min(NumberOfPeriodsCondsidered, indexOfCurrentPeriod + 1); indexofperiod++)
		if (tempSlove.serviceLevelForSeveralCustomers[1].serviceLevelLessThanThreshold(indexofperiod,thresholdOfServiceLevel_SeveralCustomers[1])  ||
			tempSlove.serviceLevelForSeveralCustomers[2].serviceLevelLessThanThreshold(indexofperiod, thresholdOfServiceLevel_SeveralCustomers[2]) ||
			tempSlove.serviceLevelForSeveralCustomers[3].serviceLevelLessThanThreshold(indexofperiod, thresholdOfServiceLevel_SeveralCustomers[3])
			)isSLBothOverLowBound = false;


	if (indexOfCurrentPeriod >= StartTimeOfNightShift && indexOfCurrentPeriod <= EndTimeOfNightShift) {
		isSLBothOverLowBound = true;
		for (int indexofperiod = StartTimeOfNightShift + 1;
			indexofperiod <= EndTimeOfNightShift; indexofperiod++)
			if (tempSlove.serviceLevelForSeveralCustomers[1].serviceLevelLessThanThreshold(indexofperiod, thresholdOfServiceLevel_SeveralCustomers[1]) ||
				tempSlove.serviceLevelForSeveralCustomers[2].serviceLevelLessThanThreshold(indexofperiod, thresholdOfServiceLevel_SeveralCustomers[2]) ||
				tempSlove.serviceLevelForSeveralCustomers[3].serviceLevelLessThanThreshold(indexofperiod, thresholdOfServiceLevel_SeveralCustomers[3]) 
				)isSLBothOverLowBound = false;
	}

	if (isSLBothOverLowBound && (this->computeSumOfServiceTime(tempSlove.currentSolutionOfServers) <= lastObjValue))
		return true;
	else
		return false;
}

void LowerServerWhenHighSL::saveNewLowBound(Solve1 *solve, Solve1 tempSolve) {
	(*solve) = tempSolve;
}





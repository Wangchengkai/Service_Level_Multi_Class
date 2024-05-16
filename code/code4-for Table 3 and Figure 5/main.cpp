#include "showSomeItems.h"
#include "lowerServerWhenHighSL.h"
#include "computeServiceLevel_Uniformization.h"
#include <time.h>
#include "approximate.h"
#include"simu_mahazhou.h"
#define onedayDivided 8

clock_t startTimeClock[10] = { 0 }, endTimeClock[10] = { 0 };
double timeForCplex = 0;
double timeForUniformAndComputeSL = 0;
double timeForUniformAndComputeSLNew = 0;

SolutionOfServers convertFromShiftToServers(SolutionOnShift solutionOnShift, Shifts shiftEnv) {
	SolutionOfServers solutionOfServers;
	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
		solutionOfServers.setSolutionOfServersForTime(indexOfPeriods,
			solutionOnShift.getSumOfServersOnPeriod(indexOfPeriods, shiftEnv));
	}
	return solutionOfServers;
}

Shifts shiftsEnv;
CplexSolve cplexSolve;

//unordered_map<string, ServiceLevel[NumberOfPriors + 1]> LowBoundVisited;
//unordered_set<string>LowBoundVisitedSet;
unordered_map<string, ServiceLevel[NumberOfPriors + 1]> SoutionOnShiftVisited;
unordered_set<string>SoutionOnShiftVisitedSet;

int countUniformTimes = 0;

void Solve1::computeServiceLevel_SeveralType() {
	//if (LowBoundVisited.count(this->lowBoundOfServers.convertLowBoundToString()) >= 1) {
	//	this->serviceLevelForSeveralCustomers[1] = LowBoundVisited[this->lowBoundOfServers.convertLowBoundToString()][1];
	//	this->serviceLevelForSeveralCustomers[2] = LowBoundVisited[this->lowBoundOfServers.convertLowBoundToString()][2];
	//	this->serviceLevelForSeveralCustomers[3] = LowBoundVisited[this->lowBoundOfServers.convertLowBoundToString()][3];
	//}
	if (SoutionOnShiftVisitedSet.count(this->currentSolutionOnShift.convertSolutionOnShiftToString()) >= 1) 
	{
		cout << "key" << this->currentSolutionOnShift.convertSolutionOnShiftToString() << endl;
		this->serviceLevelForSeveralCustomers[1] = SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][1];
		this->serviceLevelForSeveralCustomers[2] = SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][2];
		this->serviceLevelForSeveralCustomers[3] = SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][3];
		showServiceLevelByPeriod_SeveralType(this->serviceLevelForSeveralCustomers);
	}
	else {
		countUniformTimes++;
		if (IsUsingSimu == true) {
			OperSimu operSimu;
			operSimu.sim_main(shiftsEnv, this->currentSolutionOnShift, this->serviceLevelForSeveralCustomers[1], this->serviceLevelForSeveralCustomers[2], this->serviceLevelForSeveralCustomers[3]);

			showServiceLevelByPeriod_SeveralType(this->serviceLevelForSeveralCustomers);
			//system("pause");
		}
		else{
			for (int switchIndex = 1; switchIndex <= NumberOfPriors - 1; switchIndex++)
			{
				OperData::initialInitialData();
				OperData::initialData(switchIndex);

				ComputationOfServicelLevel computationOfServicelLevel;

				switch (switchIndex)
				{
				case 1:
					computationOfServicelLevel = this->serviceLevel.initialComputeServiceLevel_TwoTypes(shiftsEnv, this->currentSolutionOfServers, this->currentSolutionOnShift);
															this->serviceLevel.computeServiceLevel_forHighPriority(shiftsEnv, this->currentSolutionOfServers, this->currentSolutionOnShift, computationOfServicelLevel);
					//showServiceLevelByPeriod( this->serviceLevel);
					this->serviceLevelForSeveralCustomers[1] = this->serviceLevel;
					break;

				case 2:
					//torWaitingTime = tor_SeveralCustomers[1];
					computationOfServicelLevel = this->serviceLevel.initialComputeServiceLevel_ThreeTypes(shiftsEnv, this->currentSolutionOfServers, this->currentSolutionOnShift);
					
										for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
						Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
					}
					torWaitingTimePL = tor_SeveralCustomers[2];
					this->serviceLevelForLowPriorityCustomer.computeServiceLevel_forLowPriority(shiftsEnv, this->currentSolutionOfServers, this->currentSolutionOnShift, computationOfServicelLevel, "Mid");
					//showServiceLevelByPeriod( this->serviceLevelForLowPriorityCustomer);
					this->serviceLevelForSeveralCustomers[2] = this->serviceLevelForLowPriorityCustomer;

										for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
						Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod] + Lambda_initial[2][indexOfPeriod];
					}
					torWaitingTimePL = tor_SeveralCustomers[3];
					this->serviceLevelForLowPriorityCustomer.computeServiceLevel_forLowPriority(shiftsEnv, this->currentSolutionOfServers, this->currentSolutionOnShift, computationOfServicelLevel, "Low");
					//showServiceLevelByPeriod( this->serviceLevelForLowPriorityCustomer);
					this->serviceLevelForSeveralCustomers[3] = this->serviceLevelForLowPriorityCustomer;


					break;
				default:
					break;
				}


			}
		}
		SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][1] = this->serviceLevelForSeveralCustomers[1];
		SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][2] = this->serviceLevelForSeveralCustomers[2];
		SoutionOnShiftVisited[this->currentSolutionOnShift.convertSolutionOnShiftToString()][3] = this->serviceLevelForSeveralCustomers[3];
		SoutionOnShiftVisitedSet.insert(this->currentSolutionOnShift.convertSolutionOnShiftToString());
		cout << "save key" << this->currentSolutionOnShift.convertSolutionOnShiftToString() << endl;
	}

}



class SolveFlow{
public://test
	Solve1 solve;
	
public:
	void initialEnvironment() {
		shiftsEnv.setAllShifts();
		//cplexSolve.initialIsShiftIncludePeriod();
		cplexSolve.initial_AfterShiftEnvSetting();
		
		//showShiftType(shiftsEnv);
		//system("pause");
			}

	void initialPart() {
		solve.lowBoundOfServers.setLowBoundOfServersAsOne();
		
		showLowBoundOfServers(solve.lowBoundOfServers);
		
		solve.currentSolutionOnShift = cplexSolve.cplexModel(solve);
		//showSolutionsByShift(solve.currentSolutionOnShift);
		
		solve.currentSolutionOfServers = convertFromShiftToServers(solve.currentSolutionOnShift, shiftsEnv);
		showSolutionsByPeriod(solve.currentSolutionOfServers);
		//showSolutionsByShift(solutionOnShift);
		//showSolutionsByPeriod(solve.currentSolutionOfServers);
		

	}

	void setLowBoundAfterMinusServers() {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			solve.lowBoundOfServers.setLowBoundOfServersForTime(indexOfPeriod,
				min(solve.currentSolutionOfServers.getSolutionOfServersForTime(indexOfPeriod),
					solve.lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriod)));
	}

	void mainPart() {
		cout << "***********************start One Iteration*******************" << endl;

		
		solve.computeServiceLevel_SeveralType();
				//showServiceLevelByPeriod(solve.serviceLevel);



		solve.periodUnits.setPeriodUnits(solve.serviceLevelForSeveralCustomers);
		//showPeriodsUnitsLessThanThreshold(solve.periodUnits);
		//cout << "test4" << endl;

		//cout << "old LB£º" << endl;
		showLowBoundOfServers(solve.lowBoundOfServers);
		
		//system("pause");

				for (int indexOfUnit = solve.periodUnits.totalNumOfPeriodUnits; indexOfUnit >= 1; indexOfUnit--)			solve.forOneUnit(indexOfUnit);
		//setLowBoundAfterMinusServers();

		//cout << "new LB£º" << endl;
		showLowBoundOfServers(solve.lowBoundOfServers);

		solve.currentSolutionOnShift = cplexSolve.cplexModel(solve);		//showSolutionsByShift(solve.currentSolutionOnShift);
		solve.currentSolutionOfServers = convertFromShiftToServers(solve.currentSolutionOnShift, shiftsEnv);

		//showServiceLevelByPeriod(solve.serviceLevel);
		//showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);



		LowerServerWhenHighSL lowerServerWhenHighSL;
		//for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
		for (int indexOfPeriod = NumberOfPeriodsCondsidered; indexOfPeriod >= 1; indexOfPeriod--) {
			if (solve.serviceLevelForSeveralCustomers[1].serviceLevelMoreThanHighThreshold(indexOfPeriod,thresholdOfServiceLevelTooHigh_SeveralCustomers[1]) &&
				solve.serviceLevelForSeveralCustomers[2].serviceLevelMoreThanHighThreshold(indexOfPeriod, thresholdOfServiceLevelTooHigh_SeveralCustomers[2])&&
				solve.serviceLevelForSeveralCustomers[3].serviceLevelMoreThanHighThreshold(indexOfPeriod, thresholdOfServiceLevelTooHigh_SeveralCustomers[3])
				)
								//if (indexOfPeriod >= StartTimeOfNightShift && indexOfPeriod <= EndTimeOfNightShift) {
				//	lowerServerWhenHighSL.doIfOneServiceLeverOverHighThreshold(EndTimeOfNightShift, &solve, shiftsEnv);
				//	indexOfPeriod = EndTimeOfNightShift;
				//}
				//else {
				lowerServerWhenHighSL.doIfOneServiceLeverOverHighThreshold(indexOfPeriod, &solve, shiftsEnv);

//				}
		}
		
		//cout << "new1 LB£º" << endl;
		showLowBoundOfServers(solve.lowBoundOfServers);
		//showSolutionsByPeriod(solve.currentSolutionOfServers);


		//showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);
		//system("pause");
	}

	void show(){
		//showServiceLevelByPeriod(solve.serviceLevel);
		showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);
		showShiftType(shiftsEnv);
		showSolutionsByShift(solve.currentSolutionOnShift);
		showSolutionsByShift_ToPython_Draw(solve.currentSolutionOnShift, shiftsEnv);
		LowerServerWhenHighSL L;
		cout << "current solution" << L.computeSumOfServiceTime(solve.currentSolutionOfServers) << endl;
	}
};

#include "functionsForXJ.h"
void Uniformization(int offshift[], int p[], double SL[]);

int main_fanzhen()
//int main_testUniform() 
{
	startTimeClock[5] = clock();
	SolveFlow solveFlow;
	OperData::initialInitialData();

	showArrivalRate_initial();
	solveFlow.initialEnvironment();
	showShiftType(shiftsEnv);
	Solve1 solve;

		for (int indexOfShift = 1; indexOfShift <= shiftsEnv.totalNumOfShifts; indexOfShift++) {
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 3);
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 1 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 6)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 5 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 11)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 11 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 1 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 0);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 6)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 5 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 7 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 11 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 6)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 5 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 7 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 11 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		
		//
				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 5)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 1 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 7 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 3);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 10 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

	}
	solve.currentSolutionOfServers = convertFromShiftToServers(solve.currentSolutionOnShift, shiftsEnv);
	showSolutionsByPeriod(solve.currentSolutionOfServers);
	if (IsUsingSimu == true) {
		OperSimu operSimu;

		OperData::initialInitialData();
		showServiceRate();
		showArrivalRate_initial();

		cout << "test" << endl;
		operSimu.sim_main(shiftsEnv, solve.currentSolutionOnShift, solve.serviceLevelForSeveralCustomers[1], solve.serviceLevelForSeveralCustomers[2], solve.serviceLevelForSeveralCustomers[3]);

		showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);

		system("pause");
	}
	else {
		for (int switchIndex = 1; switchIndex <= NumberOfPriors - 1; switchIndex++)
		{
			//OperData::initialInitialData();
			OperData::initialData(switchIndex);
			
			showArrivalRate();

			ComputationOfServicelLevel computationOfServicelLevel;

			switch (switchIndex)
			{
			case 1:
				computationOfServicelLevel = solve.serviceLevel.initialComputeServiceLevel_TwoTypes(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift);
												solve.serviceLevel.computeServiceLevel_forHighPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel);
				//showServiceLevelByPeriod(solve.serviceLevel);
				solve.serviceLevelForSeveralCustomers[1] = solve.serviceLevel;
				break;

			case 2:
				//torWaitingTime = tor_SeveralCustomers[1];
				computationOfServicelLevel = solve.serviceLevel.initialComputeServiceLevel_ThreeTypes(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift);
				
								for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
					Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
				}
				torWaitingTimePL = tor_SeveralCustomers[2];
				solve.serviceLevelForLowPriorityCustomer.computeServiceLevel_forLowPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel, "Mid");
				//showServiceLevelByPeriod(solve.serviceLevelForLowPriorityCustomer);
				solve.serviceLevelForSeveralCustomers[2] = solve.serviceLevelForLowPriorityCustomer;

								for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
					Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod] + Lambda_initial[2][indexOfPeriod];
				}
				torWaitingTimePL = tor_SeveralCustomers[3];
				solve.serviceLevelForLowPriorityCustomer.computeServiceLevel_forLowPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel, "Low");
				//showServiceLevelByPeriod(solve.serviceLevelForLowPriorityCustomer);
				solve.serviceLevelForSeveralCustomers[3] = solve.serviceLevelForLowPriorityCustomer;


				break;
			default:
				break;
			}


		}

		cout << endl;
		//system("pause");

		endTimeClock[5] = clock();
		cout << "totaltime\t" << double(endTimeClock[5] - startTimeClock[5]) / double(CLOCKS_PER_SEC) << endl;
		cout << totaltimeOfHyper4d << endl;;
		showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);
		system("pause");
	}
	return 0;
}





int main_youhua()
{
	startTimeClock[0] = clock();
	OperData::initialInitialData();
	SolveFlow solveFlow;
	solveFlow.initialEnvironment();
	cout << "test1" << endl;
	showShiftType(shiftsEnv);

	//system("pause");

	solveFlow.initialPart();
	cout << "test2" << endl;

	int countOfIter = 0;
	do {
		solveFlow.mainPart();
		countOfIter++;
	} while (
		!solveFlow.solve.serviceLevelForSeveralCustomers[1].serviceLevelBothMoreThanThreshold(thresholdOfServiceLevel_SeveralCustomers[1]) ||
		!solveFlow.solve.serviceLevelForSeveralCustomers[2].serviceLevelBothMoreThanThreshold(thresholdOfServiceLevel_SeveralCustomers[2]) ||
		!solveFlow.solve.serviceLevelForSeveralCustomers[3].serviceLevelBothMoreThanThreshold(thresholdOfServiceLevel_SeveralCustomers[3])
		);
	cout << "iteration times" << countOfIter << endl;
	cout << "Uniformization times" << countUniformTimes << endl;


	endTimeClock[0] = clock();
	cout << "totalTime:\t" << double(endTimeClock[0] - startTimeClock[0]) / CLOCKS_PER_SEC << endl;
	cout << "timeForCplex:\t" << timeForCplex << endl;
	cout << "timeForUniformAndComputeSL:\t" << timeForUniformAndComputeSL << endl;
	cout << "timeForUniformAndComputeSLNew:\t" << timeForUniformAndComputeSLNew << endl;


	cout << "test3" << endl;
	//system("pause");


	solveFlow.show();

	cout << thiscase << endl;

	system("pause");
	return 0;

}

int main() {
	//By this part you can use the approximate method
	
	//APPROXIMATE approximate;

	//approximate.start();

	//system("pause");

	//return 0;

	//By this part you can use the exact method(with Fanzhen_1_OrYouhua_0 == false, thus enable the solve part)
	
	if (Fanzhen_1_OrYouhua_0)main_fanzhen(); //evolution
	else main_youhua();//solve


	system("pause");
	system("pause");
	system("pause");
	system("pause");

	return 0;
}
int LowBoundOfServers::getLowBoundOfServersForTime(int numOfPeriod) {
	testPeriodIndex(numOfPeriod);
	return lowBoundOfServers[numOfPeriod - 1];
}
void LowBoundOfServers::setLowBoundOfServersForTime(int numOfPeriod, int valueOfServers) {
	testPeriodIndex(numOfPeriod);
	lowBoundOfServers[numOfPeriod - 1] = valueOfServers;
}
void LowBoundOfServers::lowBoundMinusMinus(int numOfPeriod) {
	lowBoundOfServers[numOfPeriod - 1]--;
}
void LowBoundOfServers::lowBoundPlusPlus(int numOfPeriod) {
	lowBoundOfServers[numOfPeriod - 1]++;
}
void LowBoundOfServers::setLowBoundOfServersAsOne() {
	memset(lowBoundOfServers, 0, sizeof(lowBoundOfServers));
	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
		setLowBoundOfServersForTime(indexOfPeriod, 1);
	}
	
	
}

int SolutionOfServers::getSolutionOfServersForTime(int numOfPeriod) const {
	testPeriodIndex(numOfPeriod);
	return solutionOfServers[numOfPeriod - 1];
}
void SolutionOfServers::setSolutionOfServersForTime(int numOfPeriod, int valueOfServers) {
	testPeriodIndex(numOfPeriod);
	solutionOfServers[numOfPeriod - 1] = valueOfServers;
}

double ServiceLevel::getServiceLevelForPeriod(int numOfPeriod) {
	testPeriodIndex(numOfPeriod);
	return serviceLevel[numOfPeriod - 1];
}
void ServiceLevel::setServiceLevel(int numOfPeriod, double valueOfSL) {
	serviceLevel[numOfPeriod - 1] = valueOfSL;
}
bool ServiceLevel::serviceLevelLessThanThreshold(int numOfPeriod, double threshold) {
	if (this->getServiceLevelForPeriod(numOfPeriod) < threshold) return true;
	else return false;
}
bool ServiceLevel::serviceLevelMoreThanHighThreshold(int numOfPeriod,double thresholdOfSLTooHigh) {
	if (this->getServiceLevelForPeriod(numOfPeriod) > thresholdOfSLTooHigh) return true;
	else return false;
}
int computeOffShiftInPeriod(Shifts shiftEnv, SolutionOnShift solutionOnShift, int endTimeOfOffShift) {
	int countOfShifts = 0;
	for (int indexOfShift = 1; indexOfShift <= shiftEnv.totalNumOfShifts; indexOfShift++) {
		if (shiftEnv.shifts[indexOfShift].endTime == endTimeOfOffShift) {
			countOfShifts += solutionOnShift.getServerOnShift(indexOfShift);
		}
	}
	return countOfShifts;
}
void ServiceLevel::computeServiceLevel(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			this->setServiceLevel(indexOfPeriod, 0.0);
	}
	else {

		startTimeClock[3] = clock();
		ComputationOfServicelLevel computationOfServicelLevelInNewWay;
		computationOfServicelLevelInNewWay.initialComputationOfServicelLevel(solutionOnShift, solutionOfServers, shiftEnv);
		//for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
		//	this->setServiceLevel(indexOfPeriod,
		//		computationOfServicelLevelInNewWay.computeServiceLevelForPeriod_H(indexOfPeriod, solutionOfServers, shiftEnv));
		//
		//endTimeClock[3] = clock(); timeForUniformAndComputeSLNew += double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC);
		//cout << "2\t" << double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC) << endl;


		startTimeClock[4] = clock();
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			//cout << "test***\t" << indexOfPeriod << endl;
			this->setServiceLevel(indexOfPeriod,
				computationOfServicelLevelInNewWay.computeServiceLevelForPeriod_L(indexOfPeriod, solutionOfServers, shiftEnv));

		}
		endTimeClock[4] = clock(); timeForUniformAndComputeSLNew += double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC);
		//cout << "3\t" << double(double(endTimeClock[4] - startTimeClock[4]) / CLOCKS_PER_SEC) << endl;

		//showServiceLevelByPeriod(*this);
		//system("pause");
		
	}
}
ComputationOfServicelLevel ServiceLevel::initialComputeServiceLevel_TwoTypes(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
	}
	else {
		ComputationOfServicelLevel computationOfServicelLevelInNewWay;
		computationOfServicelLevelInNewWay.initialComputationOfServicelLevel(solutionOnShift, solutionOfServers, shiftEnv);
		return computationOfServicelLevelInNewWay;
	}
}
ComputationOfServicelLevel ServiceLevel::initialComputeServiceLevel_ThreeTypes(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
	}
	else {
		ComputationOfServicelLevel computationOfServicelLevelInNewWay;
		computationOfServicelLevelInNewWay.initialComputationOfServicelLevel_3Type(solutionOnShift, solutionOfServers, shiftEnv);
		return computationOfServicelLevelInNewWay;
	}
}



void ServiceLevel::computeServiceLevel_forHighPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift, ComputationOfServicelLevel computationOfServicelLevelInNewWay) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			this->setServiceLevel(indexOfPeriod, 0.0);
	}
	else {
		startTimeClock[3] = clock();
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			this->setServiceLevel(indexOfPeriod,
				computationOfServicelLevelInNewWay.computeServiceLevelForPeriod_H(indexOfPeriod, solutionOfServers, shiftEnv));
		
		endTimeClock[3] = clock(); timeForUniformAndComputeSLNew += double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC);
		//cout << "2\t" << double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC) << endl;
	}
}

void ServiceLevel::computeServiceLevel_forLowPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift, ComputationOfServicelLevel computationOfServicelLevelInNewWay) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			this->setServiceLevel(indexOfPeriod, 0.0);
	}
	else {
		startTimeClock[4] = clock();
//#pragma omp parallel 
//#pragma omp for
//#pragma omp parallel for schedule(dynamic, 10) 
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			//cout << "test***\t" << indexOfPeriod << endl;
			this->setServiceLevel(indexOfPeriod,
				computationOfServicelLevelInNewWay.computeServiceLevelForPeriod_L(indexOfPeriod, solutionOfServers, shiftEnv));

		}
		endTimeClock[4] = clock(); timeForUniformAndComputeSLNew += double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC);
		cout << "3\t" << double(double(endTimeClock[4] - startTimeClock[4]) / CLOCKS_PER_SEC) << endl;

		//showServiceLevelByPeriod(*this);
		//system("pause");

	}
}

void ServiceLevel::computeServiceLevel_forLowPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift, ComputationOfServicelLevel computationOfServicelLevelInNewWay,string type) {
	if (!isServiveCapcityEnough(solutionOfServers)) {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			this->setServiceLevel(indexOfPeriod, 0.0);
	}
	else {
		startTimeClock[4] = clock();
		//#pragma omp parallel 
		//#pragma omp for
		//#pragma omp parallel for schedule(dynamic, 10) 
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			//cout << "test***\t" << indexOfPeriod << endl;
			this->setServiceLevel(indexOfPeriod,
				computationOfServicelLevelInNewWay.computeServiceLevelForPeriod_L(indexOfPeriod, solutionOfServers, shiftEnv, type));

		}
		endTimeClock[4] = clock(); timeForUniformAndComputeSLNew += double(double(endTimeClock[3] - startTimeClock[3]) / CLOCKS_PER_SEC);
		cout << "3\t" << double(double(endTimeClock[4] - startTimeClock[4]) / CLOCKS_PER_SEC) << endl;

		//showServiceLevelByPeriod(*this);
		//system("pause");

	}
}



bool ServiceLevel::isServiveCapcityEnough(SolutionOfServers solutionOfServers) {
	double sumOfCapcity = 0;
	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
		sumOfCapcity += double(solutionOfServers.getSolutionOfServersForTime(indexOfPeriod)) * double(MiuForEachPeriod[indexOfPeriod]);
	double totalLamda = 0;
	for (int period = 1; period <= NumberOfPeriodsCondsidered; period++)
		totalLamda += Lamda[period];	for (int period = 1; period <= NumberOfPeriodsCondsidered; period++)
		totalLamda += LamdaLP[period];	//cout << "sumOfCapcity\t" << sumOfCapcity << "\ttotalLamda\t" << totalLamda << endl;

	if (sumOfCapcity > totalLamda)return true;
	else return false;

}



void PeriodUnitsNotMeetSolutionLevel::setPeriodUnits(ServiceLevel *SL) {
	initialPeriodUnits();
	int countNumOfPeriodUnits = 0;
	int numofperiod = 1;
	while (numofperiod <= NumberOfPeriodsCondsidered) {
		if (SL[1].serviceLevelLessThanThreshold(numofperiod, thresholdOfServiceLevel_SeveralCustomers[1]) ||
			SL[2].serviceLevelLessThanThreshold(numofperiod, thresholdOfServiceLevel_SeveralCustomers[2]) ||
			SL[3].serviceLevelLessThanThreshold(numofperiod, thresholdOfServiceLevel_SeveralCustomers[3]))

			this->DoIfServiceLevelLessThanThreshold(&countNumOfPeriodUnits, &numofperiod, SL);
		else
			numofperiod++;
	}
	totalNumOfPeriodUnits = countNumOfPeriodUnits;

	
}
void PeriodUnitsNotMeetSolutionLevel::DoIfServiceLevelLessThanThreshold(int* countNumOfPeriodUnits, int* numofperiod,ServiceLevel *SL) {
	startTimeOfPeriodUnit[*countNumOfPeriodUnits] = *numofperiod - 1;
	while (*numofperiod <= NumberOfPeriodsCondsidered) {
		if (SL[1].serviceLevelLessThanThreshold(*numofperiod, thresholdOfServiceLevel_SeveralCustomers[1]) ||
			SL[2].serviceLevelLessThanThreshold(*numofperiod, thresholdOfServiceLevel_SeveralCustomers[2]) ||
			SL[3].serviceLevelLessThanThreshold(*numofperiod, thresholdOfServiceLevel_SeveralCustomers[3]))
			(*numofperiod)++;
		else
			break;
	}
	endTimeOfPeriodUnit[*countNumOfPeriodUnits] = *numofperiod - 1;
	(*countNumOfPeriodUnits)++;
	/*cout << "numofperiod\t" << *numofperiod << endl;
	cout << "countNumOfPeriodUnits\t" << *countNumOfPeriodUnits << endl;*/
}

void testIndexOfShiftOutOfRange(int index) {
	if (index > MaxNumOfShifts) { cout << "error2: indexofshift OutOfRange" << endl; system("pause"); }
}
void Shifts::setNightShift(int index) {
	int day = 1;
	shifts[index + day - 1].setShift(index+day-1, (StartTimeOfNightShift + (day - 1) * onedayDivided), (EndTimeOfNightShift + (day - 1) * onedayDivided)) ;}
void Shifts::setAllShifts() {
	int countIndex = 1;
	setNightShift(countIndex); countIndex+=1;
	int day = 1;
	for (int startTime = 0; startTime < StartTimeOfNightShift; startTime++)
		for (int endTime = startTime + minShiftTime; endTime <= min(StartTimeOfNightShift, startTime + maxShiftTime); endTime++) {
			//cout << "starttime\t" << startTime << "\tendtime\t" << endTime << endl;
			shifts[countIndex].setShift(countIndex, startTime, endTime);
			//cout << "shift\t" << shifts[countIndex].startTime << "\t" << shifts[countIndex].endTime << endl;
			//system("pause");
			countIndex++;
			}
	testIndexOfShiftOutOfRange(countIndex);
	this->totalNumOfShifts = countIndex - 1;}

bool Solve1::bothMoreThanSLThreshold(ServiceLevel *SL, int left, int right) {
	bool mark = true;
	for (int indexOfPeriod = left; indexOfPeriod <= right; indexOfPeriod++)
		if (SL[1].serviceLevelLessThanThreshold(indexOfPeriod, thresholdOfServiceLevel_SeveralCustomers[1]) ||
			SL[2].serviceLevelLessThanThreshold(indexOfPeriod, thresholdOfServiceLevel_SeveralCustomers[2]) ||
			SL[3].serviceLevelLessThanThreshold(indexOfPeriod, thresholdOfServiceLevel_SeveralCustomers[3])
		)mark = false;
	return mark;
}
void Solve1::forOneUnit(int indexOfUnit) {
	cout << "forOneUnit\t" << indexOfUnit << endl;

	int unitLeft = periodUnits.getStartTime(indexOfUnit) + 1;	int unitRight = periodUnits.getEndTime(indexOfUnit);
	updataLowBoundForOneUnitToPlusPlus(unitLeft, unitRight);
	updataServersNumForOneUnitByLowBound(unitLeft, unitRight);
	//serviceLevel.computeServiceLevel(shiftsEnv, this->currentSolutionOfServers,this->currentSolutionOnShift);
	this->computeServiceLevel_SeveralType();


	if (bothMoreThanSLThreshold(this->serviceLevelForSeveralCustomers, unitLeft, unitRight))
		tryToMinusLowBound(unitLeft, unitRight);
}
void Solve1::updataLowBoundForOneUnitToPlusPlus(int unitLeft, int unitRight) {
	
	for (int indexOfPeriod = unitLeft; indexOfPeriod <= unitRight; indexOfPeriod++) {
				//if (StartTimeOfNightShift <= indexOfPeriod && indexOfPeriod <= EndTimeOfNightShift) {
		//	for (int nightPeriod = StartTimeOfNightShift+1; nightPeriod <= EndTimeOfNightShift; nightPeriod++) {
		//		lowBoundOfServers.lowBoundPlusPlus(nightPeriod);
		//	}
		//	indexOfPeriod = EndTimeOfNightShift;
		//}
		//else	
			lowBoundOfServers.lowBoundPlusPlus(indexOfPeriod);
	}
}
void Solve1::updataServersNumForOneUnitByLowBound(int unitLeft, int unitRight) {
	for (int indexOfPeriod = unitLeft; indexOfPeriod <= unitRight; indexOfPeriod++) {
		currentSolutionOfServers.setSolutionOfServersForTime(indexOfPeriod,
			max(
				currentSolutionOfServers.getSolutionOfServersForTime(indexOfPeriod),
				lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriod)));
	}
}
SolutionOfServers Solve1::tryToMinusServers(int unitLeft, int unitRight) {
	Solve1 tempSolve = *this;

	//SolutionOfServers tempSolutionOfServers = currentSolutionOfServers;
	//ServiceLevel tempServicelLevel;
	int operLeft = max(unitLeft - ExpendOperBound, 1);
	int operRight = min(unitRight + ExpendOperBound, NumberOfPeriodsCondsidered);
	//to debug

	//PeriodOutOfRange POR;
		for (int indexToMinus = operLeft; indexToMinus <= operRight; indexToMinus++) {
		if (this->lowBoundOfServers.getLowBoundOfServersForTime(indexToMinus) <= 1)continue;
		
		SolutionOfServers slackForBackspace = tempSolve.currentSolutionOfServers;
				//if (tempSolutionOfServers.getSolutionOfServersForTime(indexToMinus) 
		//- this->lowBoundOfServers.getLowBoundOfServersForTime(indexToMinus)<=1)
						minusServerInOnePeriodTo_LBminusOne(&(tempSolve.currentSolutionOfServers), indexToMinus);
		
		tempSolve.computeServiceLevel_SeveralType();
		

				if (!bothMoreThanSLThreshold(tempSolve.serviceLevelForSeveralCustomers, operLeft, operRight))
			tempSolve.currentSolutionOfServers = slackForBackspace;

	}
	return tempSolve.currentSolutionOfServers;
}

bool Solve1::tryToMinusLowBound(int unitLeft, int unitRight) {
	Solve1 tempSolve = *this;
	bool isChange = false;
	int operLeft = max(unitLeft - ExpendOperBound, 1);
	int operRight = min(unitRight + ExpendOperBound, NumberOfPeriodsCondsidered);
	//to debug

	//PeriodOutOfRange POR;
	for (int indexToMinus = operLeft; indexToMinus <= operRight; indexToMinus++) {
		if (tempSolve.lowBoundOfServers.getLowBoundOfServersForTime(indexToMinus) <= 1)continue;
		

		Solve1 slackForBackspace = tempSolve;
		if (StartTimeOfNightShift <= indexToMinus && indexToMinus <= EndTimeOfNightShift) {
			for (int nightPeriod = StartTimeOfNightShift+1; nightPeriod <= EndTimeOfNightShift; nightPeriod++) {
				lowBoundOfServers.lowBoundMinusMinus(indexToMinus);
			}
			indexToMinus = EndTimeOfNightShift;
		}
		else
			tempSolve.lowBoundOfServers.lowBoundMinusMinus(indexToMinus);
		tempSolve.currentSolutionOnShift = cplexSolve.cplexModel(tempSolve);
		tempSolve.currentSolutionOfServers = convertFromShiftToServers(tempSolve.currentSolutionOnShift, shiftsEnv);
		tempSolve.computeServiceLevel_SeveralType();
		//tempSolve.serviceLevel.computeServiceLevel(shiftsEnv, tempSolve.currentSolutionOfServers, tempSolve.currentSolutionOnShift);
		
		
		if (!bothMoreThanSLThreshold(tempSolve.serviceLevelForSeveralCustomers, operLeft, operRight))
			tempSolve = slackForBackspace;
		else isChange = true;
	}
	if (isChange)*this = tempSolve;
	return isChange;

}


void Solve1::minusServerInOnePeriodTo_LBminusOne(SolutionOfServers* solutionOfMinusOneServer, int period) {
	(*solutionOfMinusOneServer).setSolutionOfServersForTime(period,
		(this->lowBoundOfServers.getLowBoundOfServersForTime(period) - 1));
}



void CplexSolve::initialTotalNumOfShifts() {
	totalNumOfShifts = shiftsEnv.totalNumOfShifts;
}
void CplexSolve::initialIsShiftIncludePeriod() {
	for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
		for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
			if (shiftsEnv.shifts[indexOfShifts].periodInShift(indexOfPeriods))
				isShiftIncludePeriod[indexOfShifts][indexOfPeriods] = true;
			else
				isShiftIncludePeriod[indexOfShifts][indexOfPeriods] = false;
		}
}

void CplexSolve::initial_AfterShiftEnvSetting() {
	initialTotalNumOfShifts();
	initialIsShiftIncludePeriod();
	//cout << "!!@" << endl;
	//cout << this->totalNumOfShifts << endl;
	//system("pause");
}



SolutionOnShift CplexSolve::cplexModel(Solve1 solve1)
{
	startTimeClock[1] = clock();

	SolutionOnShift solutionOnShift;
	//input
	ILOSTLBEGIN
		IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray varNumOfServerInShift(env, MaxNumOfShifts + 1, 0, TotalNumberOfServers, ILOINT);

		IloExpr obj(env);
		for (int indexOfShift = 1; indexOfShift <= totalNumOfShifts; indexOfShift++) {
			obj += varNumOfServerInShift[indexOfShift] * double(shiftsEnv.shifts[indexOfShift].length);
		}
		for (int indexOfShift = 1; indexOfShift <= totalNumOfShifts; indexOfShift++) {
			obj += varNumOfServerInShift[indexOfShift]* PunishForOneMoreServer;
		}
		model.add(IloMinimize(env, obj));
		//obj.end();
		
		for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
			IloExpr sumOfServers(env);
			for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
				sumOfServers += isShiftIncludePeriod[indexOfShifts][indexOfPeriods] * varNumOfServerInShift[indexOfShifts];
			model.add(sumOfServers >= solve1.lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriods));
		}

		IloExpr sumOfServers1(env);		for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) 
			for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
				sumOfServers1 += double(isShiftIncludePeriod[indexOfShifts][indexOfPeriods]) * varNumOfServerInShift[indexOfShifts] * MiuForEachPeriod[indexOfPeriods];
		model.add(sumOfServers1  >= totalLamda() + 0.1);
		
		IloExpr sumOfServers2(env);
		for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
			sumOfServers2 += varNumOfServerInShift[indexOfShifts];
		model.add(sumOfServers2 <= TotalNumberOfServers);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.solve();
		
		if (cplex.getStatus() == IloAlgorithm::Infeasible) {
			env.out() << "no solution" << endl;
			system("pause");
		}
		else {
			for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
				solutionOnShift.setServerOnShift(indexOfShifts,
					cplex.getValue(varNumOfServerInShift[indexOfShifts]));
						cout << "current cplex value\t" << cplex.getValue(obj) << endl;
		}
	}
	catch (IloException& e)
	{
		cerr << " ERROR: " << e << endl;
	}
	catch (...)
	{
		cerr << " ERROR\n";
	};
	env.end();
	endTimeClock[1] = clock(); timeForCplex += double(double(endTimeClock[1] - startTimeClock[1]) / CLOCKS_PER_SEC);
	return solutionOnShift;
}

double CplexSolve::totalLamda(){
	double totalLamda = 0;
	for(int indexOfPrior=1;indexOfPrior<= NumberOfPriors;indexOfPrior++)
		for (int period = 1; period <= NumberOfPeriodsCondsidered; period++)
			totalLamda += Lambda_initial[indexOfPrior][period];
	return totalLamda;
}


HistoryOfVectorOfUniformState historyOfVectorOfUniformStateInAppro[25];

#define MaxNumOfZ 40

double APPROXIMATE::serviceLevel(int l, Shifts shiftEnv, Solve1 solve, double maxOfTao, HistoryOfVectorOfUniformState P_iq) {
	//initial solve
	solveForApproximate = solve;
	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
		offShift[indexOfPeriod] = computeOffShiftInPeriod(shiftEnv, solveForApproximate.currentSolutionOnShift, indexOfPeriod);
	}
	historyOfVectorOfUniformStateInAppro[l] = P_iq;


	double sum = 0;
	for (int z = solve.currentSolutionOfServers.getSolutionOfServersForTime(l); z <= MaxNumOfZ; z++)	{
				if (integralPai(z, double(l - 1), double(l)) <= 1e-5 
			&& z >= solve.currentSolutionOfServers.getSolutionOfServersForTime(l) + 5) break;

		double tempdouble = foreachZ(l, z, maxOfTao);
		//cout << "l\t" << l << "\tz\t" << z <<"\tvalue\t"<<tempdouble<< endl;
		sum += tempdouble;// foreachZ(l, z, maxOfTao);
	}
	return 1.0 - sum;// / double(MaxNumOfZ);

}


unordered_map<string, double> valueOfPai_Z_T_MAP;

double APPROXIMATE::valueOfPai_Z_T(int z, double t)
{
	//find map
	string tempKey = to_string(z) + "inter" + to_string(t);
	if (valueOfPai_Z_T_MAP.count(tempKey)) return valueOfPai_Z_T_MAP[tempKey];


	int indexOfState = z;
	double timeInThisPeriod = t - int(t);
	int indexOfPeriod = int(t) + 1;


	double gama = double(solveForApproximate.currentSolutionOfServers.getSolutionOfServersForTime(indexOfPeriod)) * double(MiuForEachPeriod[indexOfPeriod]) + Lamda[indexOfPeriod] + LamdaLP[indexOfPeriod];
	double gamaMultiT = gama * timeInThisPeriod;
	double sumProbOfState = 0;
	for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInIntegral; indexOfTran++) {
		sumProbOfState += util.getPoisson(gamaMultiT, indexOfTran) * historyOfVectorOfUniformStateInAppro[int(t-(1e-5))+1].getProbOfTotalQueueLengthAfterTran(indexOfTran, indexOfState);
	}
	//save map
	valueOfPai_Z_T_MAP[tempKey] = sumProbOfState;

	return sumProbOfState;
}

int APPROXIMATE::start() {
	cout << "hello world" << endl;
	this->readfileOfSkellam();

	//cout << valueOfA(2.5, 2.5) << endl;

	SolveFlow solveFlow;
	OperData::initialInitialData();

	showArrivalRate_initial();
	solveFlow.initialEnvironment();
	showShiftType(shiftsEnv);
	Solve1 solve;

		for (int indexOfShift = 1; indexOfShift <= shiftsEnv.totalNumOfShifts; indexOfShift++) {
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 5)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 4 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 6 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 11)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 10 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 15)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 6)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 5 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 7 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 11 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);


				//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 4)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 6)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 7)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 3 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 5 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 10)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 7 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 14)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 9 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 11 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 12 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 1);
		//if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
		//	solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);





	}
	solve.currentSolutionOfServers = convertFromShiftToServers(solve.currentSolutionOnShift, shiftsEnv);
	showSolutionsByPeriod(solve.currentSolutionOfServers);

	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
		this->offShift[indexOfPeriod] = computeOffShiftInPeriod(shiftsEnv, solveForApproximate.currentSolutionOnShift, indexOfPeriod);
	}
	this->solveForApproximate = solve;
		HistoryOfVectorOfUniformState P_iq1;
	for (int switchIndex = 2; switchIndex <= NumberOfPriors - 1; switchIndex++)
	{
		//OperData::initialInitialData();
		OperData::initialData(switchIndex);

		showArrivalRate();

		ComputationOfServicelLevel computationOfServicelLevel;

		switch (switchIndex)
		{
		case 1:
			computationOfServicelLevel = solve.serviceLevel.initialComputeServiceLevel_TwoTypes(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift);
									solve.serviceLevel.computeServiceLevel_forHighPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel);
			//showServiceLevelByPeriod(solve.serviceLevel);
			solve.serviceLevelForSeveralCustomers[1] = solve.serviceLevel;
			break;

		case 2:
			//torWaitingTime = tor_SeveralCustomers[1];



			computationOfServicelLevel = solve.serviceLevel.initialComputeServiceLevel_ThreeTypes(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift);
			

			for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
			}
			torWaitingTimePL = tor_SeveralCustomers[2];
			this->solveForApproximate = solve;

			startTimeClock[4] = clock();
			//#pragma omp parallel for schedule(dynamic, 10) 
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				HistoryOfVectorOfUniformState P_iq;
				P_iq.initialHistoryOfVectorOfUniformState(indexOfPeriod, "Mid");
				solve.serviceLevelForLowPriorityCustomer.setServiceLevel(indexOfPeriod,
					this->serviceLevel(indexOfPeriod, shiftsEnv, solve, torWaitingTimePL, P_iq)
				);
				cout << "test0:sl\t" << indexOfPeriod << "\t" << solve.serviceLevelForLowPriorityCustomer.getServiceLevelForPeriod(indexOfPeriod) << endl;
			}
			endTimeClock[4] = clock(); 
			cout << "3\t" << double(double(endTimeClock[4] - startTimeClock[4]) / CLOCKS_PER_SEC) << endl;
			solve.serviceLevelForSeveralCustomers[2] = solve.serviceLevelForLowPriorityCustomer;


			this->clearMap();
			valueOfPai_Z_T_MAP.clear();

			for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod] + Lambda_initial[2][indexOfPeriod];
			}
			torWaitingTimePL = tor_SeveralCustomers[3];
			this->solveForApproximate = solve;
			startTimeClock[5] = clock();
			//#pragma omp parallel 
			//#pragma omp for
			//#pragma omp parallel for schedule(dynamic, 10) 
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				HistoryOfVectorOfUniformState P_iq;
				P_iq.initialHistoryOfVectorOfUniformState(indexOfPeriod, "Low");
				solve.serviceLevelForLowPriorityCustomer.setServiceLevel(indexOfPeriod,
					this->serviceLevel(indexOfPeriod, shiftsEnv, solve, torWaitingTimePL, P_iq)
					);
				cout << "test1:sl\t" << indexOfPeriod << "\t" << solve.serviceLevelForLowPriorityCustomer.getServiceLevelForPeriod(indexOfPeriod) << endl;
			}
			endTimeClock[5] = clock(); 
			cout << "4\t" << double(double(endTimeClock[5] - startTimeClock[5]) / CLOCKS_PER_SEC) << endl;
			solve.serviceLevelForSeveralCustomers[3] = solve.serviceLevelForLowPriorityCustomer;

			break;
		default:
			break;
		}
	}
	showServiceLevelByPeriod_SeveralType(solve.serviceLevelForSeveralCustomers);
	return 0;
}

double APPROXIMATE::skellam_dealwith0(int indexOfArrivalRate, int indexOfServiceRate, double arrivalRate, double serviceRate, int value) {
	//cout << "000:\t" << indexOfArrivalRate << "\t" << indexOfServiceRate << "\t" << arrivalRate << "\t" << serviceRate << "\t" << value << endl;

	if (indexOfArrivalRate == 0 && indexOfServiceRate == 0) { return double(value>=0);/*double(value == 0);*/ }
	else if (indexOfArrivalRate == 0) {
		//cout << "111:\t" << indexOfArrivalRate << "\t" << indexOfServiceRate << "\t" << arrivalRate << "\t" << serviceRate << "\t" << value << endl;
		if (value < 0)return 0;
		//cout << "value\t" << util.getPoisson(serviceRate, value) << endl;
		return util.getPoissonAccum(serviceRate, value);
	}
	else if (indexOfServiceRate == 0) {
		//cout <<"222:\t"<< indexOfArrivalRate << "\t" << indexOfServiceRate << "\t" << arrivalRate << "\t" << serviceRate << "\t" << value << endl;
		if (value > 0)return 0;
		//cout << "value\t" << util.getPoisson(serviceRate, value) << endl;
		return 1.0-util.getPoissonAccum(arrivalRate, -value);
	}
}
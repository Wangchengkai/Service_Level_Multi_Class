#pragma once
#include "CplexSolve.h"
using namespace std;

void showSolutionsByShift(SolutionOnShift solutionOnShift) {
	cout << "showSolutionsByShift" << endl;
	for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
		cout << "indexOfShift: " << indexOfShifts << " value: " << solutionOnShift.getServerOnShift(indexOfShifts) << endl;
}

void showSolutionsByShift_ToPython_Draw(SolutionOnShift solutionOnShift,Shifts shiftEnv) {
	cout << "showSolutionsByShift_ToPY" << endl;
	for (int indexOfShifts = 1; indexOfShifts < shiftEnv.totalNumOfShifts; indexOfShifts++)
		cout << solutionOnShift.getServerOnShift(indexOfShifts) << ",";
	cout << solutionOnShift.getServerOnShift(shiftEnv.totalNumOfShifts) << endl;
}

void showSolutionsByPeriod(SolutionOfServers solutionOfServers) {
	cout << "showSolutionsByPeriod" << endl;
	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
		cout << "indexOfPeriod: " << indexOfPeriods 
			<< " value: " << solutionOfServers.getSolutionOfServersForTime(indexOfPeriods) << endl;
	}
}

void showServiceLevelByPeriod(ServiceLevel serviceLevel) {
	//cout << "showServiceLevelByPeriod" << endl;
	//for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
	//	cout << "indexOfPeriod: " << indexOfPeriods
	//		<< " value: " << serviceLevel.getServiceLevelForPeriod(indexOfPeriods) << endl;
	//}
	cout << "showServiceLevelByPeriod" << endl;
	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
		cout << serviceLevel.getServiceLevelForPeriod(indexOfPeriods) << endl;
	}

}

void showServiceLevelByPeriod_SeveralType(ServiceLevel *serviceLevel) {
	//for (int indexOfType = 1; indexOfType <= 4; indexOfType++)
	//{
	//	cout << "showServiceLevelOfType\t"<< indexOfType<<"\t:" <<endl;
	//	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
	//		cout << serviceLevel[indexOfType].getServiceLevelForPeriod(indexOfPeriods) << endl;
	//	}
	//}
	cout << "showServiceLevelOfType\t:" << endl;
	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
		for (int indexOfType = 1; indexOfType <= NumberOfPriors; indexOfType++)
		{
			cout << serviceLevel[indexOfType].getServiceLevelForPeriod(indexOfPeriods) << ",";
		}
		cout << endl;
	}
}



void showShiftType(Shifts shifts) {
	cout << "showShiftType" << endl;
	for (int indexOfShifts = 1; indexOfShifts <= shifts.totalNumOfShifts; indexOfShifts++) {
		cout << "indexOfShifts: " << indexOfShifts
			<< " from: " << shifts.shifts[indexOfShifts].startTime
			<< " to: " << shifts.shifts[indexOfShifts].endTime
			<< " length " << shifts.shifts[indexOfShifts].length << endl;
	}
}

void showLowBoundOfServers(LowBoundOfServers lowBoundOfServers) {
	cout << "showLowBoundOfServers" << endl;
	for (int indexOfPeriods = 1; indexOfPeriods <= NumberOfPeriodsCondsidered; indexOfPeriods++) {
		cout << lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriods) << "\t";
		//cout << "indexOfPeriods: " << indexOfPeriods
			//<< " LBvalue " << lowBoundOfServers.getLowBoundOfServersForTime(indexOfPeriods) << endl;
	}
	cout << endl;
}

void showPeriodsUnitsLessThanThreshold(PeriodUnitsNotMeetSolutionLevel PUs) {
	cout << "showPeriodsUnitsLessThanThreshold" << endl;
	for (int indexOfUnits = 1; indexOfUnits <= PUs.totalNumOfPeriodUnits; indexOfUnits++) {
		cout << "indexOfUnits: " << indexOfUnits
			<< " PeriodsUnits StartAt: " << PUs.getStartTime(indexOfUnits) 
			<< " EndAt: " << PUs.getEndTime(indexOfUnits) << endl;
	}

}

void showServiceRate(){
	cout << "show the service rate of the low priority customers\n";
	for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
		cout << MiuForEachPeriod[indexOfPeriod] << " ";
	}
	cout << endl;
}

void showArrivalRate(string str="all") {
	if (str == "low") {
		cout << "show the arrival rate of the low priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << LamdaLP[indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "mid") {
		cout << "show the arrival rate of the mid priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << LamdaMP[indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "high") {
		cout << "show the arrival rate of the high priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << Lamda[indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "all") {
		showArrivalRate("low");
		showArrivalRate("mid");
		showArrivalRate("high");
	}
	else { cout << "showArrivalRate::input error" << endl;; }
}


void showArrivalRate_initial(string str = "all") {
	if (str == "low") {
		cout << "show the arrival rate of the low priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << Lambda_initial[3][indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "mid") {
		cout << "show the arrival rate of the mid priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << Lambda_initial[2][indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "high") {
		cout << "show the arrival rate of the high priority customers\n";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << Lambda_initial[1][indexOfPeriod] << " ";
		}
		cout << endl;
	}
	else if (str == "all") {
		showArrivalRate_initial("low");
		showArrivalRate_initial("mid");
		showArrivalRate_initial("high");
	}
	else { cout << "showArrivalRate::input error" << endl;; }
}
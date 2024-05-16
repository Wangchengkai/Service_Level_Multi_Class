#pragma once
#include "data.h";

using namespace std;



void testIndexOfShiftOutOfRange(int index);



class QueueLength {
public:
	double initialLength;
};

class PeriodOutOfRange {
public:
	void testPeriodIndex(int numOfPeriod) const {
		if (numOfPeriod > NumberOfPeriodsCondsidered) { cout << "error:hour outofrange\t" <<numOfPeriod<< endl; system("pause"); }
	}
	bool isPeriodIndexOutOfRange(int numOfPeriod) {
		if (numOfPeriod > NumberOfPeriodsCondsidered) { cout << "pause:hour outofrange1\t"<<numOfPeriod << endl; system("pause"); return true; }
		else return false;
	}
};
#include<string>
class LowBoundOfServers : public PeriodOutOfRange {
	int lowBoundOfServers[NumberOfPeriodsCondsidered];
public:
	int getLowBoundOfServersForTime(int numOfPeriod);
	void setLowBoundOfServersForTime(int numOfPeriod, int valueOfServers);
	void lowBoundMinusMinus(int numOfPeriod);
	void lowBoundPlusPlus(int numOfPeriod);
	void setLowBoundOfServersAsOne();
	string convertLowBoundToString() {
		string out = "";
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			out.append(std::to_string(this->getLowBoundOfServersForTime(indexOfPeriod)));
		return out;
	}

};



class SolutionOfServers : public PeriodOutOfRange {
	int solutionOfServers[NumberOfPeriodsCondsidered];
public:
	int getSolutionOfServersForTime(int numOfPeriod) const;
	void setSolutionOfServersForTime(int numOfPeriod, int valueOfServers);
};



class InformationOfShift {
public:
	int startTime;
	int endTime;
	int length;
};

class Shift:public InformationOfShift, PeriodOutOfRange{
	int indexOfShift;
public:
	void setShift(int index,int startTime,int endTime) {
		indexOfShift = index;
		this->startTime = startTime;
		this->endTime = endTime;
		this->length = endTime - startTime;	}
	bool periodInShift(int indexOfPeriod) {
		testIndexOfShiftOutOfRange(indexOfShift);
		PeriodOutOfRange P1;
		P1.testPeriodIndex(indexOfPeriod);
		if (indexOfPeriod <= this->endTime
			&& indexOfPeriod > this->startTime)return true;
		else return false;
	}
};

class Shifts: public PeriodOutOfRange {
public:
	int totalNumOfShifts;
	Shift shifts[MaxNumOfShifts + 1];	void setNightShift(int index);
	void setAllShifts();
};

class SolutionOnShift {
	int serverOnShift[MaxNumOfShifts + 1] = { 0 };
public:
	int getServerOnShift(int index) const {
		testIndexOfShiftOutOfRange(index);
		return serverOnShift[index - 1];
	}
	void setServerOnShift(int index, int value) {
		testIndexOfShiftOutOfRange(index);
		serverOnShift[index - 1] = value;
	}
	int getSumOfServersOnPeriod(int indexOfPeriod, Shifts shiftEnv) {
		int count = 0;
		for (int indexOfShifts = 1; indexOfShifts <= MaxNumOfShifts; indexOfShifts++)
			count += getServerOnShift(indexOfShifts) * shiftEnv.shifts[indexOfShifts].periodInShift(indexOfPeriod);
		return count;
	}
	string convertSolutionOnShiftToString() {
		string out = "";
		for (int indexOfShift = 1; indexOfShift <= MaxNumOfShifts; indexOfShift++)
			out.append(std::to_string(this->getServerOnShift(indexOfShift)));
		return out;
	}
};
class ComputationOfServicelLevel;

class ServiceLevel : public PeriodOutOfRange {
	double serviceLevel[NumberOfPeriodsCondsidered];
public:
	void computeServiceLevel(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShifts);
	ComputationOfServicelLevel initialComputeServiceLevel_TwoTypes(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShifts);
	ComputationOfServicelLevel initialComputeServiceLevel_ThreeTypes(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShifts);
	void computeServiceLevel_forHighPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShifts, ComputationOfServicelLevel computationOfServicelLevelInNewWay);
	void computeServiceLevel_forLowPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShifts, ComputationOfServicelLevel computationOfServicelLevelInNewWay);
	void computeServiceLevel_forLowPriority(Shifts shiftEnv, SolutionOfServers solutionOfServers, SolutionOnShift solutionOnShift, ComputationOfServicelLevel computationOfServicelLevelInNewWay, string type);

	//double computeServiceLevelForPeriod(Shifts shiftEnv, SolutionOfServers solutionOfServers, int indexOfPeriod);
	double getServiceLevelForPeriod(int numOfPeriod);
	void setServiceLevel(int numOfPeriod, double valueOfSL);
	bool serviceLevelLessThanThreshold(int numOfPeriod,double threshold);
	bool serviceLevelMoreThanHighThreshold(int numOfPeriod, double thresholdOfSLTooHigh);
	bool serviceLevelBothMoreThanThreshold(double threshold) {
		int isBothMoreThanThreshold = true;
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
			if (this->serviceLevelLessThanThreshold(indexOfPeriod, threshold))isBothMoreThanThreshold = false;
		}
		return isBothMoreThanThreshold;
	}
	double getTotalGap() {
		double sum = 0;
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			sum += max(0.0, (this->getServiceLevelForPeriod(indexOfPeriod) - ThresholdOfServiceLevelList[indexOfPeriod]));
		return sum;
	}
	bool isServiveCapcityEnough(SolutionOfServers solutionOfServers);
};


class PeriodUnitsNotMeetSolutionLevel {
	int startTimeOfPeriodUnit[NumberOfPeriodsCondsidered];
	int endTimeOfPeriodUnit[NumberOfPeriodsCondsidered];public:
	int totalNumOfPeriodUnits = 0;
	void setPeriodUnits(ServiceLevel *SL);
	void DoIfServiceLevelLessThanThreshold(int* countNumOfPeriodUnits, int* numofperiod, ServiceLevel *SL);
	int getStartTime(int numberOfPeriodUnit) {
		PeriodOutOfRange P1;
		P1.testPeriodIndex(numberOfPeriodUnit);
		return startTimeOfPeriodUnit[numberOfPeriodUnit - 1];
	}
	int getEndTime(int numberOfPeriodUnit) {
		PeriodOutOfRange P1;
		P1.testPeriodIndex(numberOfPeriodUnit);
		return endTimeOfPeriodUnit[numberOfPeriodUnit - 1];
	}
	void initialPeriodUnits() {
		for (int indexOfUnit = 0; indexOfUnit < NumberOfPeriodsCondsidered; indexOfUnit++) {
			startTimeOfPeriodUnit[indexOfUnit] = 0;
			endTimeOfPeriodUnit[indexOfUnit] = 0;
		}
	}
};



class Solve1{
public:	
	ServiceLevel serviceLevel;
	ServiceLevel serviceLevelForLowPriorityCustomer;
	
	ServiceLevel serviceLevelForSeveralCustomers[5];
	
	
	PeriodUnitsNotMeetSolutionLevel periodUnits;
	SolutionOfServers currentSolutionOfServers;
	SolutionOnShift currentSolutionOnShift;
	//QueueLength queueLength;
	LowBoundOfServers lowBoundOfServers;



//functions:
	void forOneUnit(int indexOfUnit);
	void updataLowBoundForOneUnitToPlusPlus(int unitLeft,int unitRight);
	void updataServersNumForOneUnitByLowBound(int unitLeft, int unitRight);
	SolutionOfServers tryToMinusServers(int unitLeft, int unitRight);
	void minusServerInOnePeriodTo_LBminusOne(SolutionOfServers* solutionOfOneServer, int period);
	bool bothMoreThanSLThreshold(ServiceLevel *SL, int left, int right);
	bool tryToMinusLowBound(int unitLeft, int unitRight);

	void computeServiceLevel_SeveralType();
};

int computeOffShiftInPeriod(Shifts shiftEnv, SolutionOnShift solutionOnShift, int endTimeOfOffShift);
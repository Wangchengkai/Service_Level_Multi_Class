#pragma once
#ifndef _TWOTYPESL_

#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include "data.h"

using namespace std;

#define MaxNumOfK 500
#define MaxNumOfState MaxValueOfStatesInUniform
#define MaxNumOfPeriod NumberOfPeriodsCondsidered
#define MaxNumOfServer TotalNumberOfServers
//#define EPSILON 1.0e-5
#define tao torWaitingTime
//#define LengthOfOnePeriod 1.0


//int stateOfWaitingAnd1;

//API
//{ 0,//dayOne
//	9,2,3,5,4,1,3,16,25,29,14,20,15,11,14,9,11,17,14,11,15,24,29,14
//};
//#define Miu 5.9113
int numOfServers[MaxNumOfPeriod + 1] =
{ 0,
	1,2,3,4,5,1,2,3//,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4
};
int onShift[MaxNumOfPeriod + 1] = { 0,
	5,6,1,1,1,0,1,1,//1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,0
};






double Q_matrix[MaxNumOfState + 1][MaxNumOfState + 1];

class Q_matrixManage {
	void setQ_martix(int fromTotalState, int toTotalState, double value) {
		Q_matrix[fromTotalState][toTotalState] = value;
	}
public:
	double static getQ_martix(int fromTotalState, int toTotalState) { return Q_matrix[fromTotalState][toTotalState]; }
	void static initialQ_martix(int indexOfPeriod) {
		for (int fromTotalState = 0; fromTotalState <= MaxNumOfState; fromTotalState++)
			for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++)
				Q_matrix[fromTotalState][toTotalState] = 0;

		Q_matrix[0][0] = 0;
		Q_matrix[0][1] = 0;
		Q_matrix[MaxNumOfState][MaxNumOfState - 1] = double(numOfServers[indexOfPeriod]) * double(MiuForEachPeriod[indexOfPeriod]);
		Q_matrix[MaxNumOfState][MaxNumOfState] = 0.0 - Q_matrix[MaxNumOfState][MaxNumOfState - 1];
		for (int fromTotalState = 1; fromTotalState < MaxNumOfState; fromTotalState++)
			for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++) {
				if (fromTotalState - 1 == toTotalState)Q_matrix[fromTotalState][toTotalState] = double(numOfServers[indexOfPeriod]) * double(MiuForEachPeriod[indexOfPeriod]);
				if (fromTotalState + 1 == toTotalState)Q_matrix[fromTotalState][toTotalState] = Lamda[indexOfPeriod];
			}
		for (int fromTotalState = 1; fromTotalState <= MaxNumOfState; fromTotalState++)
			Q_matrix[fromTotalState][fromTotalState] = 0.0
			- Q_matrix[fromTotalState][fromTotalState - 1]
			- Q_matrix[fromTotalState][fromTotalState + 1];
	}
	void static showQ_martix() {
		cout << "showQ_matrix" << endl;
		for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++)
		{
			for (int fromTotalState = 0; fromTotalState <= MaxNumOfState; fromTotalState++)
				cout << setw(8) << setprecision(4) << getQ_martix(fromTotalState, toTotalState);
			cout << endl;
		}
	}
};

class VectorOfStateProb {
	double stateProb[MaxNumOfState + 1] = { 0 };
public:
	double getProbOfState(int indexOfState) {
		return stateProb[indexOfState];
	}
	void setProbOfState(int indexOfState, double value) {
		stateProb[indexOfState] = value;
	}
	void initialProbOfState() {//(int& prob) {
		for (int indexOfState = 0; indexOfState < 10; indexOfState++) {
			setProbOfState(indexOfState, 0.1);
		}
	}
	VectorOfStateProb multiQ_martix(int multiTime_k, double tInThePeriod) {
		VectorOfStateProb outVectorOfStateProb;

		int tempToState = 0;
		outVectorOfStateProb.setProbOfState(tempToState,
			Q_matrixManage::getQ_martix(0, tempToState) * getProbOfState(0) +
			Q_matrixManage::getQ_martix(1, tempToState) * getProbOfState(1));

		for (tempToState = 1; tempToState < MaxNumOfState; tempToState++)
		{
			outVectorOfStateProb.setProbOfState(tempToState,
				Q_matrixManage::getQ_martix(tempToState - 1, tempToState) * getProbOfState(tempToState - 1) +
				Q_matrixManage::getQ_martix(tempToState + 0, tempToState) * getProbOfState(tempToState + 0) +
				Q_matrixManage::getQ_martix(tempToState + 1, tempToState) * getProbOfState(tempToState + 1)
			);
		}

		tempToState = MaxNumOfState;
		outVectorOfStateProb.setProbOfState(tempToState,
			Q_matrixManage::getQ_martix(tempToState - 1, tempToState) * getProbOfState(tempToState - 1) +
			Q_matrixManage::getQ_martix(tempToState + 0, tempToState) * getProbOfState(tempToState + 0));


		for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++)
			outVectorOfStateProb.setProbOfState(toTotalState,
				outVectorOfStateProb.getProbOfState(toTotalState) / double(multiTime_k) * double(tInThePeriod));

		return outVectorOfStateProb;
	}

	void add(VectorOfStateProb addedVector) {
		for (int numOfState = 0; numOfState <= MaxNumOfState; numOfState++)
			this->setProbOfState(numOfState,
				this->getProbOfState(numOfState) + addedVector.getProbOfState(numOfState));
	}

	VectorOfStateProb cumMultiQ_martix(double tInThisPeriod, VectorOfStateProb initialVector) {
		//cout << "test2\n";

		VectorOfStateProb tempMultiVector = initialVector;
		VectorOfStateProb cumMultiVector = initialVector;
		for (int multiTime_k = 1; multiTime_k <= MaxNumOfK; multiTime_k++) {
			tempMultiVector = tempMultiVector.multiQ_martix(multiTime_k, tInThisPeriod);
			cumMultiVector.add(tempMultiVector);

			double sum = 0;
			for (int numOfState = 0; numOfState <= MaxNumOfState; numOfState++) {
				sum += abs(tempMultiVector.getProbOfState(numOfState));
				//cout << "sum\t" << numOfState << "\t" << sum << endl;
			}
			if (sum <= EPSILON) {
								return cumMultiVector;
			}
		}
		system("pause");
		cout << "warning 2\t xuyaozengjia K" << endl;

		return cumMultiVector;
	}

	VectorOfStateProb changeStateBetweenPeriods(int indexOfPeriod, VectorOfStateProb beforeChangeVector) {
				while (indexOfPeriod > MaxNumOfPeriod)indexOfPeriod -= MaxNumOfPeriod;

		VectorOfStateProb afterChangeVector;
		for (int indexOfState = 0; indexOfState <= MaxNumOfState; indexOfState++) {
			afterChangeVector.setProbOfState(
				max(0, (indexOfState - onShift[indexOfPeriod])),
				afterChangeVector.getProbOfState(max(0, (indexOfState - onShift[indexOfPeriod])))
				+ beforeChangeVector.getProbOfState(indexOfState)
			);
		}
		return afterChangeVector;
	}

	VectorOfStateProb stepByStep(double nowTime, VectorOfStateProb initialState) {
		double remainingTime = tao;
		while (remainingTime > 0) {
			int indexOfPeriod = int(nowTime) + 1; if (nowTime == int(nowTime))indexOfPeriod--;
			double considerTimeInThisPeriod = min(remainingTime, (LengthOfOnePeriod - fmod(nowTime, LengthOfOnePeriod)));
						
			Q_matrixManage::initialQ_martix(indexOfPeriod);
			//Q_matrixManage::showQ_martix();


			initialState = initialState.cumMultiQ_martix(considerTimeInThisPeriod, initialState);



			remainingTime -= considerTimeInThisPeriod;
			if (remainingTime > 0) {
				initialState = initialState.changeStateBetweenPeriods(indexOfPeriod + 1, initialState);			}
			nowTime += considerTimeInThisPeriod;


		}

		return initialState;
	}

	double getSLByVectorForArringAtTimeT() {
		return this->getProbOfState(0);
	}

	void showVectorUtill(int numOfShowStates = MaxNumOfState) {
		cout << "showPartOfStateVector" << endl;
		for (int i = 0; i <= numOfShowStates; i++)
			cout << i << ":\t" << this->getProbOfState(i) << endl;

	}

};





#define PrecisionOfIntegralOfTInOnePeriod 200

int main() {
	int indexOfPeriod = 1;
	long double sumOfIntegral = 0;
	for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++) {
		double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;

		VectorOfStateProb vectorOfStateProb;


		vectorOfStateProb.initialProbOfState();

		vectorOfStateProb = vectorOfStateProb.stepByStep(arrivingTime, vectorOfStateProb);
		
		
		cout <<"arrivingTime: "<<arrivingTime << " SL:\t" << vectorOfStateProb.getSLByVectorForArringAtTimeT() << endl;
		
		cout << vectorOfStateProb.getSLByVectorForArringAtTimeT() << endl;

		sumOfIntegral += vectorOfStateProb.getSLByVectorForArringAtTimeT();

		//vectorOfStateProb.showVectorUtill(20);
	}
	sumOfIntegral = sumOfIntegral / double(PrecisionOfIntegralOfTInOnePeriod);

	cout << " SL:\t" << sumOfIntegral << endl;


	system("pause");
	system("pause");

}



#endif // !_TWOTYPESL_

#pragma once
#ifndef _computeLowPriority_
#define _computeLowPriority_

#include <vector>
#include <iostream>
#include <iomanip>
#include <Eigen\Dense>
//#include <Eigen/MatrixFunctions>
#include <unsupported/Eigen/MatrixFunctions>
#include <time.h>
#include "data.h"


clock_t startTimeClock_TwoType[10] = { 0 }, endTimeClock_TwoType[10] = { 0 };

using namespace Eigen;

using namespace std;




//int stateOfWaitingAnd1;

//API
//{ 0,//dayOne
//	9,2,3,5,4, 1,3,16,25,29, 14,20,15,11,14,9,11,17,14,11,15,24,29,14
//	//4,1,1,2,2, 1,1,8,12,14, 7,10,7,5,7,4,5,8,7,5,7,12,14,7
//};
//#define Miu 5.9113
int numOfServers[MaxNumOfPeriod + 1] =
{ 0,
	1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4
};
int onShift[MaxNumOfPeriod + 1] = { 0,
	5,6,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,0
};




double Q_matrix[MaxNumOfState + 1][MaxNumOfState + 1];
//Matrix<double, MaxNumOfState + 1, MaxNumOfState + 1> ExpQ_eigen;
//Matrix<double, MaxNumOfState + 1, MaxNumOfState + 1> Q_eigen;
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ExpQ_eigen;
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q_eigen;

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
		Q_matrix[MaxNumOfState][MaxNumOfState - 1] = double(numOfServers[indexOfPeriod]) * MiuForEachPeriod[indexOfPeriod];
		Q_matrix[MaxNumOfState][MaxNumOfState] = 0.0 - Q_matrix[MaxNumOfState][MaxNumOfState - 1];
		for (int fromTotalState = 1; fromTotalState < MaxNumOfState; fromTotalState++)
			for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++) {
				if (fromTotalState - 1 == toTotalState)Q_matrix[fromTotalState][toTotalState] = double(numOfServers[indexOfPeriod]) * MiuForEachPeriod[indexOfPeriod];
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

	void static convertToEigen() {
		Q_eigen = Eigen::MatrixXd::Zero(MaxNumOfState + 1, MaxNumOfState + 1);
		for (int toTotalState = 0; toTotalState <= MaxNumOfState; toTotalState++)
		{
			for (int fromTotalState = 0; fromTotalState <= MaxNumOfState; fromTotalState++)
				Q_eigen(toTotalState, fromTotalState) = getQ_martix(fromTotalState, toTotalState);
		}
	}
};

int isUsingEigenTemp = IsUsingEigen;


class VectorOfStateProb {
	long double stateProb[MaxNumOfState + 1] = { 0 };
public:
	double getProbOfState(int indexOfState) {
		return stateProb[indexOfState];
	}
	void setProbOfState(int indexOfState, double value) {
		stateProb[indexOfState] = value;
	}

	double confirmSumOfProb() {
		double sum = 0;
		for (int tempState = 0; tempState <= MaxNumOfState; tempState++)
			sum += this->getProbOfState(tempState);
		return sum;

	}

	void initialProbOfState() {//(int& prob) {
		for (int indexOfState = 0; indexOfState <= MaxNumOfState; indexOfState++) {
			setProbOfState(indexOfState, 0.0);
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
			}

			if (sum <= EPSILON) {
				return cumMultiVector;
			}
			else if (sum > 1e15) {
								isUsingEigenTemp = true;
				return cumMultiVector;
			}
		}
		system("pause");
		cout << "warning 2\t 需要增加 K\t" << tInThisPeriod << "\t" << endl;
		initialVector.showVectorUtill();
		//for (int multiTime_k = 1; multiTime_k <= MaxNumOfK; multiTime_k++) {
		//	double sum = 0;
		//	for (int numOfState = 0; numOfState <= MaxNumOfState; numOfState++) {
		//		sum += abs(tempMultiVector.getProbOfState(numOfState));
		//	}
		//	cout << sum << endl;
		//}
		cout << Q_matrix[100][100 - 1] << endl;

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

	VectorOfStateProb computeByEigen(double nowTime, VectorOfStateProb initialState) {
		double remainingTime = torWaitingTimePL;
		while (remainingTime > 0) {

			int indexOfPeriod = int(nowTime) + 1; //if (nowTime == int(nowTime))indexOfPeriod--;
			double considerTimeInThisPeriod = min(remainingTime, (LengthOfOnePeriod - fmod(nowTime, LengthOfOnePeriod)));
			while (indexOfPeriod > NumberOfPeriodsCondsidered)indexOfPeriod -= NumberOfPeriodsCondsidered;
			/*	cout << "剩余时长：" << remainingTime << "\t当前时间: " << nowTime << endl;
				cout << "因此本阶段: " << indexOfPeriod << " 考虑的时长为: " << considerTimeInThisPeriod << endl;*/

			Q_matrixManage::initialQ_martix(indexOfPeriod);
			Q_matrixManage::convertToEigen();
			Q_eigen = Q_eigen * considerTimeInThisPeriod;
			//cout << Q_eigen << endl;

			ExpQ_eigen = Q_eigen.exp();

			Matrix<double, MaxNumOfState + 1, 1> vectorOfState;
			for (int indexOfState = 0; indexOfState <= MaxNumOfState; indexOfState++) {
				vectorOfState(indexOfState) = initialState.getProbOfState(indexOfState);
			}
			vectorOfState = ExpQ_eigen * vectorOfState;
			for (int indexOfState = 0; indexOfState <= MaxNumOfState; indexOfState++) {
				initialState.setProbOfState(indexOfState, vectorOfState(indexOfState));
			}


			remainingTime -= considerTimeInThisPeriod;
			if (remainingTime > 0) {
				initialState = initialState.changeStateBetweenPeriods(indexOfPeriod + 1, initialState);			}
			nowTime += considerTimeInThisPeriod;
		}
		//initialState.showVectorUtill(1);
		//system("pause");
		return initialState;
	}

	VectorOfStateProb stepByStep(double nowTime, VectorOfStateProb initialState) {
		VectorOfStateProb oldInitialState = initialState;
		double oldNowTime = nowTime;

		double remainingTime = torWaitingTimePL;
		while (remainingTime > 0) {
			int indexOfPeriod = int(nowTime) + 1;// if (nowTime == int(nowTime))indexOfPeriod--;
			double considerTimeInThisPeriod = min(remainingTime, (LengthOfOnePeriod - fmod(nowTime, LengthOfOnePeriod)));
			while (indexOfPeriod > NumberOfPeriodsCondsidered)indexOfPeriod -= NumberOfPeriodsCondsidered;
			//cout << "test\t" << indexOfPeriod << endl;
			/*cout << "剩余时长：" << remainingTime << "\t当前时间: " << nowTime << endl;
			cout << "因此本阶段: " << indexOfPeriod << " 考虑的时长为: " << considerTimeInThisPeriod << endl;*/

			Q_matrixManage::initialQ_martix(indexOfPeriod);
			//Q_matrixManage::showQ_martix();

			
			initialState = initialState.cumMultiQ_martix(considerTimeInThisPeriod, initialState);
			if (isUsingEigenTemp == true) {
				isUsingEigenTemp = IsUsingEigen;
								return computeByEigen(oldNowTime, oldInitialState);
			}


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

	
	void static API_numOfServers(int indexOfPeriod, int valueOfServers) {
		numOfServers[indexOfPeriod] = valueOfServers;
	}
	void static API_numOfOnShifts(int indexOfPeriod, int valueOfServers) {
		onShift[indexOfPeriod] = valueOfServers;
	}

};





#endif // !_computeLowPriority_



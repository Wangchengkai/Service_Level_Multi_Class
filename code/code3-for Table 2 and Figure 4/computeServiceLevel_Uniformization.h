#pragma once
#ifndef _computeServiceLevel_Uniformization_
#define _computeServiceLevel_Uniformization_

#include "OperFunctions.h"
#include "functionsForXJ.h"
#include <vector>
#include "computeLowPriority.h"
#include <omp.h>
const int numProcs = 16;

long double F_N[PrecisionOfPoissonParameter][MaxNumberOfTranstionTimesInUniform + 1] = { 0 };
long double F_N_Accumulate[PrecisionOfPoissonParameter][MaxNumberOfTranstionTimesInUniform + 1] = { 0 };
class Utils {
	static const int N_combineNum = 400;
	long double comb[N_combineNum][N_combineNum];	double maxValueOfGama = 10e-5;
public:
	void initialCombineNum() {
		for (int i = 0; i < N_combineNum; i++) {
			comb[i][0] = comb[i][i] = 1;
			for (int j = 1; j < i; j++) {
				comb[i][j] = comb[i - 1][j] + comb[i - 1][j - 1];
			}
		}
	}
	long double getCombineNum(int select, int from)const {
		if (select > from) {
			cout << "getCombineNum::error1" << select << "\t" << from << endl; system("pause");
		}
		else if (select >199|| from>199) {
			cout << "getCombineNum::error2 orr" << select << "\t" << from << endl; system("pause");
		}
		return comb[from][select];
	}


	long double getHyperGeometric(int busySelectServers, int totalSelectServers, int totalServers, int busyServers)const {
		if (busySelectServers > busyServers || totalSelectServers > totalServers || (totalSelectServers - busySelectServers) > (totalServers - busyServers)) {
			cout << busySelectServers << "\t" << totalSelectServers << "\t" << totalServers << "\t" << busyServers << endl;
			system("pause");
		}

		
		
		
		if (totalSelectServers > totalServers) {
			cout << "getHyperGeometric::error1" << totalSelectServers << "\t" << totalServers << endl; system("pause");
		}
		if (busySelectServers > totalSelectServers) { cout << "error1:" << busySelectServers << "\t" << totalSelectServers << endl; }
		if ((totalSelectServers - busySelectServers) > (totalServers - busyServers)) {
			cout << "error2:" << (totalSelectServers - busySelectServers) << "\t" << (totalServers - busyServers) << endl;
		}
		if (totalSelectServers > totalServers) { cout << "error3:" << totalSelectServers << "\t" << totalServers << endl; }

		if (getCombineNum(busySelectServers, busyServers)
			* getCombineNum((totalSelectServers - busySelectServers), (totalServers - busyServers))
			/ getCombineNum(totalSelectServers, totalServers) > 1) {
			cout << "error4:\t" <<
				getCombineNum(busySelectServers, busyServers) << "\t" <<
				getCombineNum((totalSelectServers - busySelectServers), (totalServers - busyServers)) << "\t" <<
				getCombineNum(totalSelectServers, totalServers)
				<< endl;
		}

		return getCombineNum(busySelectServers, busyServers)
			* getCombineNum((totalSelectServers - busySelectServers), (totalServers - busyServers))
			/ getCombineNum(totalSelectServers, totalServers);
	}

	void initialPoisson() {
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			maxValueOfGama = max(maxValueOfGama, double(TotalNumberOfServers) * MiuForEachPeriod[indexOfPeriod] + Lamda[indexOfPeriod] + LamdaLP[indexOfPeriod]);
				double maxValueOfMiu = 0;
		for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)
			maxValueOfMiu = max(maxValueOfMiu, MiuForEachPeriod[indexOfPeriod]);

		double maxValueOfA = torWaitingTime * double(TotalNumberOfServers) * maxValueOfMiu;
		maxValueOfGama = max(maxValueOfGama, maxValueOfA);

		maxValueOfGama += 1.0;		cout << "maxValueOfGama\t" << maxValueOfGama << endl;

		for (int indexOfPara = 0; indexOfPara < PrecisionOfPoissonParameter; indexOfPara++) {
			double gamaNow = double(indexOfPara) / double(PrecisionOfPoissonParameter - 1) * maxValueOfGama;
			double ExpFuGama = exp(-1.0 * (gamaNow));
			F_N[indexOfPara][0] = ExpFuGama;
			for (int Nnum = 1; Nnum <= MaxNumberOfTranstionTimesInUniform; Nnum++)
				F_N[indexOfPara][Nnum] = F_N[indexOfPara][Nnum - 1] * (gamaNow) / double(Nnum);
		}
	}

	void initialPoissonAccum() {
		for (int indexOfPara = 0; indexOfPara < PrecisionOfPoissonParameter; indexOfPara++) {
			long double accum = 0;
			for (int Nnum = 0; Nnum <= MaxNumberOfTranstionTimesInUniform; Nnum++) {
				accum += F_N[indexOfPara][Nnum];
				F_N_Accumulate[indexOfPara][Nnum] = accum;
			}
		}
	}

	long double getPoisson(double gama, int index) const {
		int indexOfGama = int(gama / (maxValueOfGama) * double(PrecisionOfPoissonParameter));
		if (indexOfGama < 0 || indexOfGama >= PrecisionOfPoissonParameter) {
			cout << "Poisson error:indexOfGama OutOfRange\t" << indexOfGama << endl; system("pause");
		}
		else if (index<0 || index>MaxNumberOfTranstionTimesInUniform) {
			cout << "Poisson error:indexOfTran OutOfRange\t" << index << endl; system("pause");
		}
		else
			return F_N[indexOfGama][index];
	}
	long double getPoissonAccum(double para, int index)const {
		int indexOfpara = int(para / (maxValueOfGama) * double(PrecisionOfPoissonParameter));
		if (indexOfpara < 0 || indexOfpara >= PrecisionOfPoissonParameter) {
			cout << "Poisson error1:indexOfpama OutOfRange\t" << indexOfpara << endl; system("pause");
		}
		else if (index<0 || index>MaxNumberOfTranstionTimesInUniform) {
			cout << "Poisson error1:indexOfserver OutOfRange\t" << index << endl; system("pause");
		}
		else
			return F_N_Accumulate[indexOfpara][index];
	}

	void initial() {
		initialCombineNum();
		initialPoisson();
		initialPoissonAccum();
	}
};

Utils util;

vector<vector<vector<vector<double>>>> FXY_uniform_Piq;
//double FXY_uniform_Piq[NumberOfPeriodsCondsidered][MaxNumberOfTranstionTimesInUniform + 1] [MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1];

vector < vector<vector<double>>> FXY_uniform_Piq_Mid;
vector < vector<vector<double>>> FXY_uniform_Piq_Low;


//tiqq
class GetFXY {
public:
	double static get(int indexOfPeriod, int numberOfTran, int stateOfLow, int stateOfHigh) {
		if (indexOfPeriod > NumberOfPeriodsCondsidered || indexOfPeriod <= 0) { cout << "error_getFXY" << endl; system("pause"); }
		return FXY_uniform_Piq[indexOfPeriod - 1][numberOfTran][stateOfLow][stateOfHigh];
	}
	double static get3_sys_sys_serve(int indexOfPeriod, int numberOfTran, int stateOfSum) {
		if (indexOfPeriod > NumberOfPeriodsCondsidered || indexOfPeriod <= 0) { cout << "error_getFXY" << endl; system("pause"); }
		return FXY_uniform_Piq_Mid[indexOfPeriod - 1][numberOfTran][stateOfSum];//todo
	}
	double static get3_sys_sys_sys(int indexOfPeriod, int numberOfTran, int stateOfSum) {
		if (indexOfPeriod > NumberOfPeriodsCondsidered || indexOfPeriod <= 0) { cout << "error_getFXY" << endl; system("pause"); }
		return FXY_uniform_Piq_Low[indexOfPeriod - 1][numberOfTran][stateOfSum];//todo
	}


	void showPiq_sym_sym_serve(int indexOfPeriod, int indexOfTran = 10) {
		cout << "show Piq s_s_serve in period: " << indexOfPeriod << "\tafter" << indexOfTran << "\t transform" << endl;
		for (int state = 0; state = QMAXSUMS; state++) {
			cout << "state: " << state
				<< " prob: " << GetFXY::get3_sys_sys_serve(indexOfPeriod, indexOfTran, state) << endl;
		}
	}

	void showPiq_sym_sym_sym(int indexOfPeriod, int indexOfTran = 10) {
		cout << "show Piq s_s_serve in period: " << indexOfPeriod << "\tafter" << indexOfTran << "\t transform" << endl;
		for (int state = 0; state = QMAXSUM; state++) {
			cout << "state: " << state
				<< " prob: " << GetFXY::get3_sys_sys_sys(indexOfPeriod, indexOfTran, state) << endl;
		}
	}

};

 
class HistoryOfVectorOfUniformState {//p_iq
	double P[MaxNumberOfTranstionTimesInUniform + 1][MaxValueOfStatesInUniform + 1] = { 0 };
	double P_total[MaxNumberOfTranstionTimesInUniform + 1][MaxValueOfStatesInUniform + 1] = { 0 };

public:
	void setProbOfHighProrityQueueLengthAfterTran(int indexOfTran, int indexOfState, double value) {
		if (indexOfState < 0 || indexOfState>MaxValueOfStatesInUniform) {
			cout << "error 1234\t" << indexOfState << endl; system("pause");
			throw "HistoryOfVectorOfUniformState::error:indexOfState OutOfRange";
		}
		else if (indexOfTran<0 || indexOfTran>MaxNumberOfTranstionTimesInUniform) {
			cout << "error 1235\t" << indexOfTran << endl; system("pause");
			throw "HistoryOfVectorOfUniformState::error:indexOfTran OutOfRange";
		}
		else {
			P[indexOfTran][indexOfState] = value; 
		}
	}
	void setProbOfTotalQueueLengthAfterTran(int indexOfTran, int indexOfState, double value) {
		if (indexOfState < 0 || indexOfState>MaxValueOfStatesInUniform)throw "HistoryOfVectorOfUniformState::error:indexOfState OutOfRange";
		else if (indexOfTran<0 || indexOfTran>MaxNumberOfTranstionTimesInUniform)throw "HistoryOfVectorOfUniformState::error:indexOfTran OutOfRange";
		else P_total[indexOfTran][indexOfState] = value;
	}

	double getProbOfHighProrityQueueLengthAfterTran(int indexOfTran, int indexOfState) const { return P[indexOfTran][indexOfState]; }
	double getProbOfTotalQueueLengthAfterTran(int indexOfTran, int indexOfState) const { return P_total[indexOfTran][indexOfState]; }
	
	double computeProbForTotalQueueLength(int indexOfPeriod, int indexOfTran, int valueOfState) {
		double totalProb = 0; //GetFXY getFxy;
		for (int valueOfHP = 0; valueOfHP <= valueOfState; valueOfHP++)
			totalProb += GetFXY::get(indexOfPeriod, indexOfTran, valueOfState - valueOfHP, valueOfHP);
		return totalProb;
	}
	double computeProbForHighProrityQueueLength(int indexOfPeriod, int indexOfTran, int valueOfHP) {
		double sumProb = 0; //GetFXY getFxy;
		for (int valueOfLP = 0; valueOfLP <= MaxValueOfStatesInUniform; valueOfLP++)
			sumProb += GetFXY::get(indexOfPeriod, indexOfTran, valueOfLP , valueOfHP);
		return sumProb;
	}
	double computeProbForTotalQueueLength_Mid(int indexOfPeriod, int indexOfTran, int valueOfState) {
		if (valueOfState == MaxValueOfStatesInUniform) {
			double sum = 0;
			for (int state = MaxValueOfStatesInUniform; state <= QMAXSUMS; state++) {
				sum += GetFXY::get3_sys_sys_serve(indexOfPeriod, indexOfTran, state);
			}
			return sum;
		}
		if (valueOfState > QMAXSUMS) { return 0; }
		return GetFXY::get3_sys_sys_serve(indexOfPeriod, indexOfTran, valueOfState);
	}
	double computeProbForTotalQueueLength_Low(int indexOfPeriod, int indexOfTran, int valueOfState) {
		if (valueOfState == MaxValueOfStatesInUniform) {
			double sum = 0;
			for (int state = MaxValueOfStatesInUniform; state <= QMAXSUM; state++)
				sum += GetFXY::get3_sys_sys_sys(indexOfPeriod, indexOfTran, state);
			return sum;
		}
		if (valueOfState > QMAXSUM) { return 0; }
		return GetFXY::get3_sys_sys_sys(indexOfPeriod, indexOfTran, valueOfState);

	}
	void initialHistoryOfVectorOfUniformState(int indexOfPeriod) {
				//GetFXY getFxy;

		for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInUniform; indexOfTran++) {
			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {
				setProbOfHighProrityQueueLengthAfterTran(indexOfTran, indexOfState,
					computeProbForHighProrityQueueLength(indexOfPeriod,indexOfTran,indexOfState));
				setProbOfTotalQueueLengthAfterTran(indexOfTran, indexOfState,
					computeProbForTotalQueueLength(indexOfPeriod,indexOfTran,indexOfState));
			}
		}

	}

	void initialHistoryOfVectorOfUniformState(int indexOfPeriod, string type) {
				//GetFXY getFxy;

		if (type == "Mid")
		{
			for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInUniform; indexOfTran++) {
				for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {
					setProbOfTotalQueueLengthAfterTran(indexOfTran, indexOfState,
						computeProbForTotalQueueLength_Mid(indexOfPeriod, indexOfTran, indexOfState));
				}
			}
		}
		else if (type == "Low")
		{
			for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInUniform; indexOfTran++) {
				for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {
					setProbOfTotalQueueLengthAfterTran(indexOfTran, indexOfState,
						computeProbForTotalQueueLength_Low(indexOfPeriod, indexOfTran, indexOfState));
				}
			}
		}
		else { cout << "error xxx" << endl; system("pause"); }

	}
};

class ComputeUsingOldFunctions {
public:
	double computeA_ServiceCapcity(const SolutionOfServers solutionOfServers, double startTime) {
		double endTime = startTime + torWaitingTime;
		int firstPeriod = int(startTime) + 1;
		int lastPeriod = int(endTime) + 1;
		double sum = 0;
		for (int indexOfPeriod = firstPeriod; indexOfPeriod <= lastPeriod; indexOfPeriod++) {
			int tempIndexOfPeriod = indexOfPeriod;
			while (tempIndexOfPeriod > NumberOfPeriodsCondsidered)tempIndexOfPeriod = tempIndexOfPeriod - NumberOfPeriodsCondsidered;
			sum += (min(endTime, double(indexOfPeriod)) - max(startTime, double(indexOfPeriod) - 1.0))
				* double(MiuForEachPeriod[tempIndexOfPeriod]) * double(solutionOfServers.getSolutionOfServersForTime(tempIndexOfPeriod));
		}
		return sum;
	}

	int computeXi_CustomerNeedToBeServed(int q, const SolutionOfServers solutionOfServers, double startTime, Shifts shiftEnv, int offshift[]) {
		double endTime = startTime + torWaitingTime;
		int firstPeriod = int(startTime) + 1;
		int lastPeriod = int(endTime) + 1;
		int tempLastPeriod = lastPeriod;
		while (tempLastPeriod > NumberOfPeriodsCondsidered)tempLastPeriod = tempLastPeriod - NumberOfPeriodsCondsidered;

		int numberOfServersInTheLastPeriod = solutionOfServers.getSolutionOfServersForTime(tempLastPeriod);
		int sumOfLeaveServers = 0;
		for (int indexOfPeriod = firstPeriod; indexOfPeriod <= lastPeriod - 1; indexOfPeriod++) {
			while (indexOfPeriod > NumberOfPeriodsCondsidered)indexOfPeriod -= NumberOfPeriodsCondsidered;
			sumOfLeaveServers += offshift[indexOfPeriod];
		}
		return max(q - numberOfServersInTheLastPeriod - sumOfLeaveServers, -1);
	}
};

//historyInThisPeriod   P_iq
double totaltimeOfHyper4d = 0;
class ComputationOfServicelLevel {
public:
	int offShift[NumberOfPeriodsCondsidered + 1] = { 0 };
	//UniformProcess uniformProcess;
public:
		double computeServiceLevelForPeriod_H(const int indexOfPeriod, const SolutionOfServers solutionOfServers, const Shifts shiftEnv) {
				try {
			HistoryOfVectorOfUniformState historyInThisPeriod;//P_iq
			historyInThisPeriod.initialHistoryOfVectorOfUniformState(indexOfPeriod);//(indexOfPeriod, uniformProcess.getVectorOfUniformStateInStart(indexOfPeriod), solutionOfServers);
						

			//double sumOfUnServed = 0;

			double sumOfUnServedWeighted = 0;
			double oldsumtest = 0;
			double gama = double(solutionOfServers.getSolutionOfServersForTime(indexOfPeriod)) * double(MiuForEachPeriod[indexOfPeriod]) + Lamda[indexOfPeriod] + LamdaLP[indexOfPeriod];

			ComputeUsingOldFunctions CUOF;

			

			//double thetaWaitingTime = torWaitingTime - double(int(torWaitingTime));
			//double valueSet[2] = { (0.5 - thetaWaitingTime / 2.0), (1.0 - thetaWaitingTime / 2.0) };
			//double gapSet[2] = { 1 - thetaWaitingTime,thetaWaitingTime };

			//for (int indexOfT = 0; indexOfT < 2; indexOfT++) {
			//	if (abs(thetaWaitingTime) < BreakThreshold && indexOfT == 1)continue;

			//	double valueOfT = double(indexOfPeriod) - 1.0 + valueSet[indexOfT];

			//	double A_t = CUOF.computeA_ServiceCapcity(solutionOfServers, valueOfT);

			//	double sumOfUnServedOld = sumOfUnServed;

			//	double sumOfUnservedTemp = 0;

			//	for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInIntegral; indexOfState++) {
			//		int Xi = CUOF.computeXi_CustomerNeedToBeServed(indexOfState, solutionOfServers, valueOfT, shiftEnv, offShift);

			//		double fie;
			//		if (Xi < 0)
			//			fie = 0;
			//		else
			//			fie = util.getPoissonAccum(A_t, Xi);

			//		double sumOfPF = 0;
			//		for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInIntegral; indexOfTran++) {
			//			double valueOfPF = (historyInThisPeriod.getProbOfHighProrityQueueLengthAfterTran(indexOfTran, indexOfState)
			//				* util.getPoisson(gama * (valueOfT - (double(indexOfPeriod) - 1.0)), indexOfTran));
			//			sumOfPF += valueOfPF;
			//			//if (valueOfPF <= BreakThreshold)break;
			//		}
			//		//cout << "indexofstate:\t" << indexOfPeriod << "\t" << fie * sumOfPF << endl;
			//		sumOfUnServed += fie * sumOfPF;

			//		sumOfUnservedTemp += fie * sumOfPF;
			//		//if (fie * sumOfPF < BreakThreshold)break;
			//	}
			//	if (indexOfPeriod == 1 || indexOfPeriod == 2)
			//		cout << valueOfT << "," << 1 - (sumOfUnServed - sumOfUnServedOld) << endl;
			//	sumOfUnServedWeighted += sumOfUnservedTemp * gapSet[indexOfT];
			//}
			double thetaWaitingTime = torWaitingTime - double(int(torWaitingTime));
			double gapSet[2] = { 1 - thetaWaitingTime,thetaWaitingTime };

			double SL1_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
			double SL1_AccForBataOfState[MaxValueOfStatesInIntegral + 1] = { 0 };
			double SL2_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
			double SL2_AccForBataOfState[MaxValueOfStatesInIntegral + 1] = { 0 };
			int countSL1Acc = 0;
			int countSL2Acc = 0;

			for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++) 			{
				double valueOfT = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod);
				
				if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTime))countSL1Acc++;
				else countSL2Acc++;

				double A_t = CUOF.computeA_ServiceCapcity(solutionOfServers, valueOfT);
				//double sumOfUnServedOld = sumOfUnServed;
				for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInIntegral; indexOfState++) {
					int Xi = CUOF.computeXi_CustomerNeedToBeServed(indexOfState, solutionOfServers, valueOfT, shiftEnv, offShift);
					double fie;
					if (Xi < 0)
						fie = 0;
					else
						fie = util.getPoissonAccum(A_t, Xi);
					double sumOfPF = 0;
					for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInIntegral; indexOfTran++) {
						double valueOfPF = (historyInThisPeriod.getProbOfHighProrityQueueLengthAfterTran(indexOfTran, indexOfState)
							* util.getPoisson(gama * (valueOfT - (double(indexOfPeriod) - 1.0)), indexOfTran));
						
						sumOfPF += valueOfPF;
						//if (valueOfPF <= BreakThreshold)break;
					}
					//cout << "indexofstate:\t" << indexOfPeriod << "\t" << fie * sumOfPF << endl;
					
										//sumOfUnServed += fie * sumOfPF;

					if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTime)) {
						SL1_AccForStateProb[indexOfState] += sumOfPF;
						SL1_AccForBataOfState[indexOfState] += fie;
					}
					else
					{
						SL2_AccForStateProb[indexOfState] += sumOfPF;
						SL2_AccForBataOfState[indexOfState] += fie;
					}
					//if (fie * sumOfPF < BreakThreshold)break;
				}
				/*if (indexOfPeriod == 1|| indexOfPeriod == 2)
				cout << valueOfT <<","<< 1-(sumOfUnServed-sumOfUnServedOld) << endl;*/
			}

			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInIntegral; indexOfState++) {
				SL1_AccForStateProb[indexOfState] = SL1_AccForStateProb[indexOfState] / double(countSL1Acc);
				SL1_AccForBataOfState[indexOfState] = SL1_AccForBataOfState[indexOfState] / double(countSL1Acc);
				SL2_AccForStateProb[indexOfState] = SL2_AccForStateProb[indexOfState] / double(countSL2Acc);
				SL2_AccForBataOfState[indexOfState] = SL2_AccForBataOfState[indexOfState] / double(countSL2Acc);
			}
			double sumOfUnServed_SL1 = 0;
			double sumOfUnServed_SL2 = 0;
			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInIntegral; indexOfState++) {
				sumOfUnServed_SL1 += SL1_AccForStateProb[indexOfState] * SL1_AccForBataOfState[indexOfState];
				sumOfUnServed_SL2 += SL2_AccForStateProb[indexOfState] * SL2_AccForBataOfState[indexOfState];
			}
			sumOfUnServedWeighted = sumOfUnServed_SL1 * gapSet[0] + sumOfUnServed_SL2 * gapSet[1];
			cout << indexOfPeriod << "SL1: " << 1.0 - sumOfUnServed_SL1 << " ,SL2: " << 1.0 - sumOfUnServed_SL2 << " ,weigthed sum:" << 1.0 - sumOfUnServedWeighted << endl;

			//sumOfUnServed = sumOfUnServed / double(PrecisionOfIntegralOfTInOnePeriod);
			return (1.0 - sumOfUnServedWeighted);
		}
		catch (char* e)		{
			cout << e << endl;
		}

		return -1;

	}
	
		int convertStateToStateOfWaitingAnd1(int totalCustomersInSystem, int indexOfPeriod, SolutionOfServers solutionOfServers) {
		int numOfServersNow = solutionOfServers.getSolutionOfServersForTime(indexOfPeriod);
		return max(0, (totalCustomersInSystem - numOfServersNow + 1));
	}

	double computeInitialProbOfOneStateForSomeTime(int indexOfState, double timeInThisPeriod, int indexOfPeriod, SolutionOfServers solutionOfServers, HistoryOfVectorOfUniformState& historyOfVectorOfUniformState) {
		double gama = double(solutionOfServers.getSolutionOfServersForTime(indexOfPeriod)) * double(MiuForEachPeriod[indexOfPeriod]) + Lamda[indexOfPeriod] + LamdaLP[indexOfPeriod];
		double gamaMultiT = gama * timeInThisPeriod;
		double sumProbOfState = 0;
		for (int indexOfTran = 0; indexOfTran <= MaxNumberOfTranstionTimesInIntegral; indexOfTran++) {
			sumProbOfState += util.getPoisson(gamaMultiT, indexOfTran) * historyOfVectorOfUniformState.getProbOfTotalQueueLengthAfterTran(indexOfTran, indexOfState);
		}
		return sumProbOfState;
	}

	double API_computeSLForLowPriorityCustomer(int indexOfPeriod, SolutionOfServers solutionOfServers, HistoryOfVectorOfUniformState& historyOfVectorOfUniformState) {
		if (ApproIndex == 0) {
			return API_computeSLForLowPriorityCustomer_old(indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);
		}
		else if (ApproIndex == 1) {
			return API_computeSLForLowPriorityCustomer_Appr1(indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);
		}
		else if (ApproIndex == 2) {
			return API_computeSLForLowPriorityCustomer_Appr2(indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);
		}
		else {
			cout << "wrong dhajkh " << endl;
			system("pause");
			return 0;
		}
	}

	double API_computeSLForLowPriorityCustomer_old(int indexOfPeriod, SolutionOfServers solutionOfServers, HistoryOfVectorOfUniformState & historyOfVectorOfUniformState) {
		long double sumOfIntegral = 0;
//#pragma omp parallel for num_threads(2*numProcs-1)
		
		
		{
			//for (int indexOfT = 0; indexOfT < 2; indexOfT++)
			//{
			//	if (abs(thetaWaitingTimePL) < BreakThreshold && indexOfT == 1)continue;
			//	double arrivingTime = double(indexOfPeriod) - 1.0 + valueSet[indexOfT] + 1e-9;

			//	double timeInThisPeriod = arrivingTime - double(int(arrivingTime));
			//	VectorOfStateProb vectorOfStateProb;
			//	vectorOfStateProb.initialProbOfState();
			//	for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {

			//		int API_state_waitingAnd1 = convertStateToStateOfWaitingAnd1(indexOfState, indexOfPeriod, solutionOfServers);

			//		double probOfThisSystemState = computeInitialProbOfOneStateForSomeTime(indexOfState, timeInThisPeriod, indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);

			//		vectorOfStateProb.setProbOfState(API_state_waitingAnd1,
			//			vectorOfStateProb.getProbOfState(API_state_waitingAnd1) + probOfThisSystemState);
			//	}

			//	//if (abs(vectorOfStateProb.confirmSumOfProb() - 1) > EPSILON) {
						//	//	cout << "indexOfPeriod\t" << indexOfPeriod << "\t" << vectorOfStateProb.confirmSumOfProb() << endl;
			//	//}
			//	//cout << endl;
			//	//vectorOfStateProb.showVectorUtill(20);

						//	if (IsUsingEigen == true) {
			//		vectorOfStateProb = vectorOfStateProb.computeByEigen(arrivingTime, vectorOfStateProb);
			//	}
			//	else
			//	{
			//		vectorOfStateProb = vectorOfStateProb.stepByStep(arrivingTime, vectorOfStateProb);
			//	}




			//	double tempsum = vectorOfStateProb.getSLByVectorForArringAtTimeT();
			//	sumOfIntegral += tempsum * gapSet[indexOfT];

			//	if (sumOfIntegral == 0) {
			//		cout << "error\t" << indexOfPeriod << endl;
			//		vectorOfStateProb.showVectorUtill(20);
			//		//system("pause");
			//	}

			//	/*if (abs(vectorOfStateProb.confirmSumOfProb() - 1) > EPSILON) {
						//		cout << "indexOfPeriod\t" << indexOfPeriod << "\t" << vectorOfStateProb.confirmSumOfProb() << endl;
			//	}*/
			//	/*if (indexOfPeriod == 10)
			//		cout << arrivingTime << "\t"<< vectorOfStateProb.getSLByVectorForArringAtTimeT() / double(PrecisionOfIntegralOfTInOnePeriod) <<endl;*/
			//	if (indexOfPeriod == 1 || indexOfPeriod == 2)
			//		cout << arrivingTime << "," << tempsum << endl;




			//}
			//sumOfIntegral = sumOfIntegral / double(PrecisionOfIntegralOfTInOnePeriod);
		}

		double thetaWaitingTimePL = torWaitingTimePL - double(int(torWaitingTimePL));
		double gapSet[2] = { 1 - thetaWaitingTimePL,thetaWaitingTimePL };

		double SL1_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		double SL2_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		int countSL1Acc = 0;
		int countSL2Acc = 0;

		double sumOfSLWeighted = 0;

		
		for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++)
		{

			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;

			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL))countSL1Acc++;
			else countSL2Acc++;

			double timeInThisPeriod = arrivingTime - double(int(arrivingTime));
			//VectorOfStateProb vectorOfStateProb;
			//vectorOfStateProb.initialProbOfState();
			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {

				int API_state_waitingAnd1 = convertStateToStateOfWaitingAnd1(indexOfState, indexOfPeriod, solutionOfServers);

				double probOfThisSystemState = computeInitialProbOfOneStateForSomeTime(indexOfState, timeInThisPeriod, indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);

				//vectorOfStateProb.setProbOfState(API_state_waitingAnd1,
				//	vectorOfStateProb.getProbOfState(API_state_waitingAnd1) + probOfThisSystemState);


				if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
					SL1_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;				}
				else
				{
					SL2_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;
				}
			}

			//if (abs(vectorOfStateProb.confirmSumOfProb() - 1) > EPSILON) {
						//	cout << "indexOfPeriod\t" << indexOfPeriod << "\t" << vectorOfStateProb.confirmSumOfProb() << endl;
			//}
			//cout << endl;
			//vectorOfStateProb.showVectorUtill(20);
		}

		for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
			SL1_AccForStateProb[API_state_waitingAnd1] = SL1_AccForStateProb[API_state_waitingAnd1] / double(countSL1Acc);
			SL2_AccForStateProb[API_state_waitingAnd1] = SL2_AccForStateProb[API_state_waitingAnd1] / double(countSL2Acc);
		}
		
		double SL1_AccForBataOfTheStateVector =  0 ;
		double SL2_AccForBataOfTheStateVector =  0 ;


				for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++){
			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;

			double out = 0;
			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
				VectorOfStateProb vectorOfStateProbSL1;
				vectorOfStateProbSL1.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL1.setProbOfState(API_state_waitingAnd1, SL1_AccForStateProb[API_state_waitingAnd1]);
				}
				if (IsUsingEigen == true) {
					vectorOfStateProbSL1 = vectorOfStateProbSL1.computeByEigen(arrivingTime, vectorOfStateProbSL1);
				}
				else
				{
					vectorOfStateProbSL1 = vectorOfStateProbSL1.stepByStep(arrivingTime, vectorOfStateProbSL1);
				}
				double tempsum = vectorOfStateProbSL1.getSLByVectorForArringAtTimeT();
				out = tempsum;
				SL1_AccForBataOfTheStateVector += tempsum;
			}
			else if (abs(thetaWaitingTimePL) >= BreakThreshold) {
				VectorOfStateProb vectorOfStateProbSL2;
				vectorOfStateProbSL2.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL2.setProbOfState(API_state_waitingAnd1, SL2_AccForStateProb[API_state_waitingAnd1]);
				}
				if (IsUsingEigen == true) {
					vectorOfStateProbSL2 = vectorOfStateProbSL2.computeByEigen(arrivingTime, vectorOfStateProbSL2);
				}
				else
				{
					vectorOfStateProbSL2 = vectorOfStateProbSL2.stepByStep(arrivingTime, vectorOfStateProbSL2);
				}
				double tempsum = vectorOfStateProbSL2.getSLByVectorForArringAtTimeT();
				out = tempsum;

				SL2_AccForBataOfTheStateVector += tempsum;
			}
		
				cout <<"***peroid: "<< indexOfPeriod << " ***time " << arrivingTime << " service level: " << out << endl;
			/*if (abs(vectorOfStateProb.confirmSumOfProb() - 1) > EPSILON) {
				cout << "out ": error sum is not 1 LP_Customer\t" << endl;
				cout << "indexOfPeriod\t" << indexOfPeriod << "\t" << vectorOfStateProb.confirmSumOfProb() << endl;
			}*/
			/*if (indexOfPeriod == 10)
				cout << arrivingTime << "\t"<< vectorOfStateProb.getSLByVectorForArringAtTimeT() / double(PrecisionOfIntegralOfTInOnePeriod) <<endl;*/

		}
		SL1_AccForBataOfTheStateVector = SL1_AccForBataOfTheStateVector / double(countSL1Acc);
		SL2_AccForBataOfTheStateVector = SL2_AccForBataOfTheStateVector / double(countSL2Acc);


		if (abs(thetaWaitingTimePL) >= BreakThreshold) {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0] + SL2_AccForBataOfTheStateVector * gapSet[1];
			cout << indexOfPeriod << "£¨Lower£©SL1: " << SL1_AccForBataOfTheStateVector << " ,SL2: " << SL2_AccForBataOfTheStateVector << " ,weighted sum:" << sumOfSLWeighted << endl;
		}
		else {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0];
			cout << indexOfPeriod << "£¨Lower£©SL1: " << SL1_AccForBataOfTheStateVector << " ,no SL2, weighted sum:" << sumOfSLWeighted << endl;
		}

		//cout << " SL:\t" << sumOfIntegral << endl;
		
		return sumOfSLWeighted;
	}

	double API_computeSLForLowPriorityCustomer_Appr1(int indexOfPeriod, SolutionOfServers solutionOfServers, HistoryOfVectorOfUniformState& historyOfVectorOfUniformState) {
		long double sumOfIntegral = 0;

		double thetaWaitingTimePL = torWaitingTimePL - double(int(torWaitingTimePL));
		double gapSet[2] = { 1 - thetaWaitingTimePL,thetaWaitingTimePL };

		double SL1_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		double SL2_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		int countSL1Acc = 0;
		int countSL2Acc = 0;

		double sumOfSLWeighted = 0;


		for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++)
		{

			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;

			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL))countSL1Acc++;
			else countSL2Acc++;

			double timeInThisPeriod = arrivingTime - double(int(arrivingTime));
			//VectorOfStateProb vectorOfStateProb;
			//vectorOfStateProb.initialProbOfState();
			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {

				int API_state_waitingAnd1 = convertStateToStateOfWaitingAnd1(indexOfState, indexOfPeriod, solutionOfServers);

				double probOfThisSystemState = computeInitialProbOfOneStateForSomeTime(indexOfState, timeInThisPeriod, indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);

				//vectorOfStateProb.setProbOfState(API_state_waitingAnd1,
				//	vectorOfStateProb.getProbOfState(API_state_waitingAnd1) + probOfThisSystemState);


				if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
					SL1_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;				}
				else
				{
					SL2_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;
				}
			}

		}

		for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
			SL1_AccForStateProb[API_state_waitingAnd1] = SL1_AccForStateProb[API_state_waitingAnd1] / double(countSL1Acc);
			SL2_AccForStateProb[API_state_waitingAnd1] = SL2_AccForStateProb[API_state_waitingAnd1] / double(countSL2Acc);
		}

		double SL1_AccForBataOfTheStateVector = 0;
		double SL2_AccForBataOfTheStateVector = 0;

				int right = NumberOfPiecesInAppro / 2;
		int left = NumberOfPiecesInAppro - right;

		for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++) {
			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;
			cout << "origon arrival time: " << arrivingTime << " ";


			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
				
				int indexNew = ceil( 1e-9+
					double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) /
					((1.0 - thetaWaitingTimePL) / double(left)));				if (indexOfT == 0)indexNew = 1;

				arrivingTime = double(indexOfPeriod) - 1.0 + (1.0 - thetaWaitingTimePL) * double(2 * indexNew - 1) / (2.0 * left)+ 1e-9;
				cout << "renew arrival time (left): " << arrivingTime << endl;

				VectorOfStateProb vectorOfStateProbSL1;
				vectorOfStateProbSL1.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL1.setProbOfState(API_state_waitingAnd1, SL1_AccForStateProb[API_state_waitingAnd1]);
				}
				if (IsUsingEigen == true) {
					vectorOfStateProbSL1 = vectorOfStateProbSL1.computeByEigen(arrivingTime, vectorOfStateProbSL1);
				}
				else
				{
					vectorOfStateProbSL1 = vectorOfStateProbSL1.stepByStep(arrivingTime, vectorOfStateProbSL1);
				}
				double tempsum = vectorOfStateProbSL1.getSLByVectorForArringAtTimeT();
				SL1_AccForBataOfTheStateVector += tempsum;
			}
			else if (abs(thetaWaitingTimePL) >= BreakThreshold) {
				
				int indexNew = ceil(1e-9 +
					(double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) - (1.0 - thetaWaitingTimePL)) /
					(double(thetaWaitingTimePL) / double(right)));				if (indexOfT == 0)indexNew = 1;
				cout << " indexofNew: " << indexNew << "  ";
				arrivingTime = double(indexOfPeriod) - 1.0 + (1 - thetaWaitingTimePL) + (thetaWaitingTimePL) * double(2 * indexNew - 1) / (2.0 * right) + 1e-9;
				cout << "renew arrival time(right): " << arrivingTime << endl;

				VectorOfStateProb vectorOfStateProbSL2;
				vectorOfStateProbSL2.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL2.setProbOfState(API_state_waitingAnd1, SL2_AccForStateProb[API_state_waitingAnd1]);
				}
				if (IsUsingEigen == true) {
					vectorOfStateProbSL2 = vectorOfStateProbSL2.computeByEigen(arrivingTime, vectorOfStateProbSL2);
				}
				else
				{
					vectorOfStateProbSL2 = vectorOfStateProbSL2.stepByStep(arrivingTime, vectorOfStateProbSL2);
				}
				double tempsum = vectorOfStateProbSL2.getSLByVectorForArringAtTimeT();
				SL2_AccForBataOfTheStateVector += tempsum;
			}


			/*if (abs(vectorOfStateProb.confirmSumOfProb() - 1) > EPSILON) {
				cout << "out : error not 1 LP_Customer\t" << endl;
				cout << "indexOfPeriod\t" << indexOfPeriod << "\t" << vectorOfStateProb.confirmSumOfProb() << endl;
			}*/
			/*if (indexOfPeriod == 10)
				cout << arrivingTime << "\t"<< vectorOfStateProb.getSLByVectorForArringAtTimeT() / double(PrecisionOfIntegralOfTInOnePeriod) <<endl;*/

		}
		SL1_AccForBataOfTheStateVector = SL1_AccForBataOfTheStateVector / double(countSL1Acc);
		SL2_AccForBataOfTheStateVector = SL2_AccForBataOfTheStateVector / double(countSL2Acc);


		if (abs(thetaWaitingTimePL) >= BreakThreshold) {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0] + SL2_AccForBataOfTheStateVector * gapSet[1];
			cout << indexOfPeriod << "Lower SL1: " << SL1_AccForBataOfTheStateVector << " ,SL2: " << SL2_AccForBataOfTheStateVector << " ,weigthed sum:" << sumOfSLWeighted << endl;
		}
		else {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0];
			cout << indexOfPeriod << "Lower SL1: " << SL1_AccForBataOfTheStateVector << " ,no SL2£¬weigthed sum:" << sumOfSLWeighted << endl;
		}

		//cout << " SL:\t" << sumOfIntegral << endl;

		return sumOfSLWeighted;
	}

	double API_computeSLForLowPriorityCustomer_Appr2(int indexOfPeriod, SolutionOfServers solutionOfServers, HistoryOfVectorOfUniformState& historyOfVectorOfUniformState) {
		long double sumOfIntegral = 0;

		double thetaWaitingTimePL = torWaitingTimePL - double(int(torWaitingTimePL));
		double gapSet[2] = { 1 - thetaWaitingTimePL,thetaWaitingTimePL };

		double SL1_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		double SL2_AccForStateProb[MaxValueOfStatesInIntegral + 1] = { 0 };
		int countSL1Acc = 0;
		int countSL2Acc = 0;

		double sumOfSLWeighted = 0;


		for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++)
		{

			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;

			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL))countSL1Acc++;
			else countSL2Acc++;

			double timeInThisPeriod = arrivingTime - double(int(arrivingTime));
			//VectorOfStateProb vectorOfStateProb;
			//vectorOfStateProb.initialProbOfState();
			for (int indexOfState = 0; indexOfState <= MaxValueOfStatesInUniform; indexOfState++) {

				int API_state_waitingAnd1 = convertStateToStateOfWaitingAnd1(indexOfState, indexOfPeriod, solutionOfServers);

				double probOfThisSystemState = computeInitialProbOfOneStateForSomeTime(indexOfState, timeInThisPeriod, indexOfPeriod, solutionOfServers, historyOfVectorOfUniformState);

				//vectorOfStateProb.setProbOfState(API_state_waitingAnd1,
				//	vectorOfStateProb.getProbOfState(API_state_waitingAnd1) + probOfThisSystemState);


				if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
					SL1_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;				}
				else
				{
					SL2_AccForStateProb[API_state_waitingAnd1] += probOfThisSystemState;
				}
			}

		}

		for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
			SL1_AccForStateProb[API_state_waitingAnd1] = SL1_AccForStateProb[API_state_waitingAnd1] / double(countSL1Acc);
			SL2_AccForStateProb[API_state_waitingAnd1] = SL2_AccForStateProb[API_state_waitingAnd1] / double(countSL2Acc);
		}

		double SL1_AccForBataOfTheStateVector = 0;
		double SL2_AccForBataOfTheStateVector = 0;

				int right = NumberOfPiecesInAppro / 2;
		int left = NumberOfPiecesInAppro - right;

		for (int indexOfT = 0; indexOfT < PrecisionOfIntegralOfTInOnePeriod; indexOfT++) {
			double arrivingTime = double(indexOfPeriod) - 1.0 + double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) + 1e-9;


			if (double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) < (1.0 - thetaWaitingTimePL)) {
								//arrivingTime = double(indexOfPeriod) - 1.0 + (1.0 / 2.0 - thetaWaitingTimePL / 2.0) + 1e-9;
				int indexNew = ceil(1e-9+
					double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) /
					((1.0 - thetaWaitingTimePL) / double(left)));				if (indexOfT == 0)indexNew = 1;

				arrivingTime = double(indexOfPeriod) - 1.0 + (1.0 - thetaWaitingTimePL) * double(2.0 * indexNew - 1) / (2.0 * left) + 1e-9;


				VectorOfStateProb vectorOfStateProbSL1;
				vectorOfStateProbSL1.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL1.setProbOfState(API_state_waitingAnd1, SL1_AccForStateProb[API_state_waitingAnd1]);
				}
				double tempsum = 0;
				if (IsUsingEigen == true) {
					vectorOfStateProbSL1 = vectorOfStateProbSL1.computeByEigen(arrivingTime, vectorOfStateProbSL1);
				}
				else
				{
					double sum = 0;
					int count = 10;
					for (int i = 1; i <= count; i++) {
						double remainT = torWaitingTimePL;
						double var = min(torWaitingTimePL, (1 - thetaWaitingTimePL )/ 2.0 / double(left));//length
						remainT += ((i - (double)(1 + count) / 2.0) / ((double)(1 + count) / 2.0)) * var;
						sum += vectorOfStateProbSL1.stepByStep(arrivingTime, vectorOfStateProbSL1, remainT).getSLByVectorForArringAtTimeT();


					}
					tempsum = sum / (double)count;
				}
				//double tempsum = vectorOfStateProbSL1.getSLByVectorForArringAtTimeT();
				SL1_AccForBataOfTheStateVector += tempsum;
			}
			else if (abs(thetaWaitingTimePL) >= BreakThreshold) {
								//arrivingTime = double(indexOfPeriod) - 1.0 + (1 - thetaWaitingTimePL / 2) + 1e-9;

				int indexNew = ceil(1e-9 +
					(double(indexOfT) / double(PrecisionOfIntegralOfTInOnePeriod) - (1.0 - thetaWaitingTimePL)) /
					(double(thetaWaitingTimePL) / double(right)));				if (indexOfT == 0)indexNew = 1;

				arrivingTime = double(indexOfPeriod) - 1.0 + (1 - thetaWaitingTimePL) + (thetaWaitingTimePL) * double(2.0 * indexNew - 1) / (2.0 * right) + 1e-9;


				VectorOfStateProb vectorOfStateProbSL2;
				vectorOfStateProbSL2.initialProbOfState();
				for (int API_state_waitingAnd1 = 0; API_state_waitingAnd1 <= MaxValueOfStatesInUniform; API_state_waitingAnd1++) {
					vectorOfStateProbSL2.setProbOfState(API_state_waitingAnd1, SL2_AccForStateProb[API_state_waitingAnd1]);
				}
				double tempsum = 0;
				if (IsUsingEigen == true) {
					vectorOfStateProbSL2 = vectorOfStateProbSL2.computeByEigen(arrivingTime, vectorOfStateProbSL2);
				}
				else
				{
					double sum = 0;
					int count = 10;
					for (int i = 1; i <= count; i++) {
						double remainT = torWaitingTimePL;
						double var = min(torWaitingTimePL, (thetaWaitingTimePL / 2.0)/double(right));//length
						remainT += ((i - (double)(1 + count) / 2.0) / ((double)(1 + count) / 2.0)) * var;
						sum += vectorOfStateProbSL2.stepByStep(arrivingTime, vectorOfStateProbSL2, remainT).getSLByVectorForArringAtTimeT();
					}
					tempsum = sum / (double)count;
				}
				SL2_AccForBataOfTheStateVector += tempsum;
			}


		}
		SL1_AccForBataOfTheStateVector = SL1_AccForBataOfTheStateVector / double(countSL1Acc);
		SL2_AccForBataOfTheStateVector = SL2_AccForBataOfTheStateVector / double(countSL2Acc);


		if (abs(thetaWaitingTimePL) >= BreakThreshold) {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0] + SL2_AccForBataOfTheStateVector * gapSet[1];
			cout << indexOfPeriod << "Lower SL1: " << SL1_AccForBataOfTheStateVector << " ,SL2: " << SL2_AccForBataOfTheStateVector << " ,weigthed sum:" << sumOfSLWeighted << endl;
		}
		else {
			sumOfSLWeighted = SL1_AccForBataOfTheStateVector * gapSet[0];
			cout << indexOfPeriod << "Lower SL1: " << SL1_AccForBataOfTheStateVector << " ,no SL2£¬weigthed sum:" << sumOfSLWeighted << endl;
		}

		//cout << " SL:\t" << sumOfIntegral << endl;

		return sumOfSLWeighted;
	}


	double computeServiceLevelForPeriod_L(int indexOfPeriod, const SolutionOfServers solutionOfServers, const Shifts shiftEnv) {
		try {
			HistoryOfVectorOfUniformState historyInThisPeriod;//P_iq
			historyInThisPeriod.initialHistoryOfVectorOfUniformState(indexOfPeriod);//(indexOfPeriod, uniformProcess.getVectorOfUniformStateInStart(indexOfPeriod), solutionOfServers);
						
			double SLinThisPeriod = API_computeSLForLowPriorityCustomer(indexOfPeriod, solutionOfServers, historyInThisPeriod);

			return SLinThisPeriod;
		}
		catch (char* e) {
			cout << e << endl;
		}

		return -1;

	}

	double computeServiceLevelForPeriod_L(int indexOfPeriod, const SolutionOfServers solutionOfServers, const Shifts shiftEnv, string type) {
		try {
			HistoryOfVectorOfUniformState historyInThisPeriod;//P_iq
			historyInThisPeriod.initialHistoryOfVectorOfUniformState(indexOfPeriod, type);//(indexOfPeriod, uniformProcess.getVectorOfUniformStateInStart(indexOfPeriod), solutionOfServers);
			
			double SLinThisPeriod = API_computeSLForLowPriorityCustomer(indexOfPeriod, solutionOfServers, historyInThisPeriod);

			return SLinThisPeriod;
		}
		catch (char* e) {
			cout << e << endl;
		}

		return -1;

	}


	int computeUpShiftServerAtTheStartOfOnePeriod(int indexOfPeriod, const SolutionOnShift solutionOnShift, const Shifts shiftEnv) {
		int countOfShifts = 0;
		for (int indexOfShift = 1; indexOfShift <= shiftEnv.totalNumOfShifts; indexOfShift++) {
			if (shiftEnv.shifts[indexOfShift].startTime == indexOfPeriod - 1) {
				countOfShifts += solutionOnShift.getServerOnShift(indexOfShift);
			}
		}
		return countOfShifts;
	}

	void initialComputationOfServicelLevel(const SolutionOnShift solutionOnShift, const SolutionOfServers solutionOfServers, const Shifts shiftEnv) {
		try {
			util.initial();
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				offShift[indexOfPeriod] = computeOffShiftInPeriod(shiftEnv, solutionOnShift, indexOfPeriod);
			}

						
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)cout << offShift[indexOfPeriod] << "\t";
			cout << endl;
			
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)cout << solutionOfServers.getSolutionOfServersForTime(indexOfPeriod) << "\t";
			cout << endl;

			Uniformization_priority_FXY(solutionOnShift, solutionOfServers, shiftEnv);

			
						for (int indexOfPeriod = 1; indexOfPeriod <= MaxNumOfPeriod; indexOfPeriod++) {
				VectorOfStateProb::API_numOfServers(indexOfPeriod, solutionOfServers.getSolutionOfServersForTime(indexOfPeriod));
				VectorOfStateProb::API_numOfOnShifts(indexOfPeriod, computeUpShiftServerAtTheStartOfOnePeriod(indexOfPeriod, solutionOnShift, shiftEnv));
			}

		}
		catch (char * e) {
			cout << e << endl;
		}
	}

	void initialComputationOfServicelLevel_3Type(const SolutionOnShift solutionOnShift, const SolutionOfServers solutionOfServers, const Shifts shiftEnv) {
		try {
			util.initial();
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
				offShift[indexOfPeriod] = computeOffShiftInPeriod(shiftEnv, solutionOnShift, indexOfPeriod);
			}

						
			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)cout << offShift[indexOfPeriod] << "\t";
			cout << endl;

			for (int indexOfPeriod = 1; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++)cout << solutionOfServers.getSolutionOfServersForTime(indexOfPeriod) << "\t";
			cout << endl;

			//Uniformization_priority_FXY(solutionOnShift, solutionOfServers, shiftEnv);//**************************************todo
			Uniformization_priority_FXY_forthreetype(solutionOnShift, solutionOfServers, shiftEnv);//**************************************todo


						
						for (int indexOfPeriod = 1; indexOfPeriod <= MaxNumOfPeriod; indexOfPeriod++) {
				VectorOfStateProb::API_numOfServers(indexOfPeriod, solutionOfServers.getSolutionOfServersForTime(indexOfPeriod));
				VectorOfStateProb::API_numOfOnShifts(indexOfPeriod, computeUpShiftServerAtTheStartOfOnePeriod(indexOfPeriod, solutionOnShift, shiftEnv));
			}
		}
		catch (char* e) {
			cout << e << endl;
		}
	}


	//***************************************************FXY**********************************************
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x)< (y)) ? (x) : (y))
	double fac(int i)
	{
		if (i <= 1)
			return 1;
		else if (i == 2)
			return 2;
		else if (i == 3)
			return 6;
		else if (i == 4)
			return 24;
		else
			return fac(i - 1) * i;
	}
		double HyperGeo(int k, int n, int q, int N_d)	{
		if (q == 0 || n == 0)
			return 0;
		else if (N_d - q - n + k < 0 || q < k)
			return 0;
		else
		{
			return fac(q) * fac(N_d - q) * fac(n) * fac(N_d - n) / fac(k) / fac(q - k) / fac(n - k) / fac(N_d - q - n + k) / fac(N_d);
		}
	}

	long factorial(long number) {
		if (number <= 1)
			return 1;
		else
			return number * factorial(number - 1);
	}

	int combinator(int n, int m) {
		int temp;
		if (n < m) {
			temp = n;
			n = m;
			m = temp;
		}
		return factorial(n) / (factorial(m) * factorial(n - m));
	}

	double fai(double a, int b)//calculate (a^b)/(b!)
	{
		double result = 1;
		for (int i = 1; i <= b; i++) result *= (a / i);
		return result;
	}

	double HyperGeo_fan(int ddeta, int u, int q, int s)	{
		return util.getHyperGeometric(ddeta, u, s, q);


		//return double(combinator(q, ddeta)) * double(combinator(s - q, u - ddeta)) / double(combinator(s, u));
	}



	void Uniformization_priority_FXY(const SolutionOnShift solutionOnShift, const SolutionOfServers solutionOfServers, const Shifts shiftEnv)
	{
		//double FXY_uniform_Piq[NumberOfPeriodsCondsidered][MaxNumberOfTranstionTimesInUniform + 1] [MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1];
		FXY_uniform_Piq = vector<vector<vector<vector<double>>>>(NumberOfPeriodsCondsidered, vector<vector<vector<double>>>(MaxNumberOfTranstionTimesInUniform + 1,
			vector<vector<double>>(MaxValueOfStatesInUniform + 1, vector<double>(MaxValueOfStatesInUniform + 1, 0))));


		//#define muu(x) (MIN(p[t], x)*mu)
		long double pi[MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1]{},
			npi[MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1]{},
			newPi[MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1]{};		long double w[NumberOfPeriodsCondsidered], u[NumberOfPeriodsCondsidered], q1[NumberOfPeriodsCondsidered], q2[NumberOfPeriodsCondsidered], sum_u, sum_w, sum_q, W[MaxValueOfStatesInUniform + 1], x, y, z;
		double lamda_1[NumberOfPeriodsCondsidered], lamda_2[NumberOfPeriodsCondsidered], p[NumberOfPeriodsCondsidered];
		long double WT[NumberOfPeriodsCondsidered];
		//static long double pi_store[NumberOfPeriodsCondsidered][nmax + 1][qmax + 1], f_N[NumberOfPeriodsCondsidered][nmax + 1];
		int t, j, n;

		//for (j = 0; j <= qmax; j++) pi[j] = 0;
		//pi[insertOfClass1][insertOfClass2] = 1.0;
		pi[insertOfClass2][insertOfClass1] = 1.0;


		double xunhuan = 99;
		int num_xunhuan = 0;
		while (xunhuan > 1e-4)
		{
			long double temp_pi[MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1];
			for (int i = 0; i < MaxValueOfStatesInUniform; i++)
				for (j = 0; j < MaxValueOfStatesInUniform; j++)
					temp_pi[i][j] = pi[i][j];
			for (t = 0; t < NumberOfPeriodsCondsidered; t++)
			{
				//int temp = t * dt;
				//lamda[t] = lambda[temp];
				lamda_1[t] = LamdaLP[t + 1];// cout << t << " " << lamda_1[t] << endl;
				lamda_2[t] = Lamda[t+1];// cout << t << " " << lamda_2[t] << endl;
				p[t] = solutionOfServers.getSolutionOfServersForTime(t+1);//s[temp];

				w[t] = u[t] = q1[t] = q2[t] = 0;

				for (j = 0; j <= MaxValueOfStatesInUniform; j++) {
					W[j] = 0;
					for (int jj = 0; jj <= MaxValueOfStatesInUniform; jj++) newPi[j][jj] = 0;
				}
				//for (j = 0; j <= qmax; j++) pi_store[t][0][j] = pi[j];

				x = lamda_1[t] + lamda_2[t] + p[t] * MiuForEachPeriod[t+1];
				n = 0;
				y = exp(-1 * x * dt);
				//f_N[t][n] = y;
				do {
					if (y<1e-9 && n>x) break;
					for (int i = 0; i < MaxValueOfStatesInUniform; i++) {
						for (j = 0; j < MaxValueOfStatesInUniform; j++)
						{
							if (i + j < MaxValueOfStatesInUniform) W[i + j] += pi[i][j];
							newPi[i][j] += pi[i][j] * y;



							npi[i][j] = 1 / x * (pi[i][j] * (p[t] * MiuForEachPeriod[t+1] - (MIN(p[t], i + j) * MiuForEachPeriod[t + 1])) +
								pi[i][j + 1] * (MIN(p[t], j + 1) * MiuForEachPeriod[t + 1]) + pi[i + 1][j] * (MIN(MAX(0, p[t] - j), i + 1) * MiuForEachPeriod[t + 1]) +
								((i == 0) ? 0 : pi[i - 1][j] * lamda_1[t]) + ((j == 0) ? 0 : pi[i][j - 1] * lamda_2[t]));

						}
					}
					for (int i = 0; i < MaxValueOfStatesInUniform; i++) {
						for (j = 0; j < MaxValueOfStatesInUniform; j++)
						{
							pi[i][j] = npi[i][j];
							//pi_store[t][n+1][j]= npi[j];
							//if (t == 0 && npi[j] > 0) cout << "pp " << npi[j] << "  " << j << endl;
							w[t] += y * W[j] * MAX(i + j - p[t], 0) * dt / (n + 1);
						}
					}


					//system("pause");
					n++;
					y = y * (x * dt) / n;
					//f_N[t][n] = y;
				} while (n < MaxNumberOfTranstionTimesInUniform);


				z = 0;
				for (int i = 0; i < MaxValueOfStatesInUniform; i++) for (j = 0; j < MaxValueOfStatesInUniform; j++) z += newPi[i][j];

				for (int i = 0; i < MaxValueOfStatesInUniform; i++)
					for (j = 0; j < MaxValueOfStatesInUniform; j++)
					{
						pi[i][j] = newPi[i][j] / z;
						q1[t] += i * pi[i][j];
						q2[t] += j * pi[i][j];

					}
				//u[t] = lamda[t] * dt - q[t] + ((t == 0) ? 0 : q[t - 1]);

				int delta;
				delta = offShift[t+1];
								if (delta > 0)
				{
					for (int i = 0; i < MaxValueOfStatesInUniform; i++)
						for (j = 0; j < MaxValueOfStatesInUniform; j++)
							pi[i][j] = 0;
					for (int i = 0; i <= p[t]; i++)
					{
						for (int jj = 0; i + jj <= p[t]; jj++) {
							int ii = MAX(0, i + jj - delta);
							for (int j = ii; j <= MIN(i + jj, p[t] - delta); j++)
							{
								pi[j - MAX(0, j - i)][MAX(0, j - i)] += HyperGeo_fan(i + jj - j, delta, i + jj, p[t]) * newPi[i][jj] / z;
								//printf("%f\n", HyperGeo_fan(i + jj - j, delta, i + jj, p[t]));
							}
						}

					}
					for (int i = 0; i <= p[t]; i++)
					{
						for (int jj = p[t] - i + 1; jj < MaxValueOfStatesInUniform; jj++) {
							pi[i - MAX(0, delta - jj)][MAX(0, jj - delta)] += newPi[i][jj] / z;
						}
					}
					for (int i = p[t] + 1; i < MaxValueOfStatesInUniform; i++)
					{
						for (int jj = 0; jj < MaxValueOfStatesInUniform; jj++) {
							pi[i - MAX(0, delta - jj)][MAX(0, jj - delta)] += newPi[i][jj] / z;
						}

					}
				}



				//if (t == 20) {
				//	for (j = 0, printf("\n\n"); j <= qmax; j++) printf("%f;", pi[j]);
				//	printf("\n\n");
				//}

			}
			xunhuan = 0;
			for (int i = 0; i < MaxValueOfStatesInUniform; i++)
				for (int j = 0; j < MaxValueOfStatesInUniform; j++)
				{
					xunhuan += abs(temp_pi[i][j] - pi[i][j]);
				}
			num_xunhuan++;
		}


		//printf("\n\n Uniformization \n t; av-u; av-q1; av-w \n");


		//for (t = 0, sum_u = 0, sum_q = 0, sum_w = 0; t < NumberOfPeriodsCondsidered; t++)
		//{
		//	//printf("%d; %f;%f;%f\n", t + 1, u[t], q1[t], w[t]);
		//	//printf("%d; %f\n", t + 1, q[t]);
		//	sum_u += u[t];
		//	sum_w += w[t];
		//	sum_q += q1[t];
		//}
		//printf("\n\n sum; %f;%f;%f\n", sum_u, sum_q, sum_w);
		

		for (t = 0; t < NumberOfPeriodsCondsidered; t++)
		{
			//int temp = t * dt;
			//lamda[t] = lambda[temp];
			lamda_1[t] = LamdaLP[t + 1];
			lamda_2[t] = Lamda[t + 1];
			p[t] = solutionOfServers.getSolutionOfServersForTime(t + 1);//s[temp];

			w[t] = u[t] = q1[t] = q2[t] = 0;

			for (j = 0; j <= MaxValueOfStatesInUniform; j++) {
				W[j] = 0;
				for (int jj = 0; jj <= MaxValueOfStatesInUniform; jj++) newPi[j][jj] = 0;
			}
			//for (j = 0; j <= qmax; j++) pi_store[t][0][j] = pi[j];

			x = lamda_1[t] + lamda_2[t] + p[t] * MiuForEachPeriod[t + 1];
			n = 0;
			y = exp(-1 * x * dt);
			//f_N[t][n] = y;
			do {
				for (int i = 0; i < MaxValueOfStatesInUniform; i++)
					for (j = 0; j < MaxValueOfStatesInUniform; j++)
						FXY_uniform_Piq[t][n][i][j] = pi[i][j];				for (int i = 0; i < MaxValueOfStatesInUniform; i++) {
					for (j = 0; j < MaxValueOfStatesInUniform; j++)
					{
						if (i + j < MaxValueOfStatesInUniform) W[i + j] += pi[i][j];
						newPi[i][j] += pi[i][j] * y;



						npi[i][j] = 1 / x * (pi[i][j] * (p[t] * MiuForEachPeriod[t + 1] - (MIN(p[t], i + j) * MiuForEachPeriod[t + 1])) +
							pi[i][j + 1] * (MIN(p[t], j + 1) * MiuForEachPeriod[t + 1]) + pi[i + 1][j] * (MIN(MAX(0, p[t] - j), i + 1) * MiuForEachPeriod[t + 1]) +
							((i == 0) ? 0 : pi[i - 1][j] * lamda_1[t]) + ((j == 0) ? 0 : pi[i][j - 1] * lamda_2[t]));

					}
				}
				for (int i = 0; i < MaxValueOfStatesInUniform; i++) {
					for (j = 0; j < MaxValueOfStatesInUniform; j++)
					{
						pi[i][j] = npi[i][j];
						//pi_store[t][n+1][j]= npi[j];
						//if (t == 0 && npi[j] > 0) cout << "pp " << npi[j] << "  " << j << endl;
						w[t] += y * W[j] * MAX(i + j - p[t], 0) * dt / (n + 1);
					}
				}

				//system("pause");
				n++;
				y = y * (x * dt) / n;
				//f_N[t][n] = y;
			} while (n < MaxNumberOfTranstionTimesInUniform);


			z = 0;
			for (int i = 0; i < MaxValueOfStatesInUniform; i++) for (j = 0; j < MaxValueOfStatesInUniform; j++) z += newPi[i][j];

			for (int i = 0; i < MaxValueOfStatesInUniform; i++)
				for (j = 0; j < MaxValueOfStatesInUniform; j++)
				{
					pi[i][j] = newPi[i][j] / z;
					q1[t] += i * pi[i][j];
					q2[t] += j * pi[i][j];

				}
			//u[t] = lamda[t] * dt - q[t] + ((t == 0) ? 0 : q[t - 1]);


			int delta;
			delta = offShift[t + 1];
						if (delta > 0)
			{
				for (int i = 0; i < MaxValueOfStatesInUniform; i++)
					for (j = 0; j < MaxValueOfStatesInUniform; j++)
						pi[i][j] = 0;
				for (int i = 0; i <= p[t]; i++)
				{
					for (int jj = 0; i + jj <= p[t]; jj++) {
						int ii = MAX(0, i + jj - delta);
						for (int j = ii; j <= MIN(i + jj, p[t] - delta); j++)
						{
							pi[j - MAX(0, j - i)][MAX(0, j - i)] += HyperGeo_fan(i + jj - j, delta, i + jj, p[t]) * newPi[i][jj] / z;
							//printf("%f\n", HyperGeo_fan(i + jj - j, delta, i + jj, p[t]));
						}
					}

				}
				for (int i = 0; i <= p[t]; i++)
				{
					for (int jj = p[t] - i + 1; jj < MaxValueOfStatesInUniform; jj++) {
						pi[i - MAX(0, delta - jj)][MAX(0, jj - delta)] += newPi[i][jj] / z;
					}
				}
				for (int i = p[t] + 1; i < MaxValueOfStatesInUniform; i++)
				{
					for (int jj = 0; jj < MaxValueOfStatesInUniform; jj++) {
						pi[i - MAX(0, delta - jj)][MAX(0, jj - delta)] += newPi[i][jj] / z;
					}

				}
			}


			//if (t == 20) {
			//	for (j = 0, printf("\n\n"); j <= qmax; j++) printf("%f;", pi[j]);
			//	printf("\n\n");
			//}

		}

		//cout << "test1" << endl;
		//for (int i = 0; i <= 20; i++) {
		//	for (int j = 0; j <= 20; j++)
		//	{
		//		cout <<i<<"  "<<j<<"  " << FXY_uniform_Piq[0][1][i][j] << "\n";
		//	}
		//	cout << endl;
		//}
		//FXY_uniform_Piq[NumberOfPeriodsCondsidered][MaxNumberOfTranstionTimesInUniform + 1][MaxValueOfStatesInUniform + 1][MaxValueOfStatesInUniform + 1];

		for (t = 0, sum_u = 0, sum_q = 0, sum_w = 0; t < NumberOfPeriodsCondsidered; t++)
		{
			printf("%d; %f;%f;%f\n", t + 1, u[t], q1[t], w[t]);
			//printf("%d; %f\n", t + 1, q[t]);
			sum_u += u[t];
			sum_w += w[t];
			sum_q += q1[t];
		}
		printf("\n\n sum; %f;%f;%f\n", sum_u, sum_q, sum_w);
		cout << "uniformization times:" << num_xunhuan << endl;
	}

	double timetemp[2];
	double HyperGeo_4d(int q1, int q2, int q3, int idle, int n1, int n2, int n3, int n)	{
		timetemp[0] = clock();
		int q = q1 + q2 + q3 + idle;
		double res = 1;
		double da = q, xiao = n;
		while (xiao >= 1) {
			res *= (xiao / da);
			--xiao; --da;
		}
		da = q1; xiao = n1;
		while (xiao >= 1) {
			res *= (da / xiao);
			--xiao; --da;
		}
		da = q2; xiao = n2;
		while (xiao >= 1) {
			res *= (da / xiao);
			--xiao; --da;
		}
		da = q3; xiao = n3;
		while (xiao >= 1) {
			res *= (da / xiao);
			--xiao; --da;
		}
		da = idle; xiao = n - n1 - n2 - n3;
		while (xiao >= 1) {
			res *= (da / xiao);
			--xiao; --da;
		}
		timetemp[1] = clock();
		totaltimeOfHyper4d += double(timetemp[1] - timetemp[0]) / double(CLOCKS_PER_SEC);
		return res;
	}


	void Uniformization_priority_FXY_forthreetype(const SolutionOnShift solutionOnShift, const SolutionOfServers solutionOfServers, const Shifts shiftEnv)
	{
		FXY_uniform_Piq_Mid = vector<vector<vector<double>>>(NumberOfPeriodsCondsidered, vector<vector<double>>(MaxNumberOfTranstionTimesInUniform + 1,
			vector<double>(QMAXSUMS + 1, 0)));
		FXY_uniform_Piq_Low = vector<vector<vector<double>>>(NumberOfPeriodsCondsidered, vector<vector<double>>(MaxNumberOfTranstionTimesInUniform + 1,
			vector<double>(QMAXSUM + 1, 0)));


		long double pi[QMAX1 + 1][QMAX2 + 1][QMAX3 + 1][SMAX + 1]{}, npi[QMAX1 + 1][QMAX2 + 1][QMAX3 + 1][SMAX + 1]{}, newPi[QMAX1 + 1][QMAX2 + 1][QMAX3 + 1][SMAX + 1]{};		long double w[NumberOfPeriodsCondsidered], u[NumberOfPeriodsCondsidered], q1[NumberOfPeriodsCondsidered], q2[NumberOfPeriodsCondsidered], q3[NumberOfPeriodsCondsidered], sum_u, sum_w, sum_q1, sum_q2, sum_q3, W[QMAXSUM + 1], x, y, z;
		double lamda_1[NumberOfPeriodsCondsidered], lamda_2[NumberOfPeriodsCondsidered], lamda_3[NumberOfPeriodsCondsidered], p[NumberOfPeriodsCondsidered];
		long double WT[NumberOfPeriodsCondsidered];
		//static long double pi_store[NumberOfPeriodsCondsidered][nmax + 1][QMAX + 1], f_N[NumberOfPeriodsCondsidered][nmax + 1];
		int t, j, n;

		//for (j = 0; j <= QMAX; j++) pi[j] = 0;
		pi[insertOfClass1][insertOfClass2][insertOfClass3][min(insertOfClass3,int(solutionOfServers.getSolutionOfServersForTime(1) +1e-6))] = 1.0;
		for (t = 0; t < NumberOfPeriodsCondsidered; t++) p[t] = solutionOfServers.getSolutionOfServersForTime(t+1);
		double xunhuan = 99;		int num_xunhuan = 0;

		while (xunhuan > 1)
		//while (xunhuan > 98)
		{
			FXY_uniform_Piq_Mid = vector<vector<vector<double>>>(NumberOfPeriodsCondsidered, vector<vector<double>>(MaxNumberOfTranstionTimesInUniform + 1,
				vector<double>(QMAXSUMS + 1, 0)));
			FXY_uniform_Piq_Low = vector<vector<vector<double>>>(NumberOfPeriodsCondsidered, vector<vector<double>>(MaxNumberOfTranstionTimesInUniform + 1,
				vector<double>(QMAXSUM + 1, 0)));
			long double temp_pi[QMAX1 + 1][QMAX2 + 1][QMAX3 + 1][SMAX + 1];
			for (int i = 0; i < QMAX1; i++)
				for (j = 0; j < QMAX2; j++)
					for (int k = 0; k < QMAX3; k++)
						for (int s3 = 0; s3 <= SMAX; ++s3)
							temp_pi[i][j][k][s3] = pi[i][j][k][s3];
			for (t = 0; t < NumberOfPeriodsCondsidered; t++)
			{
				cout << "t " << t << endl;
				int temp = t * dt;
				//lamda[t] = lambda[temp];
				lamda_2[t] = LamdaMP[t + 1];//lambda_1[t];
				lamda_1[t] = Lamda[t + 1];//lambda_1[t]+lambda_2[t];
				lamda_3[t] = LamdaLP[t + 1];//lambda_1[t]+lambda_2[t];
				p[t] = solutionOfServers.getSolutionOfServersForTime(t + 1);

				w[t] = u[t] = q1[t] = q2[t] = q3[t] = 0;

				for (j = 0; j <= QMAXSUM; j++) {
					W[j] = 0;
					//for (int jj = 0; jj <= QMAX; jj++) newPi[j][jj] = 0;
				}
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= SMAX; ++s3)
								newPi[i][j][k][s3] = 0;
				//for (j = 0; j <= QMAX; j++) pi_store[t][0][j] = pi[j];

				x = lamda_1[t] + lamda_2[t] + lamda_3[t] + p[t] * MiuForEachPeriod[t + 1];
				n = 0;
				y = exp(-1 * x * dt);
				//f_N[t][n] = y;
				do {
					if (y<1e-4 && n>x) break;
					//cout << n << endl;
//#pragma omp parallel
//#pragma omp for
					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) npi[i][j][k][s3] = 0;

					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) {
									FXY_uniform_Piq_Mid[t][n][i + j + s3] += pi[i][j][k][s3];
									FXY_uniform_Piq_Low[t][n][i + j + k] += pi[i][j][k][s3];
								}

					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(p[t], SMAX); ++s3)
								{
									//if (i + j < QMAX) W[i + j] += pi[i][j];
									/*if (i == 0 && j == 0 && k == 0 && s3 == 0)
										cout << "d" << endl;*/
									newPi[i][j][k][s3] += pi[i][j][k][s3] * y;
									double check_com = 0;
																		if (i + 1 + j + k > p[t])									{
										int s_23 = MAX(p[t] - i, 0);
										if (s_23 > 0) {
											npi[i + 1][j][k][s3] += pi[i][j][k][s3] * lamda_1[t] * double(MAX(MIN(p[t] - i - s3, j), 0)) / double(s_23);											check_com += pi[i][j][k][s3] * lamda_1[t] * double(MAX(MIN(p[t] - i - s3, j), 0)) / double(s_23);
																						if (s3 > 0) {
												npi[i + 1][j + 1][k - 1][s3 - 1] += pi[i][j][k][s3] * lamda_1[t] * double(s3) / double(s_23);												check_com += pi[i][j][k][s3] * lamda_1[t] * double(s3) / double(s_23);
											}										}
										else {
											npi[i + 1][j][k][s3] += pi[i][j][k][s3] * lamda_1[t];
											check_com += pi[i][j][k][s3] * lamda_1[t];
										}
									}
									else {
										npi[i + 1][j][k][s3] += pi[i][j][k][s3] * lamda_1[t];
										check_com += pi[i][j][k][s3] * lamda_1[t];
									}
																		npi[i][j + 1][k][s3] += pi[i][j][k][s3] * lamda_2[t];
									check_com += pi[i][j][k][s3] * lamda_2[t];
																		if (i + j + k < p[t]) {
										npi[i][j][k + 1][s3 + 1] += pi[i][j][k][s3] * lamda_3[t];
										check_com += pi[i][j][k][s3] * lamda_3[t];
									}
									else {
										npi[i][j][k + 1][s3] += pi[i][j][k][s3] * lamda_3[t];
										check_com += pi[i][j][k][s3] * lamda_3[t];
									}

																		if (i > 0) {
										if (i - 1 + j + s3 < p[t]) {											if (k > s3) {
												npi[i - 1][j][k][s3 + 1] += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
											}
											else {
												npi[i - 1][j][k][s3] += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
											}
										}
										else {
											npi[i - 1][j][k][s3] += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
											check_com += pi[i][j][k][s3] * MIN(p[t], i) * MiuForEachPeriod[t + 1];
										}
									}

																		if (j > 0) {
										if (i + j - 1 + s3 < p[t]) {											if (k > s3) {
												npi[i][j - 1][k][s3 + 1] += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
											}
											else {
												npi[i][j - 1][k][s3] += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
											}
										}
										else {
											npi[i][j - 1][k][s3] += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
											check_com += pi[i][j][k][s3] * double(MAX(0, MIN(j, p[t] - i - s3))) * MiuForEachPeriod[t + 1];
										}
									}

																		if (s3 > 0) {										if (i + j + s3 - 1 < p[t])
										{
											if (k > s3) {
												npi[i][j][k - 1][s3] += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
											}
											else {
												npi[i][j][k - 1][s3 - 1] += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
												check_com += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
											}
										}
										else {
											npi[i][j][k - 1][s3 - 1] += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
											check_com += pi[i][j][k][s3] * s3 * MiuForEachPeriod[t + 1];
										}
									}



																		npi[i][j][k][s3] += MAX(p[t] - i - j - k, 0) * MiuForEachPeriod[t + 1] * pi[i][j][k][s3];
									check_com += MAX(p[t] - i - j - k, 0) * MiuForEachPeriod[t + 1] * pi[i][j][k][s3];
									check_com /= x;
									//if (abs(check_com - pi[i][j][k][s3]) > 1e-5) cout << i << " " << j << " " << k << " " << s3 << " " << check_com << " " << pi[i][j][k][s3] << endl;

									/*npi[i][j] = 1 / x * (pi[i][j] * (p[t] * mu[t] - (MIN(p[t], i + j)*mu[t])) +
										pi[i][j + 1] * (MIN(p[t], j + 1)*mu[t]) + pi[i + 1][j] * (MIN(MAX(0, p[t] - j), i + 1)*mu[t]) +
										((i == 0) ? 0 : pi[i - 1][j] * lamda_1[t]) + ((j == 0) ? 0 : pi[i][j - 1] * lamda_2[t]));*/
										//npi[i][j][k][s3] /= x;
								}

					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3)
								{
									pi[i][j][k][s3] = npi[i][j][k][s3] / x;
									//pi_store[t][n+1][j]= npi[j];
									//if (t == 0 && npi[j] > 0) cout << "pp " << npi[j] << "  " << j << endl;
									//w[t] += y * W[j] * MAX(i + j - p[t], 0)* dt / (n + 1);
								}


					//system("pause");
					double dddd = 0;
					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) {
									dddd += pi[i][j][k][s3];
									//if (pi[i][j][k][s3] > 0) cout << i << " "
										//<< j << " " << k << " " << s3 << " " << pi[i][j][k][s3] << endl;
								}
										n++;
					y = y * (x * dt) / n;
					//f_N[t][n] = y;
				} while (n < MaxNumberOfTranstionTimesInUniform);
				double dddd = 0;
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) dddd += pi[i][j][k][s3];
				if (abs(dddd - 1) > 1e-5) std::cout << "prob sum:" << dddd << endl;

				z = 0;
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) z += newPi[i][j][k][s3];

				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3)
							{
								pi[i][j][k][s3] = newPi[i][j][k][s3] / z;
								q1[t] += i * pi[i][j][k][s3];
								q2[t] += j * pi[i][j][k][s3];
								q3[t] += k * pi[i][j][k][s3];

							}

				dddd = 0;
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) dddd += pi[i][j][k][s3];
				if (abs(dddd - 1) > 1e-5) std::cout << "prob sum 222:" << dddd << endl;
				//u[t] = lamda[t] * dt - q[t] + ((t == 0) ? 0 : q[t - 1]);

				int delta;
				delta = offShift[t + 1];
				int t_next = ((t == NumberOfPeriodsCondsidered - 1) ? 0 : t + 1), p_delta = p[t] - delta;
				int n_on = int(p[t_next]) - p_delta; //cout << p[t_next] << " " << p_delta <<" "<<t_next << endl;
								if (delta > 0)
				{
					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3)
									pi[i][j][k][s3] = 0;

					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) {
									int tmp_q1 = MIN(i, p[t]), tmp_q2 = MAX(MIN(p[t] - i - s3, j), 0), tmp_q3 = s3, tmp_idle = MAX(p[t] - i - j - k, 0);
									double check = 0;
									for (int tmp_n1 = 0; tmp_n1 <= MIN(delta, tmp_q1); ++tmp_n1)
										for (int tmp_n2 = 0; tmp_n2 <= MIN(delta - tmp_n1, tmp_q2); ++tmp_n2) {
											for (int tmp_n3 = 0; tmp_n3 <= MIN(delta - tmp_n1 - tmp_n2, tmp_q3); ++tmp_n3) {
												int n_idle = delta - tmp_n1 - tmp_n2 - tmp_n3;
												if (n_idle > tmp_idle) continue;
												int no_busy = p_delta - (i - tmp_n1 + j - tmp_n2 + s3 - tmp_n3);
												if (no_busy > 0) pi[i - tmp_n1][j - tmp_n2][k - tmp_n3][s3 - tmp_n3 + MIN(no_busy, k - s3)] +=
													HyperGeo_4d(tmp_q1, tmp_q2, tmp_q3, tmp_idle, tmp_n1, tmp_n2, tmp_n3, delta) * newPi[i][j][k][s3] / z;
												else pi[i - tmp_n1][j - tmp_n2][k - tmp_n3][s3 - tmp_n3] +=
													HyperGeo_4d(tmp_q1, tmp_q2, tmp_q3, tmp_idle, tmp_n1, tmp_n2, tmp_n3, delta) * newPi[i][j][k][s3] / z;
												check += HyperGeo_4d(tmp_q1, tmp_q2, tmp_q3, tmp_idle, tmp_n1, tmp_n2, tmp_n3, delta);

											}
										}
									/*if (abs(check - 1) > 1e-5) { cout << "check " << check << endl;
									cout << i << " " << j << " " << k << " " << s3 << endl;
									}*/
								}

				}

				dddd = 0;
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= SMAX; ++s3) dddd += pi[i][j][k][s3];
								//cout <<"qian " << pi[0][0][1][0] << endl;
				//cout << n_on << endl;
				if (n_on > 0) {
					for (int i = 0; i < QMAX1; i++)
						for (j = 0; j < QMAX2; j++)
							for (int k = 0; k < QMAX3; k++)
								for (int s3 = 0; s3 <= MIN(SMAX, p[t]); ++s3) {
									//if (i == 0 && j == 0 && k == 1 && s3 == 0) cout << "d " << endl;
									if (k > s3 && i + j + s3 < p[t_next]) {
										if (pi[i][j][k][s3] <= 0) continue;
										if (i == 0 && j == 0 && k == 1 && MIN(p_delta + n_on - i - j, k) == 0) {
											cout << i << " " << j << " " << k << " " << s3 << " " << MIN(p_delta + n_on - i - j, k) << endl; system("pause");
										}
										//if (pi[i][j][k][MIN(p_delta + n_on - i - j, k)] > 0) cout << i << " " << j << " " << k << " " << s3 << " " << MIN(p_delta + n_on - i - j, k) <<" "<< pi[i][j][k][MIN(p_delta + n_on - i - j, k)]<< endl;
										pi[i][j][k][MIN(p_delta + n_on - i - j, k)] += pi[i][j][k][s3];
										pi[i][j][k][s3] = 0;
									}
								}
				}

				//cout << "hou " << pi[0][0][1][0] << endl;
				dddd = 0;
				for (int i = 0; i < QMAX1; i++)
					for (j = 0; j < QMAX2; j++)
						for (int k = 0; k < QMAX3; k++)
							for (int s3 = 0; s3 <= SMAX; ++s3) dddd += pi[i][j][k][s3];
				if (abs(dddd - 1) > 1e-5) std::cout << "prob sum on: " << dddd << endl;


				//if (t == 20) {
				//	for (j = 0, printf("\n\n"); j <= QMAX; j++) printf("%f;", pi[j]);
				//	printf("\n\n");
				//}

			}
			xunhuan = 0;
			for (int i = 0; i < QMAX1; i++)
				for (j = 0; j < QMAX2; j++)
					for (int k = 0; k < QMAX3; k++)
						for (int s3 = 0; s3 <= SMAX; ++s3)
						{
							xunhuan += abs(temp_pi[i][j][k][s3] - pi[i][j][k][s3]);
						}
			num_xunhuan++;
			cout << num_xunhuan << "\t" << xunhuan << endl;
		}
		


		printf("\n\n Uniformization \n t; av-u; av-q1;av-q2;av-q3; av-w \n");


		for (t = 0, sum_u = 0, sum_q1 = 0, sum_q2 = 0, sum_q3 = 0, sum_w = 0; t < NumberOfPeriodsCondsidered; t++)
		{
			printf("%d; %f;%f;%f;%f;%f\n", t + 1, u[t], q1[t], q2[t], q3[t], w[t]);
			//printf("%d; %f\n", t + 1, q[t]);
			sum_u += u[t];
			sum_w += w[t];
			sum_q1 += q1[t];
			sum_q2 += q2[t];
			sum_q3 += q3[t];
		}
		printf("\n\n sum; %f;%f;%f;%f;%f\n", sum_u, sum_q1, sum_q2, sum_q3, sum_w);
		
		
		cout << "test2" << endl;
		for (int i = 0; i <= 20; i++) {
			//for (int j = 0; j <= 20; j++)
			{
				cout << i << "  " << j << "  " << FXY_uniform_Piq_Mid[0][1][i] << "\n";
			}
			//cout << endl;
		}
		
		cout << "test3" << endl;
		for (int i = 0; i <= 20; i++) {
			//for (int j = 0; j <= 20; j++)
			{
				cout << i << "  " << j << "  " << FXY_uniform_Piq_Low[0][1][i] << "\n";
			}
			//cout << endl;
		}
		
		
		cout << num_xunhuan << endl;

	}

	//***************************************************end**********************************************

};




#endif // !_computeServiceLevel_Uniformization_
#include "approximate.h"
//#include "computeServiceLevel_Uniformization.h"
#define integralAccuracy 200
#define taoStarAccuracy 10
#define integralAccuracyForT 200
//#include "computeServiceLevel_Uniformization.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>


unordered_map<string, double> valueOfA_MAP;
double APPROXIMATE::valueOfA(double t, double tao) {
	//find map
	string tempKey = to_string(t) + "inter" + to_string(tao);
	if (valueOfA_MAP.count(tempKey)) return valueOfA_MAP[tempKey];
	
	double startTime = t;
	SolutionOfServers solutionOfServers = solveForApproximate.currentSolutionOfServers;

	double endTime = startTime + tao;
	int firstPeriod = int(startTime) + 1;
	int lastPeriod = int(endTime) + 1;
	double sum = 0;
	for (int indexOfPeriod = firstPeriod; indexOfPeriod <= lastPeriod; indexOfPeriod++) {
		int tempIndexOfPeriod = indexOfPeriod;
		while (tempIndexOfPeriod > NumberOfPeriodsCondsidered)tempIndexOfPeriod = tempIndexOfPeriod - NumberOfPeriodsCondsidered;
		sum += (min(endTime, double(indexOfPeriod)) - max(startTime, double(indexOfPeriod) - 1.0))
			* double(MiuForEachPeriod[tempIndexOfPeriod]) * double(solutionOfServers.getSolutionOfServersForTime(tempIndexOfPeriod));
	}

	//save map
	valueOfA_MAP[tempKey] = sum;
	return sum;
}

unordered_map<string, double> valueOfAarr_MAP;
double APPROXIMATE::valueOfAarr(double t, double tao) {
	//find map
	string tempKey = to_string(t) + "inter" + to_string(tao);
	if (valueOfAarr_MAP.count(tempKey)) return valueOfAarr_MAP[tempKey];

	double startTime = t;
	SolutionOfServers solutionOfServers = solveForApproximate.currentSolutionOfServers;

	double endTime = startTime + tao;
	int firstPeriod = int(startTime) + 1;
	int lastPeriod = int(endTime) + 1;
	double sum = 0;
	for (int indexOfPeriod = firstPeriod; indexOfPeriod <= lastPeriod; indexOfPeriod++) {
		int tempIndexOfPeriod = indexOfPeriod;
		while (tempIndexOfPeriod > NumberOfPeriodsCondsidered)tempIndexOfPeriod = tempIndexOfPeriod - NumberOfPeriodsCondsidered;
		sum += (min(endTime, double(indexOfPeriod)) - max(startTime, double(indexOfPeriod) - 1.0))
			*double(Lamda[tempIndexOfPeriod]);	}
	//save map
	valueOfAarr_MAP[tempKey] = sum;
	return sum;
}

unordered_map<string, double> integralPai_MAP;
unordered_map<string, double> integralPaiMultA_MAP;
unordered_map<string, double> integralPaiMultAarr_MAP;

double APPROXIMATE::integralPai(int z, double LB, double UB) {
	//find map
	string tempKey = to_string(z) + "inter" + to_string(LB)+ "inter1" + to_string(UB);
	if (integralPai_MAP.count(tempKey)) return integralPai_MAP[tempKey];


	double sum = 0;
	for (int index = 0; index < integralAccuracy; index++) {
		double t = (UB - LB) / double(integralAccuracy) * double(index) + LB;
		//if (t > 24)cout << "error jhala\n";
		sum += valueOfPai_Z_T(z, t);
	}
	sum = sum / double(integralAccuracy);

	//save map
	integralPai_MAP[tempKey] = sum;
	return sum;
}
double APPROXIMATE::integralPaiMultA(int z, double LB, double UB, double tao) {
	//find map
	string tempKey = to_string(z) + "inter" + to_string(LB) + "inter1" + to_string(UB) + "inter2" + to_string(tao);
	if (integralPaiMultA_MAP.count(tempKey)) return integralPaiMultA_MAP[tempKey];


	//if (UB > 24)cout << "error hadfkhkah" << endl;
	//cout << "integralpaimultA\t" << LB << "\t" << UB << "\t" << tao << endl;


	double sum = 0;
	for (int index = 0; index < integralAccuracy; index++) {
		double t = (UB - LB) / double(integralAccuracy) * double(index) + LB;
		sum += valueOfPai_Z_T(z, t) * valueOfA(t, tao);
	}
	sum = sum / double(integralAccuracy);

	//save map
	integralPaiMultA_MAP[tempKey] = sum;
	return sum;
}
double APPROXIMATE::integralPaiMultA_arr(int z, double LB, double UB, double tao) {
	//find map
	string tempKey = to_string(z) + "inter" + to_string(LB) + "inter1" + to_string(UB) + "inter2" + to_string(tao);
	if (integralPaiMultAarr_MAP.count(tempKey)) return integralPaiMultAarr_MAP[tempKey];


	double sum = 0;
	for (int index = 0; index < integralAccuracy; index++) {
		double t = (UB - LB) / double(integralAccuracy) * double(index) + LB;
		sum += valueOfPai_Z_T(z, t) * valueOfAarr(t, tao);
	}
	sum = sum / double(integralAccuracy);

	//save map

	integralPaiMultAarr_MAP[tempKey] = sum;
	return sum;
}


double APPROXIMATE::expectedServiceRate_Abar(int l, int z, double LB, double UB, double tao) {
	return integralPaiMultA(z, LB, UB, tao) / integralPai(z, LB, UB);
}
double APPROXIMATE::expectedArrivalRate_Abar_arr(int l, int z, double LB, double UB, double tao) {
	return integralPaiMultA_arr(z, LB, UB, tao) / integralPai(z, LB, UB);
}


#define accurancyOfSkellam_arrival 1000
#define accurancyOfSkellam_service 1000
const double maxValueOfArrivalRate = 30;
const double maxValueOfServiceRate = 30;

int LB_mu1_mu2[accurancyOfSkellam_service + 1][accurancyOfSkellam_arrival + 1] = { 0 };
int UB_mu1_mu2[accurancyOfSkellam_service + 1][accurancyOfSkellam_arrival + 1] = { 0 };
vector<vector<vector<double>>> valueOfSkellam;

double APPROXIMATE::readfileOfSkellam() {
	valueOfSkellam = vector<vector<vector<double>>>(accurancyOfSkellam_service + 1, vector<vector<double>>(accurancyOfSkellam_arrival + 1,
		vector<double>()));


	//ifstream inFile("C://Users//15366//try_skellam.csv", ios::in);
	//ifstream inFile("D://wangck//try_skellam1000.csv", ios::in);
	ifstream inFile("C://Users//Administrator//Desktop//skellam-1000.csv", ios::in);

	//ifstream inFile("D://wangck//try_skellam100.csv", ios::in);

	string lineStr;
	vector<vector<string>> strArray;
	int i, j;
	i = 0;
	char* end;
	if (inFile.fail())
		cout << "file read fail" << endl;
	while (getline(inFile, lineStr))
	{
		j = 0;
				//cout << lineStr << endl;
				stringstream ss(lineStr);
		string str;
		vector<string> lineArray;
				getline(ss, str, ',');
		double mu1 = static_cast<double>(atof(str.c_str()));
		int indexOfMu1 = int(0.001 + mu1 / double(maxValueOfServiceRate) * double(accurancyOfSkellam_service));

		getline(ss, str, ',');
		double mu2 = static_cast<double>(atof(str.c_str()));
		int indexOfMu2 = int(0.001 + mu2 / double(maxValueOfArrivalRate) * double(accurancyOfSkellam_arrival));

		getline(ss, str, ',');
		LB_mu1_mu2[indexOfMu1][indexOfMu2] = static_cast<int>(atoi(str.c_str()));
		getline(ss, str, ',');
		UB_mu1_mu2[indexOfMu1][indexOfMu2] = static_cast<int>(atoi(str.c_str()));
		while (getline(ss, str, ','))
		{
			valueOfSkellam[indexOfMu1][indexOfMu2].push_back(static_cast<double>(atof(str.c_str())));
		}
	}

}


double APPROXIMATE::skellam(double serviceRate, double arrivalRate, int value){
	if (arrivalRate > maxValueOfArrivalRate || serviceRate > maxValueOfServiceRate) { cout << "error ahdsjkhakljh:skellam overheap" << endl; system("pause"); return -1; }

	int indexOfArrivalRate = int(arrivalRate / double(maxValueOfArrivalRate) * double(accurancyOfSkellam_arrival));
	int indexOfServiceRate = int(serviceRate / double(maxValueOfServiceRate) * double(accurancyOfSkellam_service));
	int LB = LB_mu1_mu2[indexOfServiceRate][indexOfArrivalRate];
	int UB = UB_mu1_mu2[indexOfServiceRate][indexOfArrivalRate];

	//cout << "###:\t" << indexOfArrivalRate << "\t" << indexOfServiceRate << "\t" << arrivalRate << "\t" << serviceRate << "\t" << UB << "\t" << LB << "\t" << value << endl;


	if (indexOfArrivalRate == 0 || indexOfServiceRate == 0)return this->skellam_dealwith0(indexOfArrivalRate, indexOfServiceRate, arrivalRate, serviceRate, value);
	
	if (value < LB)return 0;
	if (value >= UB)return 1;

	/*cout << indexOfArrivalRate << "  ****   " << indexOfServiceRate << endl;
	cout << "LB\t" << LB << endl;
	cout << "UB\t" << UB << endl;*/

	return valueOfSkellam[indexOfServiceRate][indexOfArrivalRate][value - LB];
}



double APPROXIMATE::valueOfBeta2_Left(int l, int z, double tao) {
	double theta = tao - double(int(tao));//todo
	//cout << "left1\t" << l << "\t" << z << "\ttheta\t" << theta << "\t" << tao << endl;

	double Abar = expectedServiceRate_Abar(l, z, double(l) - 1.0, double(l) - theta, tao);
	//cout << Abar << endl;


	double Abar_arr = expectedArrivalRate_Abar_arr(l, z, double(l) - 1.0, double(l) - theta, tao);
	//cout << Abar_arr << endl;

	int startPeriod = l;
	int endPeriod = int(l - 1 + tao - (1e-9)) + 1;	int sumOffDuty = 0;
	for (int l1 = startPeriod; l1 < endPeriod; l1++) {		if (l1 > NumberOfPeriodsCondsidered) { cout << "error adjakhahdjhkxxx\n"; system("pause"); }
		sumOffDuty += this->offShift[l1];	}
	if (endPeriod > NumberOfPeriodsCondsidered)endPeriod = endPeriod - NumberOfPeriodsCondsidered;
	int valueOfCustomersInQueue = z - sumOffDuty - solveForApproximate.currentSolutionOfServers.getSolutionOfServersForTime(endPeriod);
	double tempSkellam = skellam(Abar, Abar_arr, valueOfCustomersInQueue);
	//cout << "left\t" << Abar << "\t" << Abar_arr << "\t" << valueOfCustomersInQueue << "\t" << tempSkellam << endl;

	return tempSkellam;
}

double APPROXIMATE::valueOfBeta2_Right(int l, int z, double tao) {	double theta = tao - double(int(tao));//todo

	//cout <<"right1\t"<< l << "\t" << z << "\ttheta\t" << theta << "\t" << tao << endl;
	double Abar = expectedServiceRate_Abar(l, z, double(l) - theta, double(l), tao);
	//cout << Abar << endl;
	
	double Abar_arr = expectedArrivalRate_Abar_arr(l, z, double(l) - theta, double(l), tao);
	//cout << Abar_arr << endl;

	int startPeriod = l;
	int endPeriod = int(l + tao - (1e-5)) + 1;	int sumOffDuty = 0;
	for (int l1 = startPeriod; l1 < endPeriod; l1++) {		if (l1 > NumberOfPeriodsCondsidered) { cout << "error adjakhahdjhk\n"; system("pause"); }
		int templ1 = l1; if (templ1 > NumberOfPeriodsCondsidered)templ1 = templ1 - NumberOfPeriodsCondsidered;
		sumOffDuty += this->offShift[templ1];	}
	if (endPeriod > NumberOfPeriodsCondsidered)endPeriod = endPeriod - NumberOfPeriodsCondsidered;
	int valueOfCustomersInQueue = z - sumOffDuty - solveForApproximate.currentSolutionOfServers.getSolutionOfServersForTime(endPeriod);
	double tempSkellam = skellam(Abar, Abar_arr, valueOfCustomersInQueue);
	//cout << "right\t" << Abar << "\t" << Abar_arr << "\t" << valueOfCustomersInQueue << "\t" << tempSkellam << endl;

	return tempSkellam;
}

double APPROXIMATE::valueOfBeta2_new(int l, int z, double tao, double left, double right) {
	//double theta = tao - double(int(tao));//todo

	//cout <<"right1\t"<< l << "\t" << z << "\ttheta\t" << theta << "\t" << tao << endl;
	double Abar = expectedServiceRate_Abar(l, z, left, right, tao);
	//cout << Abar << endl;

	double Abar_arr = expectedArrivalRate_Abar_arr(l, z, left, right , tao);
	//cout << Abar_arr << endl;

	int startPeriod = left;
	int endPeriod = int(right - (1e-5)) + 1;	int sumOffDuty = 0;
	for (int l1 = startPeriod; l1 < endPeriod; l1++) {		if (l1 > NumberOfPeriodsCondsidered) { cout << "error adjakhahdjhk\n"; system("pause"); }
		int templ1 = l1; if (templ1 > NumberOfPeriodsCondsidered)templ1 = templ1 - NumberOfPeriodsCondsidered;
		sumOffDuty += this->offShift[templ1];	}
	if (endPeriod > NumberOfPeriodsCondsidered)endPeriod = endPeriod - NumberOfPeriodsCondsidered;
	int valueOfCustomersInQueue = z - sumOffDuty - solveForApproximate.currentSolutionOfServers.getSolutionOfServersForTime(endPeriod);
	double tempSkellam = skellam(Abar, Abar_arr, valueOfCustomersInQueue);
	//cout << "right\t" << Abar << "\t" << Abar_arr << "\t" << valueOfCustomersInQueue << "\t" << tempSkellam << endl;

	return tempSkellam;
}


double APPROXIMATE::valueOfBeta2_Choose(int l, int z, double tao, double t) {
	if (int(t + tao) == int(l - 1 + tao)) {
		return valueOfBeta2_Left(l,z,tao);	}
	else if (int(t + tao) == int(l + tao)) {
		return valueOfBeta2_Right(l, z, tao);
	}
	else { cout << "error jaklhfdakhj" << endl; system("pause"); }
}

//double APPROXIMATE::valueOfBeta2_Choose_new(int l, int z, double tao, double t, ) {
//	if (int(t + tao) == int(l - 1 + tao)) {
//	}
//	else if (int(t + tao) == int(l + tao)) {
//		return valueOfBeta2_Right(l, z, tao);
//	}
//	else { cout << "error jaklhfdakhj" << endl; system("pause"); }
//}


double APPROXIMATE::valueOfBeta1_Choose(int l, int z, double maxOfTao, double t) {
	double minValue = 10000;
	for (int indexOfTaoStar = 1; indexOfTaoStar <= taoStarAccuracy; indexOfTaoStar++) {
		double valueOfTaoStar = maxOfTao / double(taoStarAccuracy) * double(indexOfTaoStar);
		double tempValue = valueOfBeta2_Choose(l, z, valueOfTaoStar, t);
		//cout << "tao\t" << valueOfTaoStar << "\tvalue\t" << tempValue << endl;
		minValue = min(minValue, tempValue);
	}
	return minValue;
}



//double APPROXIMATE::valueOfBeta1_Left(int l, int z, double maxOfTao) {
//	double minValue = 10000;
//	for (int indexOfTaoStar = 1; indexOfTaoStar <= taoStarAccuracy; indexOfTaoStar++) {
//		double valueOfTaoStar = maxOfTao / double(taoStarAccuracy) * double(indexOfTaoStar);
//		minValue = min(minValue, valueOfBeta2_Left(l, z, valueOfTaoStar));
//	}
//	return minValue;
//}
//double APPROXIMATE::valueOfBeta1_Right(int l, int z, double maxOfTao) {
//	double minValue = 10000;
//	for (int indexOfTaoStar = 1; indexOfTaoStar <= taoStarAccuracy; indexOfTaoStar++) {
//		double valueOfTaoStar = maxOfTao / double(taoStarAccuracy) * double(indexOfTaoStar);
//		minValue = min(minValue, valueOfBeta2_Right(l, z, valueOfTaoStar));
//	}
//	return minValue;
//}




double APPROXIMATE::foreachZ(int l, int z, double maxOfTao) {
	double sum = 0;
	for (int indexOfT = 0; indexOfT < integralAccuracyForT; indexOfT++) {
		double valueOfT = double(l - 1) + double(indexOfT) / double(integralAccuracyForT);



		sum += valueOfBeta1_Choose(l, z, maxOfTao, valueOfT) * valueOfPai_Z_T(z, valueOfT);
	}
	sum = sum / double(integralAccuracyForT);

		
	//sum = sum / integralPai(z, double(l - 1), double(l));

	return sum;
}

void APPROXIMATE::clearMap() {
	valueOfA_MAP.clear();
	valueOfAarr_MAP.clear();
	integralPai_MAP.clear();
	integralPaiMultA_MAP.clear();
	integralPaiMultAarr_MAP.clear();
}
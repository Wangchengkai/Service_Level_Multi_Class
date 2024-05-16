#pragma once
#ifndef _SIMU_MHZ_
#define _SIMU_MHZ_

#include"data.h"
#include"operFunctions.h"
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include<limits>
#include<string>
#include<fstream>
#include <sstream>
#include<algorithm>
#include<list>
#include<vector>
#include<ctime>


extern int num_shift;
const int T = NumberOfPeriodsCondsidered;
const double deta = LengthOfOnePeriod;
const long num_sim = 1e6;//1e6
const long num_warm = 1e4;//1e4
#define EXP(mean) (-log(1.0-  ((double) rand() / (RAND_MAX))) *mean)
using namespace std;

extern double lamda[];
extern double lamda1[];
extern double lamda2[];
extern double miu[];


//const double miu[T] = { 5.387096774,2.147465438,2.529953917,2.004032258,2.189516129,2.588709677,2.012096774,2.045698925,1.625310174,1.327956989,2.946236559,2.470046083,1.951612903,2.161290323,2.784946237,2.204301075,3.774193548,  3.677419355,2.516129032,2.822580645,2.596774194,2.419354839,2.096774194,0.983870968, };
extern double alpha[T], alpha1[T];
//extern int sch[TotalNumberOfServers][T];





struct Patient
{
	int period;
	double arrival_time;
	int days;	Patient() {};
	Patient(int p, double at)
	{
		this->period = p;
		this->arrival_time = at;
		days = 0;
	}
};

class Server
{
public:
	Server();
	Server(int sche[T], const double miu_in[T]);
	Server(int sche[T], const double miu_in[T], int idx);
	Server(vector<int> sche, const double miu_in[T]);
	Server(vector<int> sche, const double miu_in[T], int idx);
	~Server();
	int sche[T];
	double miu[T];
	double miu_now;
	bool busy;	int prior;	double service_time;
	bool cross_period(int t);
	void finish_service();
	void start_service(double t);
	void cinprior(int pr);
	int index;

private:
	void update_service_time(double t);

};


class Model_Sim
{
public:
	Model_Sim();
	~Model_Sim();
	Server servers[TotalNumberOfServers];
	list<Patient> patients1, patients2, patients3, patients4, patients5;	void run();//FCFS
	void print_result();
	void print_miu();
	void print_shift();


public:
	double av_q[T][3];
	double av_q_all[T];
	double av_w[T][3];	double av_A[T];
	double av_u[T][3];	double av_u_ser[T][3];	double av_W_Sun[T][3];	double av_W_Sun_ser[T][3];	double av_SL[T][3];	double max_lambda;
	int q0, q1, q2, q3, q4;	int q_all;
	double generate_arrival(double t, int period);
};



class OperSimu {
public:
	void initial_shift(Shifts shiftsEnv, SolutionOnShift solutionOnShift);

	void generate_alpha()
	{
		for (int i = 0; i < T; i++)
		{
			alpha[i] = lamda[i] / (lamda[i] + lamda1[i] + lamda2[i]);
			alpha1[i] = lamda1[i] / (lamda[i] + lamda1[i] + lamda2[i]);
		}

	}

	bool compare(const Patient& x, const Patient& y) {
		return x.arrival_time < y.arrival_time;
	}


	void sim_main(Shifts shiftEnv, SolutionOnShift solutionOnShift, ServiceLevel& SL_h, ServiceLevel& SL_m, ServiceLevel& SL_l)
	{
		for (int indexOfPeriod = 0; indexOfPeriod < NumberOfPeriodsCondsidered; indexOfPeriod++) {
			lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod + 1];
			lamda1[indexOfPeriod] = Lambda_initial[2][indexOfPeriod + 1];
			lamda2[indexOfPeriod] = Lambda_initial[3][indexOfPeriod + 1];
			miu[indexOfPeriod] = MiuForEachPeriod[indexOfPeriod + 1];
		}


		initial_shift(shiftEnv, solutionOnShift);
		//cout << "test1" << endl;

		/*
		cout << "show the arrival rate of the high priority customers\n";
		for (int indexOfPeriod = 0; indexOfPeriod < NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << lamda[indexOfPeriod] << " ";
		}
		cout << endl;
		cout << "show the arrival rate of the mid priority customers\n";
		for (int indexOfPeriod = 0; indexOfPeriod < NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << lamda1[indexOfPeriod] << " ";
		}
		cout << endl;
		cout << "show the arrival rate of the low priority customers\n";
		for (int indexOfPeriod = 0; indexOfPeriod < NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << lamda2[indexOfPeriod] << " ";
		}
		cout << endl;
		cout << "show the service rate of the low priority customers\n";
		for (int indexOfPeriod = 0; indexOfPeriod < NumberOfPeriodsCondsidered; indexOfPeriod++) {
			cout << miu[indexOfPeriod] << " ";
		}
		cout << endl;*/


		generate_alpha();
		Model_Sim Model;
		//Model.print_miu();
		//Model.print_shift();
		cout << "test2" << endl;

		Model.run();
		cout << "test3" << endl;

		Model.print_result();
		cout << "test4" << endl;

		for (int t = 0; t < T; t++)
		{
			//Model.av_u[t][0] /= double((num_sim - num_warm));
			//Model.av_u[t][1] /= double((num_sim - num_warm));
			//Model.av_u[t][2] /= double((num_sim - num_warm));
			SL_h.setServiceLevel(t + 1, Model.av_SL[t][0]);
			SL_m.setServiceLevel(t + 1, Model.av_SL[t][1]);
			SL_l.setServiceLevel(t + 1, Model.av_SL[t][2]);

		}

	}

};



#endif // !SIMU_MHZ

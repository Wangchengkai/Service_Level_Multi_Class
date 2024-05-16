
#pragma once
#include "simu_mahazhou.h"


//int sch[num_shift][T] = {
//{ 1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0 },
//{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 },
//{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 } };

int num_shift = TotalNumberOfServers;
double alpha[T], alpha1[T];
int sch[TotalNumberOfServers][T];

double lamda[T];
double lamda1[T];
double lamda2[T];
double miu[T];


Server::Server()
{
}

Server::Server(vector<int> sche_in, const double miu_in[T])
{
	for (int t = 0; t < T; t++)
	{
		this->sche[t] = sche_in[t];
		this->miu[t] = miu_in[t];
	}
	this->busy = false;
	this->prior = 3;
	miu_now = miu[0];
	service_time = std::numeric_limits<double>::max();

}

Server::Server(int sche_in[T], const double miu_in[T])
{
	for (int t = 0; t < T; t++)
	{
		this->sche[t] = sche_in[t];
		this->miu[t] = miu_in[t];
	}
	this->busy = false;
	this->prior = 3;
	miu_now = miu[0];
	service_time = std::numeric_limits<double>::max();

}

//Server::Server(int sche_in[T], double miu_in[T], int idx)
//{
//	for (int t = 0; t < T; t++)
//	{
//		this->sche[t] = sche_in[t];
//		this->miu[t] = miu_in[t];
//	}
//	this->busy = false;
//	this->prior = 3;
//	miu_now = miu[0];
//	service_time = std::numeric_limits<double>::max();
//	this->index = idx;
//}

Server::~Server()
{
}

void Server::update_service_time(double t)
{
	service_time = t + EXP(1 / miu_now);
}

void Server::cinprior(int pr)
{
	prior = pr;
}

bool operator<(const Server& s1, const Server& s2)
{
	return s1.service_time < s2.service_time; }

bool Server::cross_period(int t){
	bool result = false;
	int next_t;
	if (t == T - 1) next_t = 0;
	else next_t = t + 1;
	if (this->sche[t] == 0)
	{
		if (this->sche[next_t] == 1)
		{
			this->busy = false;
			this->miu_now = miu[next_t];
			this->service_time = service_time = std::numeric_limits<double>::max();
		}
	}
	else
	{
		if (this->sche[next_t] == 1)
		{
			this->miu_now = miu[next_t];
			if (this->busy)
			{
				if (next_t != 0) this->service_time = next_t * deta + EXP(1 / miu_now);
				else this->service_time = T * deta + EXP(1 / miu_now);			}
		}
		else
		{
			if (this->busy) result = true;
			this->busy = false;
			this->miu_now = miu[next_t];
			this->service_time = service_time = std::numeric_limits<double>::max();

		}
	}
	return result;
}

void Server::finish_service()
{
	this->busy = false;
	this->service_time = service_time = std::numeric_limits<double>::max();
}

void Server::start_service(double t)
{
	this->busy = true;
	update_service_time(t);
}

void OperSimu::initial_shift(Shifts shiftsEnv, SolutionOnShift solutionOnShift) {
	for (int indexOfShift = 0; indexOfShift < TotalNumberOfServers; indexOfShift++) {
		for (int t = 0; t < T; t++) {
			sch[indexOfShift][t] = 0;
		}
	}
	
	
	int countShiftNums = 0;
	for (int indexOfShift = 1; indexOfShift <= shiftsEnv.totalNumOfShifts; indexOfShift++) {
		countShiftNums += solutionOnShift.getServerOnShift(indexOfShift);
	}
	num_shift = countShiftNums;

	int countSave = 0;
	for (int indexOfShift = 1; indexOfShift <= shiftsEnv.totalNumOfShifts; indexOfShift++) {
		int countTemp = solutionOnShift.getServerOnShift(indexOfShift);
		while (countTemp >= 1) {
			for (int t = shiftsEnv.shifts[indexOfShift].startTime; t < shiftsEnv.shifts[indexOfShift].endTime; t++)
				sch[countSave][t] = 1;
			countSave++;
			countTemp--;
		}
	}
	//for (int i = 0; i < num_shift; i++) {
	//	cout << "show Sch:" << endl;
	//	for (int t = 0; t < T; t++) {
	//		cout << sch[i][t] << "\t";
	//	}
	//	cout << endl;
	//}
}



Model_Sim::Model_Sim()
{
	for (int ser = 0; ser < num_shift; ser++)
	{
		servers[ser] = Server(sch[ser], miu);
		for (int t = 0; t < T; t++) cout << servers[ser].sche[t] << " ";
		cout << endl;
	}

	q0 = q1 = q2 = q3 = q4 = q_all = 0;
	for (int t = 0; t < T; t++)
	{
		av_q[t][0] = 0;
		av_q_all[t] = 0;
		av_A[t] = 0;
		av_w[t][0] = 0;
		av_W_Sun[t][0] = 0;
		av_u[t][0] = 0;
		av_W_Sun_ser[t][0] = 0;
		av_u_ser[t][0] = 0;
		av_SL[t][0] = 0;
		av_q[t][1] = 0;
		av_w[t][1] = 0;
		av_W_Sun[t][1] = 0;
		av_u[t][1] = 0;
		av_W_Sun_ser[t][1] = 0;
		av_u_ser[t][1] = 0;
		av_SL[t][1] = 0;
		av_q[t][2] = 0;
		av_w[t][2] = 0;
		av_W_Sun[t][2] = 0;
		av_u[t][2] = 0;
		av_W_Sun_ser[t][2] = 0;
		av_u_ser[t][2] = 0;
		av_SL[t][2] = 0;
	}

	//max_lambda = -1;
	//for (int t = 0; t < T; t++)
	//{
	//	if (lamda[t] > max_lambda) max_lambda = lamda[t];
	//}
}


Model_Sim::~Model_Sim()
{
}

void Model_Sim::run()
{



	clock_t startTime, endTime;
	startTime = clock();	double t = 0, t_arrival;
	sort(servers, servers + num_shift);
	for (int d = 0; d < num_sim; d++)
	{
		//if (d % 10000 == 0) cout << double(d) / num_sim << endl;
		t = 0;
		for (int period = 0; period < T; period++)
		{
			int arri_temp = 0;
						if (!patients1.empty() || !patients2.empty() || !patients3.empty() || !patients4.empty() || !patients5.empty())
			{
				vector<int> idx_idle_sers;
				for (int s = 0; s < num_shift; s++)
				{
					if (!servers[s].busy && servers[s].sche[period] == 1)
					{
						idx_idle_sers.push_back(s);
					}
				}
				while (idx_idle_sers.size() > 0 && (!patients1.empty()))
				{
					int idx_choosen_server = rand() % idx_idle_sers.size();
					servers[idx_idle_sers[idx_choosen_server]].start_service(t);
					servers[idx_idle_sers[idx_choosen_server]].cinprior(1);
					if (d >= num_warm) {
						double temp_waiting_time = t - patients1.front().arrival_time + patients1.front().days * T * deta;
						av_W_Sun[patients1.front().period][0] += t - patients1.front().arrival_time + patients1.front().days * T * deta;
						av_u[patients1.front().period][0] += 1;
						av_W_Sun_ser[period][0] += t - patients1.front().arrival_time + patients1.front().days * T * deta;
						av_u_ser[period][0] += 1;
						if (temp_waiting_time <= tor_SeveralCustomers[1]) av_SL[patients1.front().period][0] += 1;
					}
					patients1.pop_front();
					idx_idle_sers.erase(idx_idle_sers.begin() + idx_choosen_server);
					this->q0--;
				}
				while (idx_idle_sers.size() > 0 && (!patients4.empty()))
				{
					int idx_choosen_server = rand() % idx_idle_sers.size();
					servers[idx_idle_sers[idx_choosen_server]].start_service(t);
					servers[idx_idle_sers[idx_choosen_server]].cinprior(2);
					patients4.pop_front();
					idx_idle_sers.erase(idx_idle_sers.begin() + idx_choosen_server);
					this->q3--;
				}
				while (idx_idle_sers.size() > 0 && (!patients5.empty()))
				{
					int idx_choosen_server = rand() % idx_idle_sers.size();
					servers[idx_idle_sers[idx_choosen_server]].start_service(t);
					servers[idx_idle_sers[idx_choosen_server]].cinprior(3);
					patients5.pop_front();
					idx_idle_sers.erase(idx_idle_sers.begin() + idx_choosen_server);
					this->q4--;
				}
				while (idx_idle_sers.size() > 0 && (!patients2.empty()))
				{
					int idx_choosen_server = rand() % idx_idle_sers.size();
					servers[idx_idle_sers[idx_choosen_server]].start_service(t);
					servers[idx_idle_sers[idx_choosen_server]].cinprior(2);
					if (d >= num_warm) {
						double temp_waiting_time = t - patients2.front().arrival_time + patients2.front().days * T * deta;
						av_W_Sun[patients2.front().period][1] += t - patients2.front().arrival_time + patients2.front().days * T * deta;
						av_u[patients2.front().period][1] += 1;
						av_W_Sun_ser[period][1] += t - patients2.front().arrival_time + patients2.front().days * T * deta;
						av_u_ser[period][1] += 1;
						if (temp_waiting_time <= tor_SeveralCustomers[2]) av_SL[patients2.front().period][1] += 1;
					}
					patients2.pop_front();
					idx_idle_sers.erase(idx_idle_sers.begin() + idx_choosen_server);
					this->q1--;
				}

				while (idx_idle_sers.size() > 0 && (!patients3.empty()))
				{
					int idx_choosen_server = rand() % idx_idle_sers.size();
					servers[idx_idle_sers[idx_choosen_server]].start_service(t);
					servers[idx_idle_sers[idx_choosen_server]].cinprior(3);
					if (d >= num_warm) {
						double temp_waiting_time = t - patients3.front().arrival_time + patients3.front().days * T * deta;
						av_W_Sun[patients3.front().period][2] += t - patients3.front().arrival_time + patients3.front().days * T * deta;
						av_u[patients3.front().period][2] += 1;
						av_W_Sun_ser[period][2] += t - patients3.front().arrival_time + patients3.front().days * T * deta;
						av_u_ser[period][2] += 1;
						if (temp_waiting_time <= tor_SeveralCustomers[3]) av_SL[patients3.front().period][2] += 1;
					}
					patients3.pop_front();
					idx_idle_sers.erase(idx_idle_sers.begin() + idx_choosen_server);
					this->q2--;
				}

				sort(servers, servers + num_shift);
			}


			t_arrival = generate_arrival(t, period);
			while (t < ((period + 1) * deta))
			{
				if (patients1.size() + patients2.size() + patients3.size() > 1000) { cout<<"warining1:\t"<< patients1.size() + patients2.size() + patients3.size() <<endl; system("pause"); }
				if (t_arrival < servers[0].service_time)				{
					double ra;
					ra = rand() % (1000) / (float)(1000);
					//cout << ra << endl;
					if (t_arrival > ((period + 1) * deta))
					{
						if (d >= num_warm)
						{
							this->av_w[period][0] += this->q0 * ((period + 1) * deta - t);
							this->av_w[period][1] += this->q1 * ((period + 1) * deta - t);
							this->av_w[period][2] += this->q2 * ((period + 1) * deta - t);
						}
						t = (period + 1) * deta;
						continue;
					}
					arri_temp++;
					if (d >= num_warm)
					{
						this->av_w[period][0] += this->q0 * (t_arrival - t);
						this->av_w[period][1] += this->q1 * (t_arrival - t);
						this->av_w[period][2] += this->q2 * (t_arrival - t);
					}
					if (ra >= alpha[period] + alpha1[period]) 					{
						vector<int> idx_idle_ser;
						for (int s = 0; s < num_shift; s++)
						{
							if (!servers[s].busy && servers[s].sche[period] == 1)
							{
								idx_idle_ser.push_back(s);
							}
						}
						if (idx_idle_ser.size() > 0)
						{
							int idx_choosen_server = rand() % idx_idle_ser.size();
							servers[idx_idle_ser[idx_choosen_server]].start_service(t_arrival);
							servers[idx_idle_ser[idx_choosen_server]].cinprior(3);
							if (d >= num_warm)
							{
								av_u[period][2] += 1;
								av_u_ser[period][2] += 1;
								av_SL[period][2] += 1;
							}
							sort(servers, servers + num_shift);
						}
						else
						{
							patients3.push_back(Patient(period, t_arrival));
							q2++;
						}
						this->q_all++;
						t = t_arrival;
						t_arrival = generate_arrival(t, period);
					}
					else if ((ra >= alpha[period]) && (ra < alpha1[period] + alpha[period]))					{
						vector<int> idx_idle_ser;
						for (int s = 0; s < num_shift; s++)
						{
							if (!servers[s].busy && servers[s].sche[period] == 1)
							{
								idx_idle_ser.push_back(s);
							}
						}
						if (idx_idle_ser.size() > 0)
						{
							int idx_choosen_server = rand() % idx_idle_ser.size();
							servers[idx_idle_ser[idx_choosen_server]].start_service(t_arrival);
							servers[idx_idle_ser[idx_choosen_server]].cinprior(2);
							if (d >= num_warm)
							{
								av_u[period][1] += 1;
								av_u_ser[period][1] += 1;
								av_SL[period][1] += 1;
							}
							sort(servers, servers + num_shift);
						}
						else
						{

							patients2.push_back(Patient(period, t_arrival));
							q1++;

						}
						this->q_all++;
						t = t_arrival;
						t_arrival = generate_arrival(t, period);
					}
					else 					{
						vector<int> idx_idle_ser;
						for (int s = 0; s < num_shift; s++)
						{
							if (!servers[s].busy && servers[s].sche[period] == 1)
							{
								idx_idle_ser.push_back(s);
							}
						}
						if (idx_idle_ser.size() > 0)
						{
							int idx_choosen_server = rand() % idx_idle_ser.size();
							servers[idx_idle_ser[idx_choosen_server]].start_service(t_arrival);
							servers[idx_idle_ser[idx_choosen_server]].cinprior(1);
							if (d >= num_warm)
							{
								av_u[period][0] += 1;
								av_u_ser[period][0] += 1;
								av_SL[period][0] += 1;
							}
							sort(servers, servers + num_shift);
						}
						else
						{
							vector<int> idx_idle_ser;
							for (int s = 0; s < num_shift; s++)
							{
								if (servers[s].busy && ((servers[s].prior == 2) || (servers[s].prior == 3)))
								{
									idx_idle_ser.push_back(s);
								}
							}
							if (idx_idle_ser.size() > 0)
							{
								int idx_choosen_server = rand() % idx_idle_ser.size();
								if (servers[idx_idle_ser[idx_choosen_server]].prior == 2)
								{
									patients4.push_back(Patient(period, t_arrival));
									q3++;
								}
								else if (servers[idx_idle_ser[idx_choosen_server]].prior == 3)
								{
									patients5.push_back(Patient(period, t_arrival));
									q4++;
								}
								servers[idx_idle_ser[idx_choosen_server]].start_service(t_arrival);
								servers[idx_idle_ser[idx_choosen_server]].cinprior(1);
								if (d >= num_warm)
								{
									av_u[period][0] += 1;
									av_u_ser[period][0] += 1;
									av_SL[period][0] += 1;
								}
								sort(servers, servers + num_shift);
							}
							else
							{
								patients1.push_back(Patient(period, t_arrival));
								q0++;
							}

						}
						this->q_all++;
						t = t_arrival;
						t_arrival = generate_arrival(t, period);
					}
				}
				else
				{
					if (servers[0].service_time > ((period + 1) * deta))
					{
						if (d >= num_warm)
						{
							this->av_w[period][0] += this->q0 * ((period + 1) * deta - t);
							this->av_w[period][1] += this->q1 * ((period + 1) * deta - t);
							this->av_w[period][2] += this->q2 * ((period + 1) * deta - t);
						}
						t = (period + 1) * deta;
						continue;
					}
					if (d >= num_warm)
					{
						this->av_w[period][0] += this->q0 * (servers[0].service_time - t);
						this->av_w[period][1] += this->q1 * (servers[0].service_time - t);
						this->av_w[period][2] += this->q2 * (servers[0].service_time - t);
					}
					t = servers[0].service_time;
					servers[0].finish_service();
					this->q_all--;
					if (!patients1.empty())
					{
						servers[0].start_service(t);
						servers[0].cinprior(1);
						if (d >= num_warm) {
							double temp_waiting_time = t - patients1.front().arrival_time + patients1.front().days * T * deta;
							av_W_Sun[patients1.front().period][0] += t - patients1.front().arrival_time + patients1.front().days * T * deta;
							av_u[patients1.front().period][0] += 1;
							av_W_Sun_ser[period][0] += t - patients1.front().arrival_time + patients1.front().days * T * deta;
							av_u_ser[period][0] += 1;
							if (temp_waiting_time <= tor_SeveralCustomers[1]) av_SL[patients1.front().period][0] += 1;
						}
						patients1.pop_front();
						this->q0--;
					}
					else if (!patients4.empty())
					{
						servers[0].start_service(t);
						servers[0].cinprior(2);
						patients4.pop_front();
						this->q3--;
					}
					else if (!patients5.empty())
					{
						servers[0].start_service(t);
						servers[0].cinprior(3);
						patients5.pop_front();
						this->q4--;
					}
					else if (!patients2.empty())
					{
						servers[0].start_service(t);
						servers[0].cinprior(2);
						if (d >= num_warm) {
							double temp_waiting_time = t - patients2.front().arrival_time + patients2.front().days * T * deta;
							av_W_Sun[patients2.front().period][1] += t - patients2.front().arrival_time + patients2.front().days * T * deta;
							av_u[patients2.front().period][1] += 1;
							av_W_Sun_ser[period][1] += t - patients2.front().arrival_time + patients2.front().days * T * deta;
							av_u_ser[period][1] += 1;
							if (temp_waiting_time <= tor_SeveralCustomers[2]) av_SL[patients2.front().period][1] += 1;
						}
						patients2.pop_front();
						this->q1--;
					}

					else if (!patients3.empty())
					{
						servers[0].start_service(t);
						servers[0].cinprior(3);
						if (d >= num_warm) {
							double temp_waiting_time = t - patients3.front().arrival_time + patients3.front().days * T * deta;
							av_W_Sun[patients3.front().period][2] += t - patients3.front().arrival_time + patients3.front().days * T * deta;
							av_u[patients3.front().period][2] += 1;
							av_W_Sun_ser[period][2] += t - patients3.front().arrival_time + patients3.front().days * T * deta;
							av_u_ser[period][2] += 1;
							if (temp_waiting_time <= tor_SeveralCustomers[3]) av_SL[patients3.front().period][2] += 1;
						}
						patients3.pop_front();
						this->q2--;
					}

					sort(servers, servers + num_shift);
				}
			}

						if (d >= num_warm)
			{
				for (int s = 0; s < num_shift; s++)
					if (servers[s].busy && (servers[s].prior == 1))
						av_q[period][0] += 1;
				av_q[period][0] += q0;
				for (int s = 0; s < num_shift; s++)
					if (servers[s].busy && (servers[s].prior == 2))
						av_q[period][1] += 1;
				av_q[period][1] += q1 + q3;
				for (int s = 0; s < num_shift; s++)
					if (servers[s].busy && (servers[s].prior == 3))
						av_q[period][2] += 1;
				av_q[period][2] += q2 + q4;
				av_q_all[period] += q_all;
				av_A[period] += arri_temp;
			}
			for (int s = 0; s < num_shift; s++)
			{
				if (servers[s].cross_period(period))
					this->q_all--;
			}


		}

		
		for (int s = 0; s < num_shift; s++)
		{
			if (servers[s].busy) servers[s].service_time -= T * deta;
		}
		for (auto it = patients1.begin(); it != patients1.end(); it++)
		{
			(*it).days++;
		}
		for (auto it = patients2.begin(); it != patients2.end(); it++)
		{
			(*it).days++;
		}
		for (auto it = patients3.begin(); it != patients3.end(); it++)
		{
			(*it).days++;
		}
		for (auto it = patients4.begin(); it != patients4.end(); it++)
		{
			(*it).days++;
		}
		for (auto it = patients5.begin(); it != patients5.end(); it++)
		{
			(*it).days++;
		}
	}
	endTime = clock();	//cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void Model_Sim::print_result()
{
	cout << "t,av_w1,av_w2,av_w3,av_W_Sun1,av_W_Sun2,av_W_Sun3,av_SL_1,av_SL_2,av_SL_3" << endl;
	for (int t = 0; t < T; t++)
	{
		av_q[t][0] /= double(num_sim - num_warm);
		av_q[t][1] /= double(num_sim - num_warm);
		av_q[t][2] /= double(num_sim - num_warm);
		av_q_all[t] /= double(num_sim - num_warm);
		av_w[t][0] /= double(num_sim - num_warm);
		av_w[t][1] /= double(num_sim - num_warm);
		av_w[t][2] /= double(num_sim - num_warm);
		av_A[t] /= ((num_sim - num_warm) * deta);
		av_W_Sun[t][0] /= double((num_sim - num_warm));
		av_W_Sun[t][1] /= double((num_sim - num_warm));
		av_W_Sun[t][2] /= double((num_sim - num_warm));
		av_W_Sun_ser[t][0] /= double((num_sim - num_warm));
		av_W_Sun_ser[t][1] /= double((num_sim - num_warm));
		av_W_Sun_ser[t][2] /= double((num_sim - num_warm));
		av_u[t][0] /= double((num_sim - num_warm));
		av_u[t][1] /= double((num_sim - num_warm));
		av_u[t][2] /= double((num_sim - num_warm));
		av_u_ser[t][0] /= double((num_sim - num_warm));
		av_u_ser[t][1] /= double((num_sim - num_warm));
		av_u_ser[t][2] /= double((num_sim - num_warm));
		av_SL[t][0] /= double((num_sim - num_warm) * av_u[t][0]);
		av_SL[t][1] /= double((num_sim - num_warm) * av_u[t][1]);
		av_SL[t][2] /= double((num_sim - num_warm) * av_u[t][2]);
		printf("%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", t, av_q[t][0], av_q[t][1], av_q[t][2], av_W_Sun[t][0], av_W_Sun[t][1], av_W_Sun[t][2], av_W_Sun[t][0] / av_u[t][0], av_W_Sun[t][1] / av_u[t][1], av_W_Sun[t][2] / av_u[t][2], av_SL[t][0], av_SL[t][1], av_SL[t][2]);
	}
	//system("pause");
}


void Model_Sim::print_miu()
{
	cout << endl << "miu" << endl;
	for (int ser = 0; ser < num_shift; ser++)
	{
		for (int t = 0; t < T; t++)
		{
			cout << servers[ser].miu[t] << "\t";
		}
		cout << endl;
		//for (int t = 0; t < T; t++) cout << servers[ser].sche[t] << " ";
		//cout << endl;
	}
	cout << endl;
}

void Model_Sim::print_shift()
{
	cout << endl << "schedule" << endl;
	for (int ser = 0; ser < num_shift; ser++)
	{
		for (int t = 0; t < T; t++)
		{
			cout << servers[ser].sche[t] << "\t";
		}
		cout << endl;
		//for (int t = 0; t < T; t++) cout << servers[ser].sche[t] << " ";
		//cout << endl;
	}
	cout << endl;
}

double Model_Sim::generate_arrival(double t, int period)
{
	double result = t + EXP(1 / (lamda[period] + lamda1[period] + lamda2[period]));


	//thining
	//int p = period;
	//double result = t + EXP(1 / max_lambda);
	//p = int(floor(result / deta) + 0.5) % T;
	//double random = (float)rand() / RAND_MAX;//rand() % (999 + 1) / (float)(999 + 1);
	//while (random > lamda[p] / max_lambda)
	//{
	//	result += EXP(1 / max_lambda);
	//	random = (float)rand() / RAND_MAX;//rand() % (999 + 1) / (float)(999 + 1);
	//	//cout << random << endl;
	//	p = int(floor(result / deta) + 0.5) % T;
	//}



	return result;
}



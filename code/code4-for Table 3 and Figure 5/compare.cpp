
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

	//香港排ban
	for (int indexOfShift = 1; indexOfShift <= shiftsEnv.totalNumOfShifts; indexOfShift++) {
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 0 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 8)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 8 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 16)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 4);
		if (fmod(shiftsEnv.shifts[indexOfShift].startTime, onedayDivided) == 16 && fmod(shiftsEnv.shifts[indexOfShift].endTime, onedayDivided) == 0)
			solve.currentSolutionOnShift.setServerOnShift(indexOfShift, 2);

		////冬季1
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

		////冬季2
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

		////冬季3
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

		////冬季4
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
				//完成均匀化
				//计算第一类病人的服务水平
				solve.serviceLevel.computeServiceLevel_forHighPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel);
				//showServiceLevelByPeriod(solve.serviceLevel);
				solve.serviceLevelForSeveralCustomers[1] = solve.serviceLevel;
				break;

			case 2:
				//torWaitingTime = tor_SeveralCustomers[1];
				computationOfServicelLevel = solve.serviceLevel.initialComputeServiceLevel_ThreeTypes(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift);
				//完成均匀化。

				//计算第二类病人服务水平
				for (int indexOfPeriod = 0; indexOfPeriod <= NumberOfPeriodsCondsidered; indexOfPeriod++) {
					Lamda[indexOfPeriod] = Lambda_initial[1][indexOfPeriod];
				}
				torWaitingTimePL = tor_SeveralCustomers[2];
				solve.serviceLevelForLowPriorityCustomer.computeServiceLevel_forLowPriority(shiftsEnv, solve.currentSolutionOfServers, solve.currentSolutionOnShift, computationOfServicelLevel, "Mid");
				//showServiceLevelByPeriod(solve.serviceLevelForLowPriorityCustomer);
				solve.serviceLevelForSeveralCustomers[2] = solve.serviceLevelForLowPriorityCustomer;

				//计算第三类病人服务水平
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



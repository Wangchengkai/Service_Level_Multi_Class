#include "data.h"

double Lambda_initial[4][NumberOfPeriodsCondsidered + 1] = {
	{0
	},{0,
	//0.194,0.226,0.323,0.161,0.161,0.194,0.065,0.419,0.516,0.968,0.871,0.710,0.581,0.677,0.452,0.613,1.000,0.548,0.677,0.548,0.742,0.613,0.419,0.516
	},{0,
	//2.742,2.129,1.452,2.097,1.516,1.710,1.355,2.484,3.935,6.710,7.742,6.677,8.129,7.032,6.452,7.355,6.581,6.194,6.516,5.613,5.355,5.129,4.742,3.516
	},{0,
	
	}

};


double Lamda[NumberOfPeriodsCondsidered + 1] ={ 0,//dayOne


};
double LamdaMP[NumberOfPeriodsCondsidered + 1] ={ 0,

};
double LamdaLP[NumberOfPeriodsCondsidered + 1] ={ 0,

};

double torWaitingTime = 10.0 / 60.0;
double torWaitingTimePL = 0.5;

double Miu = 4;//(5.91 * 95416 / (95416 + 33401));
double MiuForEachPeriod[NumberOfPeriodsCondsidered + 1] ={ 1,
	Miu,Miu,Miu,Miu,Miu,Miu, Miu,Miu,//Miu,Miu,Miu,Miu, Miu,Miu,Miu,Miu,Miu,Miu, Miu,Miu,Miu,Miu,Miu,Miu
};


double partA = 0.4;//0.00934 + 0.02195;
double partB = 0.3;// 0.33406;
double partC = 0.3;// 0.63464;

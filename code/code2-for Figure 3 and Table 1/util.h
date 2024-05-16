#include <cmath>
#include <iostream>
using namespace std;

//Phi(x,n)
double Poi(double x, int n)
{
	if (n < 0)
		return 0;
	double rtn = 0;
	double y = exp(-x);
	for (int j = 0; j <= n; ++j)
	{
		rtn += y;
		y *= x / (j + 1);
	}
	return rtn;
}

double BetaDist(double x, int a, int b)
{
	double y, z;
	int i;
	y = pow(x, a);
	z = a * y * (1 - x);
	for (i = 2; i <= b; i++)
	{
		y += z;
		z *= (a + i - 1) * (1 - x) / i;
	}
	return (y);
}

//factorial
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

double HyperGeo(int _delta, int _u, int q, int _s)
{
	if (_s - q < _u - _delta|| q < _delta){
		// if (q < _delta)
		// 	cout << "fucked1" << endl;
		// else
		// 	cout << "fucked2" << endl;
		return 0;
	}
	else
	{
		return fac(q) * fac(_s - q) * fac(_u) * fac(_s - _u) / fac(_delta) / fac(q - _delta) / fac(_u - _delta) / fac(_s - q - _u + _delta) / fac(_s);
	}
}
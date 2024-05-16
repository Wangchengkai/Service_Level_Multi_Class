//precise calculation of model3 SL1

#include "dataForXJ.h"
#include "util.h"
#include "functionsForXJ.h"
#include <numeric>
#include <algorithm>

using namespace std;

double mat[N + 1][N + 1]; double piq[N][2 * qmax];
double pi[2 * qmax], npi[2 * qmax], newPi[2 * qmax];

/*
i: state
t: time interval
arr_t: arrival time epoch
*/
double betaQ(int i, int t, double arr_t, int p[], int offshift[])
{
    double at, rtn;
    if (arr_t < 1 - tao)
    {
        at = double(p[t]) * double(mu) * double(tao);
        rtn = Poi(at, i - p[t]);
    }
    else
    {
        //int pt2 = (t == T - 1) ? p[t] : p[t + 1];
        int pt2 = (t == T - 1) ? p[0] : p[t + 1];        at = p[t] * mu * (1 - arr_t) + pt2 * mu * (arr_t + tao - ceil(tao));
        rtn = Poi(at, i - pt2 - offshift[t]);
    }
    return rtn;
}

double phiQ(double t, int p1, int p2)
{
    return pow(t, p1) * pow(1 - t, p2);
}



void Uniformization(int offshift[],int p[],double SL[]){
    double lamda[T] = { 0 };
    for (int t = 1; t <= T; t++)
        lamda[t - 1] = Lamda[t];
    //memcpy(lamda, Lamda, sizeof(lamda));

  /*  cout << "offshift\t" << endl;
    for (int i = 0; i < T; i++)
        cout << "offshift[" << i << "]=\t" << offshift[i] << endl;

    cout << "p\t" << endl;
    for (int i = 0; i < T; i++)
        cout << "p[" << i << "]=\t" << p[i] << endl;

    cout << "lamda\t" << endl;
    for (int i = 0; i < T; i++)
        cout << "lamda[" << i << "]=\t" << lamda[i] << endl;*/


    //system("pause");

    int t, j, n, i;
    for (j = 0; j <= N; j++)
    {
        mat[j][j] = 1;
        for (i = j + 1; i <= N; i++)
        {
            if (j == 0)
                mat[i][j] = 1;
            else
                mat[i][j] = mat[i - 1][j - 1] + mat[i - 1][j];
        } 
    }
#define muu(x) (MIN(p[t], x) * mu)

    double w[T], u[T], q[T], sum_u, sum_w, sum_q, W[2 * qmax], x, y, z, zz, std_q[T];
    double Fq1[2 * qmax], Fq2[2 * qmax], Rq[2 * qmax], Qq[2 * qmax];// , SL[T];
    double bq[2 * qmax], cq[2 * qmax];
    double A[N], D[N], J[N], H[N], fN[N];
    int delta, a, b;
    double Beta;
    for (t = 0, x = 0; t < T; t++)
        x = MAX(lamda[t] + p[t] * mu, x);
    for (i = 0, y = exp(-x); i < N; i++)
    { /*Poisson distri of mean gamma*/
        fN[i] = y;
        y *= x / (double)(i + 1); //avoid calculating factorial!
    }

    for (i = 1, D[0] = fN[0]; i < N; i++)
        D[i] = D[i - 1] + fN[i];

    for (i = 0; i < N; i++)
    {
        for (n = i + 1, J[i] = 0; n < N; n++)
            J[i] += (double)(i + 1) * fN[n] / (n + 1);
        z = 1 - tao;
        a = i + 1;
        Beta = pow(z, a);
        zz = a * Beta * (1 - z);
        for (n = i + 1, A[i] = 0; n < N; n++)
        {
            A[i] += fN[n] * Beta;
            Beta += zz;
            b = n - i + 1;
            zz *= (a + b - 1) * (1 - z) / b;
        }
        z = 1 - tao;
        a = i + 2;
        Beta = pow(z, a);
        zz = a * Beta * (1 - z);
        for (n = i + 1, H[i] = 0; n < N; n++)
        {
            H[i] += (i + 1) * fN[n] / (n + 1) * Beta;
            Beta += zz;
            b = n - i + 1;
            zz *= (a + b - 1) * (1 - z) / b;
        }
    }


    for (j = 0; j < 2 * qmax; j++)
    {
        pi[j] = 0; //avoid using matrix
    }
    pi[0] = 1.0;

    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < 2 * qmax; j++)
            piq[i][j] = 0;
    }

        int stopBool = 0;
    double piLastTime[2 * qmax] = { 0 };

    
    while (!stopBool) {
        for (j = 0; j < 2 * qmax; j++) {
            piLastTime[j] = pi[j]; //avoid using matrix
        }

        for (t = 0; t < T; t++)
        {
            //std::cout << "time: " << t << std::endl;
            w[t] = u[t] = q[t] = 0;
            for (j = 0; j <= qmax + 1; j++)
                newPi[j] = W[j] = Fq1[j] = Fq2[j] = Rq[j] = Qq[j]  = 0;
            n = 0;

            y = exp(-x);

            do
            {
                                //if (n % 49 == 0)
                    //std::cout << n << std::endl;
                for (j = 0; j <= qmax; j++)
                {

                    piq[n][j] = pi[j];

                }


                for (j = 0; j <= qmax; j++)
                {
                    W[j] += pi[j];
                    newPi[j] += pi[j] * y;
                    npi[j] = 1 / x * (pi[j] * (x - muu(j) - lamda[t]) + pi[j + 1] * muu(j + 1) + ((j == 0) ? 0 : pi[j - 1] * lamda[t]));
                }
                for (j = 0; j <= qmax; j++)
                {
                    pi[j] = npi[j];
                    w[t] += y * W[j] * MAX(j - p[t], 0) * period / (n + 1);
                }
                n++;
                y = y * x / n;
            } while (n < N);


            for (j = 0, z = 0; j <= qmax; j++)
                z += newPi[j];
            for (j = 0, zz = 0; j <= qmax; j++)
            {
                pi[j] = newPi[j] / z;

                q[t] += j * pi[j];
                zz += j * j * pi[j];
            }


                        if (offshift[t] > 0)
            {
                std::fill(pi, pi + qmax, 0);
                for (j = 0; j < qmax; j++)
                {
                    if (j >= p[t]) {
                        pi[j - offshift[t]] += newPi[j] / z;
                    }
                    else {
                        for (int delta = MAX(0, j + offshift[t] - p[t]); delta <= MIN(offshift[t], j); delta++) {
                            pi[j - delta] += HyperGeo(delta, offshift[t], j, p[t]) * newPi[j] / z;
                        }
                    }
                }
            }


        }
        //cout << "initialstate Prob \t" << endl;
        //for (int q = 0; q < 5; q++) {
        //    cout << "q: " << pi[q] << endl;
        //}
        
        double tempSum = 0;
        for (j = 0; j < 2 * qmax; j++) {
            tempSum += abs(piLastTime[j] - pi[j]); //avoid using matrix
            //if (j <= qmax)
            //    cout << "temp" << j << "\t" << pi[j] << endl;
        }


        //cout << "GAP\t" << tempSum << endl;
        if (tempSum <= 1.0e-5){
            stopBool = 1;
        }
        //system("pause");
    }
    cout << endl;

    //*****************************************************end*************************************************

    
        
    //T*N*qmax*N*N
    for (t = 0; t < T; t++)
    {
        //std::cout << "time: " << t << std::endl;
        w[t] = u[t] = q[t] = 0;
        for (j = 0; j <= qmax + 1; j++)
            newPi[j] = W[j] = Fq1[j] = Fq2[j] = Rq[j] = Qq[j] = bq[j] = cq[j] = 0;
        n = 0;

        y = exp(-x);


        do
        {
            //if (n % 49 == 0)
                //std::cout << n << std::endl;
            for (j = 0; j <= qmax; j++)
            {
                // Fq1[j] += pi[j] * A[n] / x;
                // Fq2[j] += pi[j] * (1 - D[n] - A[n]) / x;
                piq[n][j] = pi[j];

                double temp2 = 0;
                double temp3 = 0;
                double bet = betaQ(j, t, (1 - tao) / N, p, offshift);
                double bet2[N + 1];
                for (int k = 0; k <= N; ++k)
                    bet2[k] = betaQ(j, t, 1 - tao + k * tao / N, p, offshift);
                for (i = 0; i <= n; ++i)
                {
                    // integral
                    double inte2 = 0;
                    double inte3 = 0;
                    for (int k = 0; k <= N; ++k)
                    {
                        if (k == 0 || k == N)
                        {
                            inte2 += 0.5 * phiQ(k * (1 - tao) / N, i, n - i) * bet;
                            inte3 += 0.5 * phiQ((1 - tao) + k * tao / N, i, n - i) * bet2[k];
                        }
                        else
                        {
                            inte2 += phiQ(k * (1 - tao) / N, i, n - i) * bet;
                            inte3 += phiQ((1 - tao) + k * tao / N, i, n - i) * bet2[k];
                        }
                    }

                    temp2 += inte2 * mat[n][i] * piq[i][j];
                    temp3 += inte3 * mat[n][i] * piq[i][j];
                }
                bq[j] += temp2 * (1 - tao) / N * y;
                cq[j] += temp3 * tao / N * y;
            }


            for (j = 0; j <= qmax; j++)
            {
                W[j] += pi[j];
                newPi[j] += pi[j] * y;
                npi[j] = 1 / x * (pi[j] * (x - muu(j) - lamda[t]) + pi[j + 1] * muu(j + 1) + ((j == 0) ? 0 : pi[j - 1] * lamda[t]));
            }
            for (j = 0; j <= qmax; j++)
            {
                pi[j] = npi[j];
                w[t] += y * W[j] * MAX(j - p[t], 0) * period / (n + 1);
            }
            n++;
            y = y * x / n;
        } while (n < N);




        for (j = 0, z = 0; j <= qmax; j++)
            z += newPi[j];
        for (j = 0, zz = 0; j <= qmax; j++)
        {
            pi[j] = newPi[j] / z;
            // piq[0][j] = pi[j];
            q[t] += j * pi[j];
            zz += j * j * pi[j];
        }

        std_q[t] = sqrt(zz - q[t] * q[t]);
        u[t] = lamda[t] - q[t] + ((t == 0) ? 0 : q[t - 1]);
        double sumBQ = 0, sumCQ = 0;
        // fstream test;
        // test.open("test.csv", ios::in | ios::out | ios::app);
        // test << "time: " << t << endl;
        // test << "bq,cq" << endl;
        for (i = 0, SL[t] = 1; i <= qmax; i++)
        {
            // test << bq[i] << ',' << cq[i] << '\n';
            sumBQ += bq[i];
            sumCQ += cq[i];
            SL[t] -= bq[i] + cq[i];
        }
        // test << "\n\n";
        // test.close();
        //cout << "sumBQ:\t" << sumBQ << endl;
        //cout << "sumCQ:\t" << sumCQ << endl;

                // if (delta > 0)
        // {
        //     for (i = 0; i <= p[t]; i++)
        //     {
        //         int ii = MAX(0, i - delta);
        //         for (j = ii; j <= MIN(i, p[t] - delta); j++)
        //         {
        //             pi[j] += HyperGeo(i - j, delta, i, p[t]) * newPi[i] / z;
        //         }
        //     }
        //     for (i = p[t] + 1; i < qmax; i++)
        //     {
        //         pi[i - delta] = newPi[i] / z;
        //     }
        // }

        // double sum1 = 0;
        // sum1 = std::accumulate(newPi, newPi + qmax, sum1);
        // cout << sum1 / z << endl;

        // new 
        if (offshift[t] > 0)
        {
            std::fill(pi, pi+qmax, 0);
            for (j = 0; j < qmax; j++)
            {
                if (j >= p[t]){
                    pi[j-offshift[t]] += newPi[j] / z;
                } else {
                    for (int delta = MAX(0, j + offshift[t] - p[t]); delta <= MIN(offshift[t], j); delta++) {
                        pi[j-delta] += HyperGeo(delta, offshift[t], j, p[t]) * newPi[j] / z;
                    }
                }
            }
        }
        // double sum = 0;
        // sum = std::accumulate(pi, pi + qmax, sum);
        // cout << sum << endl;
        // for (j = 0; j < qmax; j++)
        //     pi[j] /= sum;
    }
    //for (t = 0; t < T; t++)
        //cout << t + 1 << "," << SL[t] << endl;
       
    
}

#include "time.h"
#include "iostream"
#include "fstream"
#include"cmath"
using namespace std;
int main()
	{

	const int DX = 5,NX = 1000,PML = 50,NX_PML = NX + 2 * PML,favg = 20,T = 210;
	int V[NX_PML];
 	const float dt = 0.001, pi = 3.1415926;
	for (int i = 0; i < NX_PML; i++)
		{
		V[i] = 2000;
		}
	const int xn = (NX_PML / 2);
	float wavelet[152];
	for (int j = 0; j < 152; j++)
		{
		int i = j - 62;
		wavelet[j] = 4 * (1 - 2 * (pi*favg*dt*(i))*pi*favg*dt*(i))*exp(-(pi*favg*dt*(i))*(pi*favg*dt*(i)));
		}

	float A1_x[NX_PML];
	float A2_x[NX_PML];
	float U1_x[NX_PML];
	float U2_x[NX_PML];
	float U[NX_PML];
	for (int i = 0; i < NX_PML; i++)
		{
		A1_x[i] = 0;
		A2_x[i] = 0;
		U1_x[i] = 0;
		U2_x[i] = 0;
		U[i] = 0;
		}
	clock_t t1, t2;
	t1 = clock();
#pragma acc data present_or_copyin(V[0:NX_PML]) present_or_copyin(wavelet[0:152]) copyin(U2_x[0:NX_PML]) copyin(A2_x[0:NX_PML]) copyin(A1_x[0:NX_PML]) copyout(U[0:NX_PML]) 
	{
	for (int n = 1; n < T; n++)
			{
			if (n < 152)
				{
				U1_x[xn] = wavelet[n];
				}
#pragma acc  parallel loop
#pragma acc data copyin(U1_x[xn])
          for (int i = 2; i < NX_PML - 2; i++)
				{
				A2_x[i] = A1_x[i]- dt / DX*(U1_x[i + 1] - U1_x[i]);
				}
#pragma acc  parallel loop
          for (int i = 2; i < NX_PML - 2; i++)
				{
                                U2_x[i] = U1_x[i]- dt*V[i] * V[i] / DX*(A2_x[i] - A2_x[i - 1]);
                                }
#pragma acc parallel loop
	for (int s = 0; s < NX_PML; s++)
				{
				A1_x[s] = A2_x[s];
				U1_x[s] = U2_x[s];
				U[s] = U2_x[s];
				}          
          
                          }
            }
	t2 = clock();
	cout<<t2-t1<<endl;
	//快照保存
	ofstream fout("/home/echo/桌面/test.txt");
      for (int s = 0; s < NX_PML; s++)
		{
		fout<<U[s]<<endl;
		}
	  fout.close();
	//getchar();
	return 0;
	}

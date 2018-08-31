#include<stdio.h>
#include<stdlib.h>
#include <iostream>
#include <complex>
#include <cmath>
#include<stdio.h>
using namespace std;
/*declearation*/
complex<float>*** allocate3floatcomplex(int page, int row, int column);
void free3floatcomplex(complex<float>*** p, int page, int row);
float** allocate2float(int row, int column);
void free2float(float** p, int row);
float*** allocate3float(int page, int row, int column);
void free3float(float*** p, int page, int row);
/*extern function*/
extern "C"
void ffdExtrapolation(complex<float>*** Data0, float*** VelModel, float ***PhaseFactor, float* Vmin,
	int Nx, int Ny, int Nz, int Nt,
	int Nw, int  Nw1, int  Nw2, int  Nw3, int  Nw4,
	float dkx, float dky, float dw,
	int StepNum, float StepDz, int NumDepth,
	float** BeamImageXoY, float** BeamImageXoZ, float** BeamImageYoZ,
	int SxGrid, int SyGrid, int wave_len);
/*main*/
int main(int argc,char**argv){
	int Nx = 800; int Ny = 450; int Nz = 300; int Nt = 3000;
	int Nw=300; int  Nw1=1; int  Nw2=5; int  Nw3=90; int  Nw4=90;
	float dkx=0.1; float dky=0.1; float dw=0.003;
	int StepNum=60; float StepDz=50; int NumDepth=60;
	float** BeamImageXoY; float** BeamImageXoZ; float** BeamImageYoZ;
	int SxGrid=128; int SyGrid=128; int wave_len=128;
	float*** VelModel; float ***PhaseFactor; float* Vmin;
	/*local varibale*/
	int i,j,l;
	complex<float>*** Data0=allocate3floatcomplex(Nx,Ny,Nt);
	BeamImageXoY = allocate2float(Nx, Ny);
	BeamImageXoZ = allocate2float(Nx, Nz);
	BeamImageYoZ = allocate2float(Ny, Nz);
	VelModel = allocate3float(Nx, Ny, Nz);
	PhaseFactor = allocate3float(Nx, Ny, Nz);
	Vmin = (float*)malloc(Nz * sizeof(float));
	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (l = 0; l < Nz; l++) {
				VelModel[i][j][l] = 2000;
				PhaseFactor[i][j][l] = 2;
				Vmin[l] = 2000;
			}
		}
	}
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(l=0;l<Nt;l++){
				Data0[i][j][l] = (sin(i*0.01+j)*cos(l*0.134+l),0);
			}
		}
	}
	ffdExtrapolation(Data0,VelModel, PhaseFactor,Vmin,
	 Nx,  Ny,  Nz,  Nt,
	 Nw,   Nw1,   Nw2,   Nw3,   Nw4,
	 dkx, dky, dw,
	 StepNum,  StepDz,  NumDepth,
	 BeamImageXoY,  BeamImageXoZ,  BeamImageYoZ,
	 SxGrid,  SyGrid,  wave_len);


	return 0;

}
complex<float>*** allocate3floatcomplex(int page, int row, int column)
{
	int i, j, k;
	complex<float>*** p;
	p = new complex<float>**[page];
	for (i = 0; i < page; i++)
	{
		p[i] = new complex<float>*[row];
		for (j = 0; j < row; j++)
		{
			p[i][j] = new complex<float>[column];
		}
	}
	if (p == 0)
	{
		free3floatcomplex(p, page, row);
	}
	else
	{
		for (i = 0; i < page; i++)
		{
			for (j = 0; j < row; j++)
			{
				for (k = 0; k < column; k++)
				{
					p[i][j][k] = complex<float>(0,0);
				}
			}
		}
	}
	return p;
}
void free3floatcomplex(complex<float>*** p, int page, int row)
{
		int i, j;
			for (i = 0; i < page; i++)
					{
								for (j = 0; j < row; j++)
											{
															delete[] p[i][j];
																	}
									}
				for (i = 0; i < page; i++)
						{
									delete[] p[i];
										}
					delete[] p;
						p = 0;
}

float** allocate2float(int row, int column)
{
	float** p;
	p = new float*[row];
	for (int i = 0; i<row; i++)
	{
		p[i] = new float[column];
	}
	if (p == 0)
	{
		free2float(p, row);
	}
	else
	{
		for (int i = 0; i<row; i++)
		{
			for (int j = 0; j<column; j++)
			{
				p[i][j] = 0.0;
			}
		}
	}
	return p;
}

float*** allocate3float(int page, int row, int column)
{
	int i, j, k;
	float*** p;
	p = new float**[page];
	for (i = 0; i<page; i++)
	{
		p[i] = new float*[row];
		for (j = 0; j<row; j++)
		{
			p[i][j] = new float[column];
		}
	}
	if (p == 0)
	{
		free3float(p, page, row);
	}
	else
	{
		for (i = 0; i<page; i++)
		{
			for (j = 0; j<row; j++)
			{
				for (k = 0; k<column; k++)
				{
					p[i][j][k] = 0.0;;
				}
			}
		}
	}
	return p;
}

void free2float(float** p, int row)
{
	int i;
	for (i = 0; i<row; i++)
	{
		delete[] p[i];
	}
	delete[] p;
	p = 0;
}

void free3float(float*** p, int page, int row)
{
	int i, j;
	for (i = 0; i<page; i++)
	{
		for (j = 0; j<row; j++)
		{
			delete[] p[i][j];
		}
	}
	for (i = 0; i<page; i++)
	{
		delete[] p[i];
	}
	delete[] p;
	p = 0;
}

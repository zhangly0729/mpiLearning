#pragma once
#include"iostream"
#include"complex"
#include"cuda_runtime.h"
#include"cuda_device_runtime_api.h"
#include"device_launch_parameters.h"
#include"cufft.h"  
#include"stdio.h"
using namespace std;

#define PI 3.14159265354
extern "C"
void ffdExtrapolation(complex<float>*** Data0, float*** VelModel, float ***PhaseFactor, float* Vmin,
	int Nx, int Ny, int Nz, int Nt,
	int Nw, int  Nw1, int  Nw2, int  Nw3, int  Nw4,
	float dkx, float dky, float dw,
	int StepNum, float StepDz, int NumDepth,
	float** BeamImageXoY, float** BeamImageXoZ, float** BeamImageYoZ,
	int SxGrid, int SyGrid, int wave_len);


extern "C"
void ffdExtrapolation2d(complex<float>** Data0, float** VelModel, float **PhaseFactor, float* Vmin,
	int Nx, int Nz, int Nt,
	int Nw, int  Nw1, int  Nw2, int  Nw3, int  Nw4,
	float dkx, float dw,
	int StepNum, float StepDz, float** BeamImage,int wave_len);

static __device__ __host__ inline cuComplex ComplexMul(cuComplex a, cuComplex b)
{
	cuComplex c;
	c.x = a.x * b.x - a.y * b.y;
	c.y = a.x * b.y + a.y * b.x;
	return c;
}

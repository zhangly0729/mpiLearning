#include"lib.h"
#include <iostream>
#include <complex>
#include <cmath>
#include<stdio.h>

/*
	检查cuda错误函数
*/
#define Check(CALL) \
{\
	if (CALL != cudaSuccess)\
	{\
		cout << "行号:" << __LINE__ << endl;\
		cout << "错误:" << cudaGetErrorString(CALL) << endl;\
	}\
}

/********************************************
****************　辅助核函数　　*************
*********************************************/

/*
	用于反FFT变换后数据的归一化
*/
__global__ void normalizing2d(cufftComplex* data, int data_len, int value)
{
	unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
	for (int i = idx* blockDim.x; i < (idx + 1)* blockDim.x; i++)
	{
		if (i < data_len)
		{
			data[i].x /= value;
			data[i].y /= value;
		}
	}

}
__global__ void normalizing3d(cufftComplex* data, int data_len, int value)
{
	unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
	for (int i = idx* blockDim.x; i < (idx + 1)* blockDim.x; i++)
	{
		if (i < data_len)
		{
			data[i].x /= value;
			data[i].y /= value;
		}
	}

}
__global__ void normalizing3d(float* data, int data_len, int value)
{
	unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
	for (int i = idx* blockDim.x; i < (idx + 1)* blockDim.x; i++)
	{
		if (i < data_len)
		{
			data[i] /= value;
		}
	}

}
/*
	汉明窗-提取有效频段数据
*/
__global__ void HAMMING_Window2d(cufftComplex* data, cufftComplex* data0, int Nw1, int Nw2, int Nw3, int Nw4, int Nw, int LT, int Nx)
{
	unsigned int ix = threadIdx.x + blockDim.x*blockIdx.x;

	unsigned int idx = LT*ix;

	unsigned int idw = Nw*ix;

	float Hammingw = 0;
	if (ix < Nx)
	{
		for (int iw = 0; iw < LT; iw++)
		{
			if (iw >= Nw1&&iw <= Nw2)
			{
				Hammingw = 0.54f + 0.46f*cos(PI*(iw - Nw1) / (Nw2 - Nw1) - PI);
				data[idx + iw].x = data[idx + iw].x * Hammingw;
				data[idx + iw].y = data[idx + iw].y * Hammingw;
			}
			else if (iw >= Nw3&&iw <= Nw4)
			{
				Hammingw = 0.54f + 0.46f*cos(PI*(Nw3 - iw) / (Nw4 - Nw3) - PI);
				data[idx + iw].x = data[idx + iw].x * Hammingw;
				data[idx + iw].y = data[idx + iw].y * Hammingw;
			}
			else if (iw<Nw1 || iw>Nw4)
			{
				data[idx + iw].x = 0;
				data[idx + iw].y = 0;
			}
		}

		for (int iw = Nw1; iw < Nw4 + 1; iw++)
		{
			data0[idw + iw - Nw1] = data[idx + iw];
		}
	}
}
__global__ void HAMMING_Window3d(cufftComplex* data, cufftComplex* data0, int Nw1, int Nw2, int Nw3, int Nw4, int Nw, int LT, int Nx, int Ny)
{
	unsigned int ix = threadIdx.x + blockDim.x*blockIdx.x;
	unsigned int iy = threadIdx.y + blockDim.y*blockIdx.y;
	unsigned int idx = LT*(iy*Ny + ix);
	unsigned int idw = Nw*(iy*Ny + ix);
	float Hammingw = 0;
	if (ix < Ny&&iy < Nx)
	{
		for (int iw = 0; iw < LT; iw++)
		{
			if (iw >= Nw1&&iw <= Nw2)
			{
				Hammingw = 0.54f + 0.46f*cos(PI*(iw - Nw1) / (Nw2 - Nw1) - PI);
				data[idx + iw].x = data[idx + iw].x * Hammingw;
				data[idx + iw].y = data[idx + iw].y * Hammingw;
			}
			else if (iw >= Nw3&&iw <= Nw4)
			{
				Hammingw = 0.54f + 0.46f*cos(PI*(Nw3 - iw) / (Nw4 - Nw3) - PI);
				data[idx + iw].x = data[idx + iw].x * Hammingw;
				data[idx + iw].y = data[idx + iw].y * Hammingw;
			}
			else if (iw<Nw1 || iw>Nw4)
			{
				data[idx + iw].x = 0;
				data[idx + iw].y = 0;
			}
		}

		for (int iw = Nw1; iw < Nw4 + 1; iw++)
		{
			data0[idw + iw - Nw1] = data[idx + iw];
		}
	}
}


/*
	将有效频段数据还原
*/
__global__ void HAMMING_Window_Inverse2d(cufftComplex* data, cufftComplex* data0, int Nw1, int Nw2, int Nw3, int Nw4, int Nw, int LT, int Nx)
{
	unsigned int ix = threadIdx.x + blockDim.x*blockIdx.x;
	unsigned int idx = LT*ix;
	unsigned int idw = Nw*ix;
	if (ix < Nx)
	{
		for (int iw = Nw1; iw < Nw4 + 1; iw++)
		{
			data[idx + iw] = data0[idw + iw - Nw1];
		}
	}
}
__global__ void HAMMING_Window_Inverse3d(cufftComplex* data, cufftComplex* data0, int Nw1, int Nw2, int Nw3, int Nw4, int Nw, int LT, int Nx, int Ny)
{
	unsigned int ix = threadIdx.x + blockDim.x*blockIdx.x;
	unsigned int iy = threadIdx.y + blockDim.y*blockIdx.y;
	unsigned int idx = LT*(iy*Ny + ix);
	unsigned int idw = Nw*(iy*Ny + ix);
	if (iy < Nx&&ix < Ny)
	{
		for (int iw = Nw1; iw < Nw4 + 1; iw++)
		{
			data[idx + iw] = data0[idw + iw - Nw1];
		}
	}
}


/*
	相移校正
*/
__global__ void PS2d(cufftComplex* Data, float* Vmin, int iw, int Nx, int Nw, int iz, float w, float dkx, float StepDz)
{
	/*按Kx波数划分区域*/
	int ix = threadIdx.x + blockDim.x*blockIdx.x;
	int idw = Nw*ix;
	float kx = 0;
	ix < Nx / 2 ? kx = ix*dkx : kx = -(Nx - ix)*dkx;/*波数转换*/


	cuComplex Factor;
	float ARG = 1.0f - kx*kx*Vmin[iz] * Vmin[iz] / (w*w);/*计算相移因子*/
	if (ARG <= 0)
	{
		ARG = 0.0f;
		Factor = make_cuComplex(0.0f, 0.0f);
	}
	else
	{
		Factor = make_cuComplex(cos(w / Vmin[iz] * sqrtf(ARG)*StepDz), sin(w / Vmin[iz] * sqrtf(ARG)*StepDz));
	}

	if (ix < Nx)
	{
		cuComplex d = Data[idw + iw];
		Data[idw + iw] = ComplexMul(d, Factor);
	}

}
__global__ void PS3d(cufftComplex* Data, float* Vmin, int iw, int Nx, int Ny, int Nw, int iz, float w, float dkx, float dky, float StepDz)
{
	/*按Kx,Ky波数划分区域*/
	int ix = threadIdx.x + blockDim.x*blockIdx.x;
	int iy = threadIdx.y + blockDim.y*blockIdx.y;
	int idw = Nw*(iy*Ny + ix);
	float kx = 0;
	float ky = 0;

	iy < Nx / 2 ? kx = iy*dkx : kx = -(Nx - iy)*dkx;/*波数转换*/
	ix < Ny / 2 ? ky = ix*dky : ky = -(Ny - ix)*dky;/*波数转换*/

	cuComplex Factor;
	float ARG = 1.0f - (kx*kx + ky*ky)*Vmin[iz] * Vmin[iz] / (w*w);/*计算相移因子*/
	if (ARG <= 0)
	{
		ARG = 0.0f;
		Factor = make_cuComplex(0.0f, 0.0f);
	}
	else
	{
		Factor = make_cuComplex(cos(w / Vmin[iz] * sqrtf(ARG)*StepDz), sin(w / Vmin[iz] * sqrtf(ARG)*StepDz));
	}

	if (ix < Ny&&iy < Nx)
	{
		cuComplex d = Data[idw + iw];
		Data[idw + iw] = ComplexMul(d, Factor);
	}
}
/*
	时移校正
*/
__global__ void SSF2d(cufftComplex* Data, float* PhaseFactor, float* Vmin, int iw, int Nx, int Nz, int Nw, int iz, float w, float dkx, float StepDz)
{
	/*按Kx,Ky波数划分区域*/
	int ix = threadIdx.x + blockDim.x*blockIdx.x;

	int idw = Nw*ix;
	int idz = Nz*ix;

	float kx = 0;

	ix < Nx / 2 ? kx = ix*dkx : kx = -(Nx - ix)*dkx;/*波数转换*/

	if (ix < Nx)
	{
		float ARG = w*PhaseFactor[idz + iz];
		cufftComplex Factor = make_cuComplex(cos(ARG*StepDz), sin(ARG*StepDz));
		cuComplex d = Data[idw + iw];
		Data[idw + iw] = ComplexMul(d, Factor);
	}
}

__global__ void SSF3d(cufftComplex* Data, float* PhaseFactor, float* Vmin, int iw, int Nx, int Ny, int Nz, int Nw, int iz, float w, float dkx, float dky, float StepDz)
{
	/*按Kx,Ky波数划分区域*/
	int ix = threadIdx.x + blockDim.x*blockIdx.x;
	int iy = threadIdx.y + blockDim.y*blockIdx.y;

	int idw = Nw*(iy*Ny + ix);
	int idz = Nz*(iy*Ny + ix);

	float kx = 0;
	float ky = 0;

	iy < Nx / 2 ? kx = iy*dkx : kx = -(Nx - iy)*dkx;/*波数转换*/
	ix < Ny / 2 ? ky = ix*dky : ky = -(Ny - ix)*dky;/*波数转换*/


	if (ix < Ny&&iy < Nx)
	{
		float ARG = w*PhaseFactor[idz + iz];
		cufftComplex Factor = make_cuComplex(cos(ARG*StepDz), sin(ARG*StepDz));
		cuComplex d = Data[idw + iw];
		Data[idw + iw] = ComplexMul(d, Factor);
	}
}


/*
三维傅里叶有限差分延拓-接口函数
*/
void ffdExtrapolation(complex<float>*** Data0, float*** VelModel, float ***PhaseFactor, float* Vmin,
	int Nx, int Ny, int Nz, int Nt,
	int Nw, int  Nw1, int  Nw2, int  Nw3, int  Nw4,
	float dkx, float dky, float dw,
	int StepNum, float StepDz,int NumDepth,
	float** BeamImageXoY,float** BeamImageXoZ,float** BeamImageYoZ,
	int SxGrid,int SyGrid,int wave_len)
{
			
	/*
	检测显存是否越界，重置处理！
	*/
	size_t freeMem, totalMem; cudaMemGetInfo(&freeMem, &totalMem);
	float free_device_memory = freeMem / (1024.0 * 1024.0 * 1024.0);
	float need_device_memory = Nx*Ny*Nt * sizeof(float)+Nx*Ny*(Nt/2+1)*sizeof(cufftComplex) +Nx*Ny*Nw*sizeof(cufftComplex)+
			Nx*Ny*Nz * 2 * sizeof(float) + (Nz) * sizeof(float);
	need_device_memory /= (1024.0 * 1024.0 * 1024.0);//GB
	size_t worksize = 0;
	cufftEstimate3d(Nx, Ny, Nt, CUFFT_R2C, &worksize); 
	need_device_memory+=worksize/(1024.0 * 1024.0 * 1024.0);
	//二维平面FFT
	const int rank = 2;			 
	int n[rank] = { Nx, Ny };
	int inembed[3] = { Nw ,Ny ,Nx }; // 输入数据的[页数，列数，行数](3维)；[列数，行数]（2维）
	int onembed[3] = { Nw, Ny,Nx }; // 输出数据的[页数，列数，行数]；[列数，行数]（2维）
	int istride = Nw; // 每个输入信号相邻两个元素的距离
	int idist = 1;   // 每两个输入信号第一个元素的距离
	int ostride = Nw; // 每个输出信号相邻两个元素的距离
	int odist = 1;   // 每两个输出信号第一个元素的距离
	int batch = Nw;  // 进行 fft 的信号个数
	cufftEstimateMany(rank,n,inembed,istride,idist, onembed,  ostride,odist, CUFFT_C2C, batch,&worksize);
	need_device_memory+=worksize/(1024.0 * 1024.0 * 1024.0);

	if (need_device_memory > 15.8)/*是否超过全局显存*/
	{
		printf("Thr program need %fG memory,It's out of device Memory!", need_device_memory);
		cin.get();
		exit(0);
	}
	if ((free_device_memory - need_device_memory) < 0.1)/*是否超过当前可用显存,门限值0.1G显存*/
	{
		cudaDeviceReset();//重置设备，以保证足够的显存可用。
	}
	/*设置由设备1执行相关处理*/
	//cudaSetDevice(0);
	//cudaDeviceReset();
	/*
	主机内存申请
	*/
	float* host_Data0; Check(cudaMallocHost((void**)&host_Data0, Nx*Ny*Nt * sizeof(float)));
	float* host_VelModel;	Check(cudaMallocHost((void**)&host_VelModel, Nx*Ny*Nz * sizeof(float)));
	float* host_PhaseFactor; Check(cudaMallocHost((void**)&host_PhaseFactor, Nx*Ny*Nz * sizeof(float)));
	float* host_Vmin = Vmin;	/*维度相同，直接传递指针即可*/
	/*
	主机数据维度转换，3D-->1D
	*/
#pragma omp parallel for //将三维结构转为一维结构
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nt; k++)
			{
				host_Data0[i*(Ny*Nt) + j*Nt + k] = Data0[i][j][k].real();
			}
			for (int k = 0; k < Nz; k++)
			{
				host_VelModel[i*(Ny*Nz) + j*Nz + k] = VelModel[i][j][k];
				host_PhaseFactor[i*(Ny*Nz) + j*Nz + k] = PhaseFactor[i][j][k];
			}
		}

	/*
	设备内存申请
	*/
	float* device_Data0_real; Check(cudaMalloc((void**)&device_Data0_real, Nx*Ny*Nt*sizeof(float)));
	cufftComplex* device_Data0; Check(cudaMalloc((void**)&device_Data0, Nx*Ny*(Nt/2+1)*sizeof(cufftComplex)));
	cufftComplex* device_dataW; Check(cudaMalloc((void**)&device_dataW, Nx*Ny*Nw * sizeof(cufftComplex)));/*存储主频段数据*/
	float* device_VelModel; Check(cudaMalloc((void**)&device_VelModel, Nx*Ny*Nz * sizeof(float)));
	float* device_PhaseFactor; Check(cudaMalloc((void**)&device_PhaseFactor, Nx*Ny*Nz * sizeof(float)));
	float* device_Vmin; Check(cudaMalloc((void**)&device_Vmin, Nz * sizeof(float)));

	/*
	初始化设备内存 : Host-->Device
	*/
	printf("Initialization*******\n");
	Check(cudaMemcpy(device_Data0_real,host_Data0, Nx*Ny*Nt * sizeof(float),
		cudaMemcpyHostToDevice));
	Check(cudaMemcpy(device_VelModel, host_VelModel, Nx*Ny*Nz * sizeof(float),
		cudaMemcpyHostToDevice));
	Check(cudaMemcpy(device_PhaseFactor, host_PhaseFactor, Nx*Ny*Nz * sizeof(float),
		cudaMemcpyHostToDevice));
	Check(cudaMemcpy(device_Vmin, host_Vmin, Nz * sizeof(float),
		cudaMemcpyHostToDevice));

	/*创建cufftplan3D及cufftplanmany句柄。*/
	cufftHandle cufft3DHandle;
	cufftResult_t cufftresult;
	cufftresult = cufftPlan3d(&cufft3DHandle, Nx, Ny, Nt, CUFFT_R2C);
	/*
	执行fft正变换  3D (x,y,t)-->(kx,ky,w)
	*/
	cufftExecR2C(cufft3DHandle, device_Data0_real, device_Data0);
	cufftDestroy(cufft3DHandle);
	cufftresult = cufftPlan3d(&cufft3DHandle, Nx, Ny, Nt, CUFFT_C2R);
	/*
//	销毁cuFFT句柄
//	*/
//	cufftDestroy(cufft3DHandle);

	//二维plan_fft
	cufftHandle plan_Nxyfft_many;
	cufftresult=cufftPlanMany(&plan_Nxyfft_many, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch);

	/*
	提取给定频率范围的数据
	*/
	dim3 g1((Ny + 15) / 16, (Nx + 15) / 16);
	dim3 b1(16, 16);
	HAMMING_Window3d << <g1, b1 >> > (device_Data0, device_dataW, Nw1, Nw2, Nw3, Nw4, Nw,(Nt/2+1), Nx, Ny);

	/*傅里叶有限差分递推*/
	dim3 g2((Ny + 15) / 16, (Nx + 15) / 16);
	dim3 b2(16, 16);
	for (int iz = 0; iz < StepNum; iz++)
	{
		printf("iz=%d\n", iz);
		/*Step1:相移校正,FK域*/
		for (int iw = 0; iw < Nw; iw++)
		{
			float w;
			(iw + Nw1) < Nt / 2 ? w = (iw + Nw1)*dw : w = -(Nt - ((iw + Nw1)))*dw;/*频率转换*/
			if (w == 0) w = 0.0000001;
			PS3d << <g2, b2 >> > (device_dataW, device_Vmin, iw, Nx, Ny, Nw, iz, w, dkx, dky, StepDz);
		}

		/*x-o-y平面反FFT变换 (x,y,w)-->(kx,ky,w)*/
		cufftExecC2C(plan_Nxyfft_many, device_dataW, device_dataW, CUFFT_INVERSE); // 执行 cuFFTplanmany，反变换
		dim3 b3 = 512;
		dim3 g3 = ((Nx*Ny*Nw) + b3.x* b3.x - 1) / (b3.x* b3.x);
		normalizing3d << <g3, b3 >> > (device_dataW, Nx*Ny*Nw, Nx*Ny);

		/*Step2:时移校正,FX域*/
		for (int iw = 0; iw < Nw; iw++)
		{
			float w;
			(iw + Nw1) < Nt / 2 ? w = (iw + Nw1)*dw : w = -(Nt - ((iw + Nw1)))*dw;/*频率转换*/
			if (w == 0) w = 0.0000001;

			SSF3d << <g2, b2 >> > (device_dataW, device_PhaseFactor, device_Vmin, iw, Nx, Ny, Nz, Nw, iz, w, dkx, dky, StepDz);
		}

		/*x-o-y平面反FFT变换 (kx,ky,w)-->(x,y,w)*/
		cufftExecC2C(plan_Nxyfft_many, device_dataW, device_dataW, CUFFT_FORWARD); // 执行 cuFFTplanmany，正变换
		
	//	/*
	//	将给定频率范围的数据反代
	//	*/
	//	HAMMING_Window_Inverse3d << <g1, b1 >> > (device_Data0, device_dataW, Nw1, Nw2, Nw3, Nw4, Nw,( Nt/2+1), Nx, Ny);
	//	/*
	//	执行fft反变换  3D (kx,ky,w)-->(x,y,t)
	//	*/
	//	cufftExecC2R(cufft3DHandle, device_Data0, device_Data0_real);
	//	/*
	//	数据传输，:Device-->Host
	//	*/
	//	dim3 b4(512);
	//	dim3 g4(((Nx*Ny*Nt) + b4.x* b4.x - 1) / (b4.x* b4.x));
	//	normalizing3d<< <g4, b4 >> > (device_Data0_real, Nx*Ny*Nt, Nx*Ny*Nt);
	//	Check(cudaMemcpy(host_Data0, device_Data0_real, Nx*Ny*Nt * sizeof(float), cudaMemcpyDeviceToHost));

	//	/*
	//	主机数据维度转换，1D-->3D
	//	*/
	//	//#pragma omp parallel for //将三维结构转为一维结构
	//	for (int i = 0; i < Nx; i++)
	//		for (int j = 0; j < Ny; j++)
	//			for (int k = 0; k < Nt; k++)
	//			{
	//				Data0[i][j][k] = (host_Data0[i*(Ny*Nt) + j*Nt + k], 0);
	//			}


	//	//保存XoY结果
	//	if (iz ==  NumDepth)
	//	{
	//		for (int i = 0; i < Nx; i++)
	//			for (int j = 0; j < Ny; j++)
	//			{
	//				BeamImageXoY[i][j] = Data0[i][j][wave_len].real();
	//			}
	//	}
	//	
	//	//保存XoZ结果
	//	for (int i = 0; i < Nx; i++)
	//	{
	//		BeamImageXoZ[i][iz] = Data0[i][SyGrid][wave_len].real();
	//	}
	//	//保存YoZ结果
	//	for (int j = 0; j <Ny; j++)
	//	{
	//		BeamImageYoZ[j][iz] = Data0[SxGrid][j][wave_len].real();
	//	}
	}

	/*
	释放主机/设备内存
	*/
	if (host_Data0 != 0)			cudaFreeHost(host_Data0);
	if (host_VelModel != 0)			cudaFreeHost(host_VelModel);
	if (host_PhaseFactor != 0)		cudaFreeHost(host_PhaseFactor);

	if (device_Data0 != 0)			cudaFree(device_Data0);
	if (device_dataW != 0)			cudaFree(device_dataW);
	if (device_VelModel != 0)		cudaFree(device_VelModel);
	if (device_PhaseFactor != 0)	cudaFree(device_PhaseFactor);
	if (device_Vmin != 0)			cudaFree(device_Vmin);
	
	cufftDestroy(cufft3DHandle);
	cufftDestroy(plan_Nxyfft_many);
}

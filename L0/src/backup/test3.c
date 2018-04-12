/*************************************************************************
  > File Name: test.c
  > Author: ma6174
  > Mail: ma6174@163.com 
  > Created Time: Sat 07 Apr 2018 02:45:59 AM PDT
 ************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{


	int rank,size;
	int left,right;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	float starttime=MPI_Wtime();

	int N=2;
	float* dataS=(float*)malloc(N*sizeof(float));
	float *dataR;
	if(rank==0)
	{
		dataR=(float*)malloc(size*N*sizeof(float));
	}
	char numstring[10];
	int tag;	

	right=rank+1;
	left=rank-1;

	int i=0;
	for(i=0;i<N;i++)
	{
		dataS[i]=rank+i;
	}

	MPI_Gather(dataS,N,MPI_FLOAT,dataR,N,MPI_FLOAT,0,MPI_COMM_WORLD);
	if(rank==0)
		{

			for(i=0;i<size*N;i++)
				printf("%f\n",dataR[i]);
		}
		float endtime=MPI_Wtime();
	if(rank==0)
		printf("Total time is %.6f\n",(endtime-starttime));
	MPI_Finalize();
	return 0;
}

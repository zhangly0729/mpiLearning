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

	int N=100;
	float* dataS=(float*)malloc(N*sizeof(float));
	float* dataR=(float*)malloc(N*sizeof(float));

	char numstring[10];
	int tag;	

	right=rank+1;
	left=rank-1;

//	do{
		if(rank==0)
		{
			int i=0;
			for(i=0;i<N;i++)
			{
				dataS[i]=1.01;
			}

		}
//	} while (dataS[0]>0);
	MPI_Bcast(dataS,N,MPI_FLOAT,0,MPI_COMM_WORLD);
	printf("%f\n",dataS[0]);
	float endtime=MPI_Wtime();
	if(rank==0)
		printf("Total time is %.6f\n",(endtime-starttime));
	MPI_Finalize();
	return 0;
}

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

	if(rank==0)
	{
		int i=0;
		for(i=0;i<N;i++)
		{
			dataS[i]=1.01;
		}
		//Right Exchange
		sprintf(numstring,"%d%d0001",rank,right);	
		tag=atoi(numstring);
		MPI_Send(dataS,N,MPI_FLOAT,right,tag,MPI_COMM_WORLD);
	    sprintf(numstring,"%d%d0001",right,rank);
		tag=atoi(numstring);
		MPI_Recv(dataR,N,MPI_FLOAT,right,tag,MPI_COMM_WORLD,&status);	
	}
	else if(rank==size-1)
	{
		int i=0;
		for(i=0;i<N;i++)
		{
			dataS[i]=3.03;
		}
		//Left Exchange
		sprintf(numstring,"%d%d0001",left,rank);	
		tag=atoi(numstring);
		MPI_Recv(dataR,N,MPI_FLOAT,left,tag,MPI_COMM_WORLD,&status);
		sprintf(numstring,"%d%d0001",rank,left);	
		tag=atoi(numstring);
		MPI_Send(dataS,N,MPI_FLOAT,left,tag,MPI_COMM_WORLD);
	}
	else
	{
		if(rank%2==1)
		{
		//Left Exchage
		sprintf(numstring,"%d%d0001",left,rank);	
		tag=atoi(numstring);
		MPI_Recv(dataR,N,MPI_FLOAT,left,tag,MPI_COMM_WORLD,&status);
		
		printf("%f\n",dataR[0]);
		sprintf(numstring,"%d%d0001",rank,left);	
		tag=atoi(numstring);
		MPI_Send(dataS,N,MPI_FLOAT,left,tag,MPI_COMM_WORLD);

		//Right Exchange
		sprintf(numstring,"%d%d0001",rank,right);	
		tag=atoi(numstring);
		MPI_Send(dataS,N,MPI_FLOAT,right,tag,MPI_COMM_WORLD);
	    sprintf(numstring,"%d%d0001",right,rank);
		tag=atoi(numstring);
		MPI_Recv(dataR,N,MPI_FLOAT,right,tag,MPI_COMM_WORLD,&status);	
		printf("%f\n",dataR[0]);
		}
		else
		{


		}
	}
	float endtime=MPI_Wtime();
	if(rank==0)
		printf("Total time is %.6f\n",(endtime-starttime));
	MPI_Finalize();
	return 0;
}

/*************************************************************************
    > File Name: pack.c
    > Author: echo
    > Mail: iai.wu@outlook.com 
    > Created Time: Thu 12 Apr 2018 04:57:35 AM PDT
 ************************************************************************/

#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
int main(int argc,char* argv[])
{
	int i,j,k,rank,size,position;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	char buffer[200];

	if(rank==0)
	{
		i=1;
		j=2;
		MPI_Pack(&i,1,MPI_INT,buffer,200,&position,MPI_COMM_WORLD);
		MPI_Pack(&j,1,MPI_INT,buffer,200,&position,MPI_COMM_WORLD);	
		MPI_Send(&position,1,MPI_INT,1,0,MPI_COMM_WORLD);
		MPI_Send(buffer,position,MPI_PACKED,1,1,MPI_COMM_WORLD);
	}
	if(rank!=0)
	{
		MPI_Recv(&position,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(buffer,position,MPI_PACKED,0,1,MPI_COMM_WORLD,&status);
		position=0;
		MPI_Unpack(buffer,200,&position,&i,1,MPI_INT,MPI_COMM_WORLD);
		MPI_Unpack(buffer,200,&position,&j,1,MPI_INT,MPI_COMM_WORLD);
		printf("%d\n",position);
		printf("%d\n",i);
		printf("%d\n",j);
	}
	MPI_Finalize(); 
}	

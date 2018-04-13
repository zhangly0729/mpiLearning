/*************************************************************************
    > File Name: type.c
    > Author: echo
    > Mail: iai.wu@outlook.com 
    > Created Time: Thu 12 Apr 2018 05:47:27 PM PDT
 ************************************************************************/

#include<stdio.h>
#include<mpi.h>

typedef struct {
	float a;
	int i;
	double k;
} DATA_INTEGRAL;

int main(int argc,char* argv[])
{
	int i,j,k,size,rank;

	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	float a[4][4]={1,2,3,4,
					2,3,3,5,
					5,5,1,5,
					7,8,5,2};
	float b[4];


	MPI_Datatype rowtype;
	MPI_Type_contiguous(4,MPI_FLOAT,&rowtype);
	MPI_Type_commit(&rowtype);

	if(rank==0)
	{
	
		//MPI_Send(&data,1,DATA_INTEGRAL,1,0,MPI_COMM_WORLD); 
		MPI_Send(a[2],1,rowtype,1,0,MPI_COMM_WORLD);
	}
	if(rank!=0)
	{
		//MPI_Recv(&data,1,DATA_INTEGRAL,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(b,1,rowtype,0,0,MPI_COMM_WORLD,&status);
		printf("b0=%f\n",b[0]);
	}

	MPI_Type_free(&rowtype);
	MPI_Finalize();
	return 0;
}

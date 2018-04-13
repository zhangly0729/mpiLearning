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

	float a[16]={1,2,3,4,
				 2,3,3,5,
				 5,5,1,5,
				 7,8,5,2};
	float b[16];

	int blocks=2;
	int blocklength[2]={1,3};
	int blockdisplacment[2]={0,10};

	MPI_Datatype indextype;
	MPI_Type_indexed(2,blocklength,blockdisplacment,MPI_FLOAT,&indextype);
	MPI_Type_commit(&indextype);

	if(rank==0)
	{

		MPI_Send(a,1,indextype,1,0,MPI_COMM_WORLD);
	}
	if(rank!=0)
	{
		MPI_Recv(b,4,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
		for(i=0;i<4;i++)
			printf("b=%f\n",b[i]);
	}

	MPI_Type_free(&indextype);
	MPI_Finalize();
	return 0;
}

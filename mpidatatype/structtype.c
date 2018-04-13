/*************************************************************************
    > File Name: structtype.c
    > Author: echo
    > Mail: iai.wu@outlook.com 
    > Created Time: Thu 12 Apr 2018 07:35:48 PM PDT
 ************************************************************************/

#include<stdio.h>
#include<mpi.h>
typedef struct {
	float x;
	int j;
	double a;
} MYTYPE;
int main(int argc,char* argv[])
{
	int i,j,k,rank,size;
	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MYTYPE mydata;

	MPI_Datatype myvar;
	MPI_Datatype old_type[3];
	MPI_Aint indice[3];

	old_type[0]=MPI_FLOAT;
	old_type[1]=MPI_INT;
	old_type[2]=MPI_DOUBLE;
	MPI_Address(&mydata.x,&indice[0]);
	MPI_Address(&mydata.j,&indice[1]);
	MPI_Address(&mydata.a,&indice[2]);
	indice[2]-=indice[1];
	indice[1]-=indice[0];
	indice[0]=0;
	int blocks[3]={1,1,1};

	MPI_Type_struct(3,blocks,indice,old_type,&myvar);
	MPI_Type_commit(&myvar);

	if(rank==0)
	{
		mydata.x=1.0;
		mydata.j=3;
		mydata.a=3.2;
		MPI_Send(&mydata,1,myvar,1,0,MPI_COMM_WORLD);
	}
	if(rank!=0)
	{
		MPI_Recv(&mydata,1,myvar,0,0,MPI_COMM_WORLD,&status);
		printf("x=%f\n",mydata.x);
		printf("j=%d\n",mydata.j);
		printf("a=%f\n",mydata.a);
	 }
	
	MPI_Type_free(&myvar);
	MPI_Finalize();
}

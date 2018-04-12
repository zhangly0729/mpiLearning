/*************************************************************************
    > File Name: mul.c
    > Author: echo
    > Mail: iai.wu@outlook.com 
    > Created Time: Thu 12 Apr 2018 12:18:35 AM PDT
 ************************************************************************/

#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
int main(int argc,char**argv)
{
	int i,j,k,L,D;
	int rank,size;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//paramter 
	L=20;D=L/size;
	float *a=NULL,*b=NULL,c[D];
	float *A=NULL,*C=NULL;
	a=(float*)malloc(sizeof(float)*L*D);
	b=(float*)malloc(sizeof(float)*L);

	for(i=0;i<D;i++)
		c[i]=0;
	//initial the martrix & vector
	if(rank==0)
	{
		A=(float*)malloc(sizeof(float)*L*L);
		C=(float*)malloc(sizeof(float)*L);
		for(i=0;i<L;i++)
			b[i]=1;
		for(i=0;i<L;i++)
			for(j=0;j<L;j++)
			{
			A[i*L+j]=j;
			}
	}	

	//BroadCast the vector b to child process
	MPI_Bcast(b,L,MPI_FLOAT,0,MPI_COMM_WORLD);
	//Scatter the martrix A to slave process  
	MPI_Scatter(A,D*L,MPI_FLOAT,a,D*L,MPI_FLOAT,0,MPI_COMM_WORLD);

	//caculate
	for(j=0;j<D;j++)
		for(i=0;i<L;i++)
		{
			c[j]+=a[i]*b[i];
		}
	
	//Gather the value c to result vector C in master process
	MPI_Gather(&c,D,MPI_FLOAT,C,D,MPI_FLOAT,0,MPI_COMM_WORLD);

	//print the result
	if(rank==0)
	{
		for(i=0;i<L;i++)
		{
			printf("%f\n",C[i]);
		}
	}
	MPI_Finalize();	
	return 0;
}


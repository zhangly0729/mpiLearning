//for jacobi 
#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<string.h>
int main(int argc,char** argv)
{
	int num_proce,my_id,i,j,k;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proce);

	//Parameter 
	int M=400,N=680,L=(N/num_proce+2);
	float *a=NULL,*b=NULL,*c=NULL, *result=NULL;
	size_t Totalsize=sizeof(float)*M*N;
	size_t Msize=sizeof(float)*M*L;
	MPI_Status status;
	//Initial 
	a=(float*)malloc(Msize);
	b=(float*)malloc(Msize);
	result=(float*)malloc(Totalsize);
	memset(a,0,Msize);
	memset(b,0,Msize);
	float halodata[M];

	for(i=1;i<N/num_proce+1;i++)	
	{
		a[1*L+i]=8;
		//a[60*L+i]=100;
		//a[90*L+i]=100;
	}

	int left=my_id-1;
	int right=my_id+1;
	if(my_id==0) left=MPI_PROC_NULL;
	if(my_id==num_proce-1)	right=MPI_PROC_NULL;


	//Jacobi itrration
	int t;
	for(t=0;t<20;t++)
	{	
		for(i=1;i<M-1;i++)
		{
			for(j=1;j<L-1;j++)
			{
				b[i*L+j]=0.25*(a[(i+1)*L+j]+a[(i-1)*L+j]+a[i*L+j+1]+a[i*L+j-1]);
			}
		}
		a=b; 

	//Halo Exchange
	//R --> L
	for(i=0;i<M;i++)
		halodata[i]=a[i*L+1];
	MPI_Sendrecv(halodata,M,MPI_FLOAT,left,1,halodata,M,MPI_FLOAT,right,1,MPI_COMM_WORLD,&status);
	for(i=0;i<M;i++)
		a[i*L+L-1]=halodata[i];
	//L --> R
	for(i=0;i<M;i++)
		halodata[i]=a[i*L+L-2];
	MPI_Sendrecv(halodata,M,MPI_FLOAT,right,2,halodata,M,MPI_FLOAT,left,2,MPI_COMM_WORLD,&status);
	for(i=0;i<M;i++)
		a[i*L]=halodata[i];

		printf("This is %d iteration\n",t);
	}
	//Gather the result
	for(i=0;i<M;i++)
		for(j=0;j<L-2;j++)
			a[i*(L-2)+j]=b[i*L+j+1];

	//printf("k=%d\n",k);	
	MPI_Gather(a,M*(L-2),MPI_FLOAT,result,M*(L-2),MPI_FLOAT,0,MPI_COMM_WORLD);

	printf("ok!\n");	
	//Write the result
	if(my_id==0)
	{	
		FILE* fp;
		fp=fopen("result.dat","wb");
		for(k=0;k<(num_proce);k++)
		{
			printf("k=%d\n",k);
			float *pR=&result[M*(L-2)*k];
			for(i=0;i<L-2;i++)
			{		
				for(j=0;j<M;j++)
				{
					fwrite(&pR[j*(L-2)+i],sizeof(float),1,fp);
				}
			}
		}
		fclose(fp);

		printf("Ok 2\n");
	}

	MPI_Finalize();
	return 0;
}

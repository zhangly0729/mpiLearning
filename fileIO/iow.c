/*************************************************************************
    > File Name: io.c
    > Author: echo
    > Mail: iai.wu@outlook.com 
    > Created Time: Thu 12 Apr 2018 11:59:43 PM PDT
 ************************************************************************/

#include<stdio.h>
#include"mpi.h"
#include"string.h"
int main(int argc,char* argv[])
{
	int i,j,l,rank,size;
	
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	MPI_File fp;
	MPI_File_open(MPI_COMM_WORLD,"data",MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fp);
	//MPI_File_seek(fp,rank*sizeof(int),MPI_SEEK_SET);
	MPI_File_write_at(fp,rank*sizeof(int),&rank,1,MPI_INT,&status);
	MPI_File_close(&fp);


	if(rank==0)
	{
		int a[size];
		FILE* f;
		f=fopen("data","r");
		fread(a,size,sizeof(int),f);
		for(i=0;i<size;i++)
		{
			//fscanf(f,"%d",&a[i]);
			printf("a[%d]=%d\n",i,a[i]);
		}
		fclose(f);
	}
	MPI_Finalize();
	return 0;
}

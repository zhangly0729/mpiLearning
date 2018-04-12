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

	int rankX, rankY,rankZ;
	int ndims = 3;
	int dims[3] = {
		3, 3,2}; 
	int periods[2] = {
		0, 0}; 
	int reorder = 0;

	int remainX[2] = {
		0, 0}; 
	int remainY[2] = {
		0, 1}; 
	int remainZ[2] = {
		0, 1}; 
	MPI_Comm comm3d;
	MPI_Comm commX, commY,commZ;


	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm3d);
	MPI_Cart_sub(comm3d, remainX, &commX);
	MPI_Cart_sub(comm3d, remainY, &commY);
	MPI_Cart_sub(comm3d, remainZ, &commZ);
	
	MPI_Comm_rank(commX, &rankX);
	MPI_Comm_rank(commY, &rankY);
	MPI_Comm_rank(commZ, &rankZ);

	printf("rank = %d;   X = %d;  Y = %d; Z = %d \n", rank, rankX, rankY,rankZ);

	MPI_Finalize();

	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
void main(int argc,char *argv[])
{
 int size,rank;
 MPI_Status status;
 struct var {
  int i,n,size,rank;
  float mytmp,ml,mr,mlen,l,len;
 } ;
 struct var b;
 //下面构建一个自定义变量，利用MPI函数构建，若用C语言的结构体
 //则在Send和Recv是不能识别数据类型
 MPI_Datatype myvar;
 MPI_Datatype old_types[2];
 MPI_Aint indices[2]; 
 //指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
 int blocklens[2];
 MPI_Init(&argc,&argv);
 blocklens[0]=4;
 blocklens[1]=6;
 //指定原来旧的数据类型
 old_types[0]=MPI_INT;
 old_types[1]=MPI_FLOAT;  
 //指定每个块中变量的偏移量，需要依靠一个结构体的实例，这里是b。
 MPI_Address(&b,&indices[0]);
 MPI_Address(&b.mytmp,&indices[1]);
 indices[1] -= indices[0];
 indices[0]=0;
 //看来indices[0]也可以一开始就应当赋值为0
 

 //创建新数据于var之中
 MPI_Type_struct(2,blocklens,indices,old_types,&myvar);
 //注册新数据
 MPI_Type_commit(&myvar);

 MPI_Comm_size(MPI_COMM_WORLD,&size);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 if(rank == 0)
 {
  struct var a={0};
  a.i=100;
  fprintf(stdout,"before sending the value i is %d\n",a.i);
  fflush(stdout);
  MPI_Send(&a,1,myvar,1,1,MPI_COMM_WORLD);
  fprintf(stdout,"after sending the value i is %d\n",a.i);
  fflush(stdout);
 }
 else if(rank == 1)
 { struct var t={0};
  fprintf(stdout,"before receiving the value i is %d\n",t.i);
  fflush(stdout);
  MPI_Recv(&t,1,myvar,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  fprintf(stdout,"after receiving the value i is %d\n",t.i);
  fflush(stdout);
 }

 MPI_Finalize();
}

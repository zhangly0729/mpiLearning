all :
	nvcc -c lib.cu -o lib.o
	g++  -c main.cpp -o main.o
	nvcc lib.o main.o -o a.out -lcufft 
run :

clean :
	rm *.o -rf

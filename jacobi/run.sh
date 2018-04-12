#########################################################################
# File Name: run.sh
# Author: ma6174
# mail: ma6174@163.com
# Created Time: Wed 11 Apr 2018 05:46:04 PM PDT
#########################################################################
#!/bin/bash

rm a.out
mpicc jacobi.c -O3 -o  a.out
mpiexec -n 4 ./a.out

python snap.py

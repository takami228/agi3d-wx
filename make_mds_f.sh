#!/bin/sh

gcc -O3 -I ../lapack-3.4.2/lapacke/include/ -c mds_cpu_f.c -o mds_cpu_f.o 
gcc -L ../lapack-3.4.2/ -L ../lapack-3.4.2/lapacke/ mds_cpu_f.o -o mds_cpu_f -lpthread -llapacke -llapack -lgfortran -lgoto2 -lm 

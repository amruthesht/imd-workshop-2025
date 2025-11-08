#!/bin/bash

cd ctypes

if [ ! -f libvdos.so ]; then
echo "compiling C libraries"
gcc -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
gcc -shared -lgomp -lgsl -lgslcblas vdos.o -o libvdos.so
fi

cd ..

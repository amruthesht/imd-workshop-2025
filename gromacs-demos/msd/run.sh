#!/bin/bash

echo "setting GROMACS version"
source /scratch/hcho96/Research/imd-v3/gromacs-imd-v3/bin/GMXRC

#if [ ! -f libvdos.so ]; then
#echo "compiling C libraries"
#gcc -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
#gcc -L/home/hcho96/.conda/envs/imdreader-test/lib -I/home/hcho96/.conda/envs/imdreader-test/include -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
#gcc -shared -lgomp -lgsl -lgslcblas vdos.o -o libvdos.so
#gcc -L/home/hcho96/.conda/envs/imdreader-test/lib -I/home/hcho96/.conda/envs/imdreader-test/include -shared -lgomp -lgsl -lgslcblas vdos.o -o libvdos.so
#fi

mkdir -p imd
cd imd
if [ ! -f topol.tpr ]; then
echo "running grompp"
gmx grompp -f ../imd.mdp -c ../start.gro -p ../topol.top -o >& grompp.out
fi
echo "starting mdrun"
echo ""
echo "**************************************************************"
echo "wait a few seconds, then:"
echo "run (in separate terminal): python demo.py"
echo "**************************************************************"
echo ""
# modify gromacs execution string according to available resources
# here we use 2 parallel threads
gmx mdrun -v -nt 12 -imdwait -imdport 8891
cd ..

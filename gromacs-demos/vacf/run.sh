#!/bin/bash

echo "setting GROMACS version"
source /scratch/hcho96/Research/imd-v3/gromacs-imd-v3/bin/GMXRC

if [ ! -f libvdos.so ]; then
echo "compiling C libraries"
#gcc -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
gcc -L/home/hcho96/.conda/envs/imdclient/lib -I/home/hcho96/.conda/envs/imdclient/include -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
#gcc -shared -lgomp -lgsl -lgslcblas vdos.o -o libvdos.so
gcc -L/home/hcho96/.conda/envs/imdclient/lib -I/home/hcho96/.conda/envs/imdclient/include -shared -lgomp -lgsl -lgslcblas vdos.o -o libvdos.so
fi

mkdir -p imd
cd imd
if [ ! -f topol.tpr ]; then
echo "running grompp"
#Original system
gmx grompp -f ../imd.mdp -c ../start.gro -p ../topol.top -o mda.tpr >& grompp.out
#Smaller system
#gmx grompp -f ../input-streaming.mdp -c ../test-start.gro -p ../test-topol.top -n ../index.ndx -o >& grompp.out
fi
echo "starting mdrun"
#echo ""
#echo "**************************************************************"
#echo "wait a few seconds, then:"
#echo "run (in separate terminal): python demo.py"
#echo "**************************************************************"
#echo ""
# modify gromacs execution string according to available resources
# here we use 2 parallel threads
#cp ../imd.tpr .

gmx mdrun -s mda.tpr -v -nt 12 -imdwait -imdport 8894
cd ..

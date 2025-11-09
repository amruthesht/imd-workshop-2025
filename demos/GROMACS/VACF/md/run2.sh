#!/bin/bash

echo "setting GROMACS version"
source /usr/local/gromacs-imd-v3/bin/GMXRC

mkdir -p imd
cd imd
echo "running grompp"
gmx grompp -f ../imd.mdp -c ../start.gro -p ../topol.top -o imd.tpr >& grompp.out
echo "starting mdrun"
if [ ! -z $1 ]; then port=$1
else port=8888; fi
gmx mdrun -s imd.tpr -v -nt 2 -pin on -imdwait -imdport ${port}
cd ..
rm -r imd

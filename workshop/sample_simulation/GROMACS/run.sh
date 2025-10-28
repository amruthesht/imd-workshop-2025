#!/bin/bash

# IMD Streaming Workshop - GROMACS Simulation Script
# System: Lysozyme in water (30,423 atoms)
# IMDv3 streaming enabled on port 8889

echo "================================================"
echo "IMD Streaming Workshop - GROMACS"
echo "System: Lysozyme in water"
echo "Atoms: 30,423 | Waters: 9,467 | Ions: 62"
echo "================================================"
echo ""

# Create output directory
mkdir -p output

# Preprocess topology (grompp)
if [ ! -f output/run.tpr ]; then
    echo "Preprocessing topology with grompp..."
    # Run grompp from output/ directory so imdgroup.gro is created there
    # (imdgroup.gro filename is hardcoded and cannot be changed)
    cd output
    gmx grompp -f ../input/input-streaming.mdp \
               -c ../input/start.gro \
               -p ../input/topol.top \
               -n ../input/index.ndx \
               -o run.tpr \
               -po mdout.mdp \
               -maxwarn 1
    cd ..
    
    if [ $? -ne 0 ]; then
        echo "Error: grompp failed!"
        exit 1
    fi
    echo "Preprocessing complete!"
    echo "  TPR file: output/run.tpr"
    echo "  Processed MDP: output/mdout.mdp"
    echo "  IMD group: output/imdgroup.gro"
    echo ""
fi

echo "Starting MD simulation with IMDv3 streaming..."
echo "Port: 8889"
echo ""
echo "Waiting for IMD connection from MDAnalysis..."
echo "When you see: 'IMD: Will wait until I have a connection and IMD_GO orders.'"
echo "Run the connection cell in your notebook!"
echo ""
echo "================================================"
echo ""

# Run simulation with IMD streaming
gmx mdrun -v -nt 1 \
          -deffnm output/run \
          -imdwait \
          -imdport 8889 #   >& output/production.out

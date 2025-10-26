#!/bin/bash
# LAMMPS IMDv3 Streaming Production Run

echo "Starting LAMMPS production run with IMDv3 streaming..."
echo "Port: 8890"
echo ""

# Run LAMMPS
lmp -in input/input-streaming.inp -log output/run.log

echo "LAMMPS simulation complete!"

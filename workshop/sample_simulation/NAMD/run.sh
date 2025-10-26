#!/bin/bash
# NAMD IMDv3 Streaming Production Run

echo "Starting NAMD production run with IMDv3 streaming..."
echo "Port: 8888"
echo ""

# Run NAMD
namd3 +p4 input/input-streaming.inp > output/run.log 2>&1

echo "NAMD simulation complete!"

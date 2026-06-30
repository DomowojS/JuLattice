#!/bin/bash

# how to use:
# copy into console: (safes crash reports and logging data in 1 file)
# nohup bash queue_LabPC_Daniel.sh > queue_LabPC_Daniel.log 2>&1 &

LOCKFILE=/tmp/JuLattice_queue.lock
exec 9>"$LOCKFILE"
if ! flock -n 9; then
    echo "ERROR: Another queue instance is already running (lock: $LOCKFILE). Exiting."
    exit 1
fi

THREADS=62
PROJECT=/home/daniel/Software/JuLattice

echo "================================"
echo "  JuLattice Simulation Queue"
echo "================================"



echo "[1/4] Starting: Reynolds Study Re=10k D/dx=27"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS Reynolds_study_10k_Ddx27.jl
echo "[1/4] Done (exit code $?): $(date)"

echo "[2/4] Starting: Reynolds Study Re=50k D/dx=27"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS Reynolds_study_50k_Ddx27.jl
echo "[2/4] Done (exit code $?): $(date)"

echo "[3/4] Starting: Reynolds Study Re=100k D/dx=27"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS Reynolds_study_100k_Ddx27.jl
echo "[3/4] Done (exit code $?): $(date)"

echo "[4/4] Starting: Reynolds Study Re=300k D/dx=27"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS Reynolds_study_300k_Ddx27.jl
echo "[4/4] Done (exit code $?): $(date)"


echo "================================="
echo "  All simulations complete"
echo "================================="
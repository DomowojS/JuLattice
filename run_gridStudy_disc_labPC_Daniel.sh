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



echo "[1/3] Starting: Grid Study D/dx=40"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS gridStudy3_Re2760_Ma0.1_Ct0.61_Ddx=40.jl
echo "[1/3] Done (exit code $?): $(date)"

echo "[2/3] Starting: Grid Study D/dx=54"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS gridStudy3_Re2760_Ma0.1_Ct0.61_Ddx=54.jl
echo "[2/3] Done (exit code $?): $(date)"

echo "[3/3] Starting: Grid Study D/dx=70"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS gridStudy3_Re2760_Ma0.1_Ct0.61_Ddx=70.jl
echo "[3/3] Done (exit code $?): $(date)"


echo "================================="
echo "  All simulations complete"
echo "================================="
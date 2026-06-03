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

THREADS=64
PROJECT=/home/daniel/Software/JuLattice

echo "================================"
echo "  JuLattice Simulation Queue"
echo "================================"

# echo "[1/6] Starting: Grid Study D/dx=25 | VF-outflow"
# JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_gridStudy_Ddeltax25_VF-outflow.jl
# echo "[1/6] Done (exit code $?): $(date)"

echo "[2/6] Starting: Grid Study D/dx=30 | VF-outflow"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_gridStudy_Ddeltax30_VF-outflow.jl
echo "[2/6] Done (exit code $?): $(date)"

echo "[3/6] Starting: Ma Study D/dx=40 Ma=0.05 | VF-outflow"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_MaStudy_Ddeltax40_VF-outflow.jl
echo "[3/6] Done (exit code $?): $(date)"

echo "[4/6] Starting: Grid Study D/dx=25 | extrapolation-outflow"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_gridStudy_Ddeltax25_extrapolation-outflow.jl
echo "[4/6] Done (exit code $?): $(date)"

echo "[5/6] Starting: Grid Study D/dx=30 | extrapolation-outflow"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_gridStudy_Ddeltax30_extrapolation-outflow.jl
echo "[5/6] Done (exit code $?): $(date)"

echo "[6/6] Starting: Ma Study D/dx=40 Ma=0.05 | extrapolation-outflow"
JULIA_EXCLUSIVE=1 stdbuf -oL julia --project=$PROJECT -t $THREADS run_MaStudy_Ddeltax40_extrapolation-outflow.jl
echo "[6/6] Done (exit code $?): $(date)"

echo "================================="
echo "  All simulations complete"
echo "================================="

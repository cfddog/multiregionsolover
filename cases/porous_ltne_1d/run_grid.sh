#!/bin/bash
# Mesh-convergence study for porous_ltne_1d.
# usage: bash run_grid.sh <nx> <max_iter>
set -e
cd "$(dirname "$0")"
NX=$1
ITER=$2
echo "=============== nx=$NX  max_iter=$ITER ==============="
python3 gen_case.py "$NX" 2 2
# patch control.ec iteration budget
sed -i "s/^  Porous_Max_Iter=.*/  Porous_Max_Iter=$ITER/" control.ec
grep -n "Porous_Max_Iter=" control.ec
mpirun -np 1 ./opencfd-ec1.16a.out > run_grid_${NX}.log 2>&1 || { echo "RUN FAILED nx=$NX"; tail -30 run_grid_${NX}.log; exit 1; }
python3 validate_ltne.py | tee validate_grid_${NX}.txt

#!/bin/bash
# 极限测试批量运行脚本:  bash run_limits.sh <dir1> [dir2 ...]
# 仓库根目录由脚本位置推导（cases/fluid_solid/limits -> ../../..），可整体拷贝到别处运行
L=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$L/../../.." && pwd)
EXE=$R/src/MultiPhysMixSolver
if [ ! -x "$EXE" ]; then
  echo "executable not found: $EXE  (run 'make' in $R/src first)"; exit 1
fi
for d in "$@"; do
  cd $L/$d || exit 1
  rm -f flow3d_node.dat iface_couple.dat iface_highlow.dat Ts_block_1.dat \
        output_para.out flow3d.dat flow3d.vtk flow3d_*.vtk field_restart.dat \
        part_grid.dat mesh-quality.dat Step_mess.dat bc3d.inc run.log
  echo "=== START $d $(date +%H:%M:%S) ===" >> $L/status.txt
  ( date +%H:%M:%S; mpirun -np 1 $EXE; echo "EXIT=$?"; date +%H:%M:%S ) > run.log 2>&1
  echo "=== END   $d $(date +%H:%M:%S) ===" >> $L/status.txt
done
echo "ALLDONE $(date +%H:%M:%S)" >> $L/status.txt

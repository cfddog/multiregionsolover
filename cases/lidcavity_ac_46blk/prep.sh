#!/bin/bash
# 建立 4/6 块运行目录（输入用完软链），供任务 7（强扩展）与任务 8（anchoring）用。
#
#   m4a_np1 m4a_np2 m4a_np4   : 2x2 块，行主序编号，np=1/2/4
#   m4b_np1                   : 2x2 块，逆序编号（全局 1 号块在对角）
#   m6a_np1 m6a_np2 m6a_np3 m6a_np6 : 3x2 块，行主序编号
#   m6b_np1                   : 3x2 块，逆序编号
#
# 用法：./prep.sh [总步数]
#   不给参数：各目录软链算例自带 control.ec（t_end=10000）
#   给步数（冒烟用 1500）：各目录写独立 control.ec，
#   t_end=Kstep_save=步数、Kstep_show=步数/5
#
# 环境变量 FRESH=1（冒烟/对比必用）：先删掉目录里上次的 field_restart.dat
# 与输出（flow3d*.dat / flow3d_block_*.vtk / Residual.dat / Step_mess.dat）。
# 否则重跑会 Iflag_init=0 从 field_restart.dat 读回上次的 Kstep，若该文件已经
# 到 t_end，程序读入后立刻结束 -> 得到的还是旧结果（2026-10-18 踩过：m2b/m6b/
# m6c 的"修复后重跑"其实什么都没算，VTK 时间戳停在上一次）。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
NSTEP=${1:-}

mk() {   # mk <rundir> <mesh> <bc> <material> "<proc of each block | ''>"
  local d=$1 mesh=$2 bc=$3 mat=$4 procs=$5
  mkdir -p "$d"
  ln -sfn "../$mesh" "$d/Mesh3d.x"
  ln -sfn "../$bc"   "$d/bc3d.inp"
  ln -sfn "../$mat"  "$d/material.in"
  rm -f "$d/control.ec"
  if [ -n "$NSTEP" ]; then
    sed -e "s/^  t_end=.*/  t_end=${NSTEP}.d0/" \
        -e "s/^  Kstep_save=.*/  Kstep_save=${NSTEP}/" \
        -e "s/^  Kstep_show=.*/  Kstep_show=$((NSTEP / 5 > 0 ? NSTEP / 5 : 1))/" \
        control.ec > "$d/control.ec.tmp"
    mv "$d/control.ec.tmp" "$d/control.ec"
  else
    ln -sfn ../control.ec "$d/control.ec"
  fi
  rm -f "$d/partation.dat"
  if [ -n "$FRESH" ]; then
    rm -f "$d/field_restart.dat" "$d/flow3d.dat" "$d/Residual.dat" "$d/Step_mess.dat"
    rm -f "$d"/flow3d_block_*.vtk
  fi
  [ -z "$procs" ] && return 0
  local nproc=0 b=0 p
  for p in $procs; do [ "$p" -gt "$nproc" ] && nproc=$p; done
  nproc=$((nproc + 1))
  {
    printf '  Total_block_number    Total_Proc_number\n'
    printf '%12d%12d\n' "$(echo $procs | wc -w)" "$nproc"
    printf '  Block_no, Proc_no\n'
    for p in $procs; do
      b=$((b + 1))
      printf '%12d%12d\n' "$b" "$p"
    done
  } > "$d/partation.dat"
}

mk m2a_np1 Mesh3d_m2a.x bc3d_m2a.inp material_m2.in ""
mk m2b_np1 Mesh3d_m2b.x bc3d_m2b.inp material_m2.in ""
mk m4a_np1 Mesh3d_m4a.x bc3d_m4a.inp material_m4.in ""
mk m4a_np2 Mesh3d_m4a.x bc3d_m4a.inp material_m4.in "0 0 1 1"
mk m4a_np4 Mesh3d_m4a.x bc3d_m4a.inp material_m4.in "0 1 2 3"
mk m4b_np1 Mesh3d_m4b.x bc3d_m4b.inp material_m4.in ""
mk m6a_np1 Mesh3d_m6a.x bc3d_m6a.inp material_m6.in ""
mk m6a_np2 Mesh3d_m6a.x bc3d_m6a.inp material_m6.in "0 0 0 1 1 1"
mk m6a_np3 Mesh3d_m6a.x bc3d_m6a.inp material_m6.in "0 0 1 1 2 2"
mk m6a_np6 Mesh3d_m6a.x bc3d_m6a.inp material_m6.in "0 1 2 3 4 5"
mk m6b_np1 Mesh3d_m6b.x bc3d_m6b.inp material_m6.in ""
mk m6c_np1 Mesh3d_m6c.x bc3d_m6c.inp material_m6.in ""
mk m6d_np1 Mesh3d_m6d.x bc3d_m6d.inp material_m6.in ""
echo "4/6-block run dirs ready（步数 ${NSTEP:-算例自带 control.ec}）"

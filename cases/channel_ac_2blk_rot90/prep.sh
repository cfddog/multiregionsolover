#!/bin/bash
# 建立 ref_np1/ p_np1/ p_np2/ n_np1/ n_np2/ 五个运行目录（输入用软链）。
#
#   ref_np1 : 非旋转的两块基准（直接引用 ../channel_ac_2blk 的网格与 bc3d），
#             用来给本算例提供"同一套 AC 设置下的逐点参照"
#   p_np1/p_np2 : i->+y、k->-z，右手 L=(2,+1,-3)，np=1 / np=2
#   n_np1/n_np2 : i->-y、k->+z，右手 L=(2,-1,+3)（走反向 i 的负号路径），np=1 / np=2
#
# 用法：./prep.sh [总步数]
#   不给参数：各目录软链算例自带 control.ec（t_end=40000）
#   给步数（冒烟用 2000）：各目录写独立 control.ec，
#   t_end=Kstep_save=步数、Kstep_show=步数/5
#
# 环境变量 FRESH=1（重跑/对比必用）：先删掉目录里上次的 field_restart.dat 与输出
# （flow3d*.dat / flow3d_block_*.vtk / Residual.dat / Step_mess.dat）。否则
# Iflag_init=0 会从 field_restart.dat 读回上次的 Kstep，文件已到 t_end 时程序
# 读入后立刻结束，得到的还是旧结果（2026-10-18 在 46 块算例踩过）。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
NSTEP=${1:-}

mk() {   # mk <rundir> <mesh> <bc> <np>
  local d=$1 mesh=$2 bc=$3 np=$4
  mkdir -p "$d"
  ln -sfn "$mesh" "$d/Mesh3d.x"
  ln -sfn "$bc"   "$d/bc3d.inp"
  ln -sfn ../material.in "$d/material.in"
  # 先删旧的（可能是软链）：直接 `> 软链` 会截断目标文件（2026-10-18 踩过）
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
  if [ "$np" -gt 1 ]; then
    { printf '  Total_block_number    Total_Proc_number\n'
      printf '%12d%12d\n' 2 "$np"
      printf '  Block_no, Proc_no\n'
      printf '%12d%12d\n' 1 0
      printf '%12d%12d\n' 2 1
    } > "$d/partation.dat"
  fi
}

mk ref_np1 ../../channel_ac_2blk/Mesh3d.x ../../channel_ac_2blk/bc3d.inp   1
mk p_np1   ../Mesh3d_p.x             ../bc3d_p.inp                    1
mk p_np2   ../Mesh3d_p.x             ../bc3d_p.inp                    2
mk n_np1   ../Mesh3d_n.x             ../bc3d_n.inp                    1
mk n_np2   ../Mesh3d_n.x             ../bc3d_n.inp                    2
echo "ref_np1/ p_np1/ p_np2/ n_np1/ n_np2/ ready（步数 ${NSTEP:-算例自带 control.ec}）"

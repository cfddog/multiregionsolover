#!/bin/bash
# 本地 <-> 集群（wuzh02）同步。规则见 .clinerules / memory-bank/constraints.md §8。
#
#   bash util/sync_cluster.sh push   # 本地 -> 集群：代码 / 输入 / 文档 / 记忆库
#   bash util/sync_cluster.sh pull   # 集群 -> 本地：只取必要记录（日志、分析输出）
#
# 可选环境变量：HOST=wuzh02  REMOTE=sundong/MultiPhysMixSolver  DRY=1（push 只演练）
set -euo pipefail
cd "$(dirname "$0")/.."

HOST=${HOST:-wuzh02}
REMOTE=${REMOTE:-sundong/MultiPhysMixSolver}
SSH_CMD="ssh -o BatchMode=yes -o ConnectTimeout=15"

# 不同步：版本库 / 编译产物 / 可视化与重启大文件 / 实测运行目录
# 注意 Mesh3d.x 是算例输入，必须同步；bc3d.inc 有历史算例入库，也不排除。
EXCLUDES=(
  --exclude '.git/' --exclude 'build/' --exclude 'bin/' --exclude '__pycache__/'
  --exclude '*.vtk' --exclude '*.o' --exclude '*.mod' --exclude '*.pyc'
  --exclude 'flow3d.dat' --exclude 'flow3d_node.dat' --exclude 'field_restart.dat'
  --exclude 'Mesh3d.dat' --exclude 'partation.dat' --exclude 'Residual.dat'
  --exclude '*.log' --exclude '*.pw' --exclude '*.plt' --exclude '*.tec'
  --exclude 'np[0-9]*/'
  # 并行验证框架的运行目录：两侧各自生成，push 绝不能拿本地旧结果覆盖集群的新结果
  # （2026-10-18 踩过：push 用 -a 保留 mtime，把集群刚算出的 analyze_np*.txt 盖回旧版）
  --exclude 'runs/'
)

case "${1:-}" in
push)
  if [ "${DRY:-0}" = "1" ]; then
    rsync -azn --stats "${EXCLUDES[@]}" -e "$SSH_CMD" ./ "$HOST:$REMOTE/" | tail -12
    echo "[dry-run] 未做任何改动"
    exit 0
  fi
  rsync -az --delete "${EXCLUDES[@]}" -e "$SSH_CMD" ./ "$HOST:$REMOTE/"
  echo "push 完成：$(pwd) -> $HOST:$REMOTE"
  ;;
pull)
  # --update：本地更新（未推送）的文件不被集群旧版覆盖
  rsync -az --update -e "$SSH_CMD" "$HOST:$REMOTE/memory-bank/" memory-bank/
  rsync -az --update -e "$SSH_CMD" "$HOST:$REMOTE/docs/" docs/
  # cases 下只回传日志/分析文本/现场控制参数，VTK、重启文件与巨大的 run_out.log 留在集群
  rsync -az --update -e "$SSH_CMD" \
    --include '*/' --exclude 'run_out.log' \
    --include 'run_err.log' --include '*.log' --include 'Residual.dat' \
    --include 'control.ec' --include 'output_para.out' \
    --include 'analyze_*.txt' --exclude '*' \
    "$HOST:$REMOTE/cases/" cases/
  echo "pull 完成：$HOST:$REMOTE -> $(pwd)（仅 memory-bank/docs/日志）"
  ;;
*)
  sed -n '2,6p' "$0"
  exit 1
  ;;
esac

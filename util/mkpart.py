#!/usr/bin/env python3
"""mkpart.py — 生成 partation.dat（全局块号 -> 进程号）—— 并行验证/算例通用工具。

程序读法见 src/sub_partation_mpi.f90::partation（Fortran list-directed，空白分隔）：

    Total_block_number    Total_Proc_number
            4            2
    Block_no, Proc_no
            1            0
            2            0
            3            1
            4            1

第 1 行是表头文字（读时跳过），第 2 行 "块数 进程数"，第 3 行也是表头，
其后每行 "块号 进程号"（**进程号从 0 起**，块号从 1 起，行序即全局块序）。

用法
----
  # 均衡连续切分：4 块分给 2 个进程 -> {1,2} | {3,4}
  python3 util/mkpart.py --nblocks 4 --np 2 --split 2,2

  # 不均衡（专测负载不均）：6 块 -> 3 / 2 / 1
  python3 util/mkpart.py --nblocks 6 --np 3 --split 3,2,1

  # 直接给每个块的进程号（1 起编号，脚本内部转 0 起）
  python3 util/mkpart.py --blocks 1,1,2

  # 块数直接从 Mesh3d.x 的第 1 条记录读（Fortran unformatted）
  python3 util/mkpart.py --mesh cases/channel_ac_2blk/Mesh3d.x --np 2

  # 写到指定位置（默认当前目录 partation.dat）
  python3 util/mkpart.py --nblocks 4 --np 4 -o runs/t02/np4/partation.dat

注意
----
* np=1 **不要**放 partation.dat：程序找不到文件时会自动分区（全部块在 0 号进程），
  而进程数 1 与文件里的 Total_proc 必须一致，否则 partation() 直接 stop。
* 每个进程至少要分到 1 个块（Num_block=0 的进程没有 my_blocks，后续会在别处出错），
  所以 split 里不允许出现 0；确实需要空进程请改代码后再说。
"""
import argparse
import os
import struct
import sys


def nblocks_from_mesh(path):
    """Mesh3d.x 第 1 条记录 = Total_block（Fortran unformatted，4 字节整数）。"""
    with open(path, "rb") as f:
        raw = f.read(4)
        if len(raw) < 4:
            raise SystemExit("mkpart: %s 太短，读不出块数" % path)
        rc = struct.unpack("<i", raw)[0]
        if rc <= 0 or rc % 4:
            raise SystemExit("mkpart: %s 不是 Fortran unformatted 文件（记录长度 %d）" % (path, rc))
        payload = f.read(rc)
    n = struct.unpack("<%di" % (rc // 4), payload)[0]
    if n <= 0:
        raise SystemExit("mkpart: %s 第 1 条记录 = %d，不是合法块数" % (path, n))
    return n


def split_to_blocks(split, nblocks):
    """连续切分 -> 每块（1 起）的进程号（0 起）。"""
    if len(split) < 1:
        raise SystemExit("mkpart: --split 不能为空")
    if sum(split) != nblocks:
        raise SystemExit("mkpart: --split 之和 %d != 块数 %d" % (sum(split), nblocks))
    if any(k < 1 for k in split):
        raise SystemExit("mkpart: --split 里出现 0（每个进程至少要 1 个块），split=%s"
                         % ",".join(str(k) for k in split))
    out = []
    for p, k in enumerate(split):
        out.extend([p] * k)
    return out


def blocks_to_split(proc_of_block):
    """每块进程号 -> (split 计数, np)。"""
    np_ = max(proc_of_block) + 1
    split = [0] * np_
    for p in proc_of_block:
        split[p] += 1
    return split, np_


def write_partation(path, nblocks, nproc, proc_of_block):
    d = os.path.dirname(os.path.abspath(path))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, "w") as f:
        f.write(" Total_block_number    Total_Proc_number \n")
        f.write("%12d%12d\n" % (nblocks, nproc))
        f.write(" Block_no, Proc_no \n")
        for m, p in enumerate(proc_of_block, start=1):
            f.write("%12d%12d\n" % (m, p))


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="mkpart.py",
        description="生成 partation.dat（全局块号 -> 进程号，进程号 0 起）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("用法\n----\n", 1)[-1])
    ap.add_argument("--nblocks", type=int, help="全局块数")
    ap.add_argument("--mesh", help="改从 Mesh3d.x 读块数（Fortran unformatted 第 1 条记录）")
    ap.add_argument("--np", type=int, help="进程数（--split 时可省，取 len(split)）")
    ap.add_argument("--split", help="每个进程的连续块数，如 '2,2' / '3,2,1'")
    ap.add_argument("--blocks", help="每个块（1 起）的进程号（1 起），如 '1,1,2'")
    ap.add_argument("-o", "--out", default="partation.dat", help="输出文件（默认 partation.dat）")
    ap.add_argument("-q", "--quiet", action="store_true", help="不打印摘要")
    a = ap.parse_args(argv)

    if a.split and a.blocks:
        raise SystemExit("mkpart: --split 与 --blocks 只能给一个")

    nblocks = a.nblocks
    if a.mesh:
        n_mesh = nblocks_from_mesh(a.mesh)
        if nblocks is not None and nblocks != n_mesh:
            raise SystemExit("mkpart: --nblocks %d 与 %s 的块数 %d 不一致"
                             % (nblocks, a.mesh, n_mesh))
        nblocks = n_mesh
    if nblocks is None and a.split:
        nblocks = sum(int(x) for x in a.split.replace(" ", "").split(",") if x)
    if nblocks is None and a.blocks:
        nblocks = len([x for x in a.blocks.replace(" ", "").split(",") if x])
    if nblocks is None:
        raise SystemExit("mkpart: 需要 --nblocks（或 --mesh / --split / --blocks 推出块数）")

    if a.blocks:
        lst = [int(x) for x in a.blocks.replace(" ", "").split(",") if x]
        if len(lst) != nblocks:
            raise SystemExit("mkpart: --blocks 给了 %d 个值，块数是 %d" % (len(lst), nblocks))
        if any(v < 1 for v in lst):
            raise SystemExit("mkpart: --blocks 的进程号从 1 起（内部转 0 起），收到 %s" % a.blocks)
        proc_of_block = [v - 1 for v in lst]
        split, np_file = blocks_to_split(proc_of_block)
    elif a.split:
        split = [int(x) for x in a.split.replace(" ", "").split(",") if x]
        if a.np is not None and a.np != len(split):
            raise SystemExit("mkpart: --np %d 与 --split 的项数 %d 不一致"
                             % (a.np, len(split)))
        np_file = len(split)
        proc_of_block = split_to_blocks(split, nblocks)
    else:
        if a.np is None:
            raise SystemExit("mkpart: 需要 --np（或给 --split / --blocks）")
        if a.np < 1:
            raise SystemExit("mkpart: --np 至少为 1")
        if a.np > nblocks:
            raise SystemExit("mkpart: --np %d > 块数 %d（会有进程分不到块）" % (a.np, nblocks))
        # 均衡连续切分：前 rem 个进程多 1 块
        base, rem = divmod(nblocks, a.np)
        split = [base + (1 if p < rem else 0) for p in range(a.np)]
        proc_of_block = split_to_blocks(split, nblocks)
        np_file = a.np

    write_partation(a.out, nblocks, np_file, proc_of_block)
    if not a.quiet:
        print("mkpart: %s  (blocks=%d, np=%d)" % (a.out, nblocks, np_file))
        print("        split   = %s" % ",".join(str(k) for k in
                                               blocks_to_split(proc_of_block)[0]))
        print("        blocks  = %s" % ",".join(str(p + 1) for p in proc_of_block))
    return 0


if __name__ == "__main__":
    sys.exit(main())

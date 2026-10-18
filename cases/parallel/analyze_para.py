#!/usr/bin/env python3
"""analyze_para.py — 并行一致性断言：把测试目录与 np=1 参考目录逐点比较。

用法
----
  python3 analyze_para.py <ref_dir(np=1)> <test_dir(np>1)> [选项]

选项
----
  --tol-vel 1e-5   速度 max|Δ| 容差
  --tol-rms 1e-6   速度 RMS|Δ| 容差
  --tol-p   1e-6   压力**去常数后**的相对容差（max|Δp-mean(Δp)| / max|p|）
  --tol-T   1e-6   温度 T/Ts 的相对容差（max|ΔT| / max|T|）
  --no-vtk         不比较 vtk 压力（只比 flow3d_node.dat）
  --quiet          只打印结论
退出码：0 = PASS，1 = FAIL，2 = 输入缺失/格式错。

比较的两个文件
--------------
1. `flow3d_node.dat`（`$flow_ec: Iflag_flow_node=1`）——由
   `src/sub_IO.f90::output_flow_node` 写出：rank 0 按**全局块序**聚合，文件布局
   与进程数无关，可直接逐点相减。变量顺序 d,u,v,w,T,Ts（SI 单位）。
2. `flow3d_block_<m>.vtk` 的 `pressure`（单元中心、ASCII）——flow3d_node.dat
   里没有压力，而 AC 的 p 是独立数组。封闭域的压力水平只是 gauge
   （`init_ac_parallel_flags` 用全局 1 号块锚定），所以判据用**去常数**后的
   最大差，同时报告原始 max|Δp| 与常数偏移供人工判断。
"""
import argparse
import glob
import os
import re
import struct
import sys

# 集群上 LANG=C 时 stdout 的编码是 ascii，打印中文会抛 UnicodeEncodeError
# （2026-10-18 在 wuzh02 上第一次跑就踩到，导致 4 个比对结果全部变成 FAIL）
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

NVAR = 6
VARNAMES = ["d", "u", "v", "w", "T", "Ts"]


# ---------------------------------------------------------------- flow3d_node
def _records(path):
    """Fortran unformatted 文件的逐记录 payload。"""
    recs = []
    with open(path, "rb") as f:
        while True:
            raw = f.read(4)
            if len(raw) < 4:
                break
            rc = struct.unpack("<I", raw)[0]
            recs.append(f.read(rc))
            f.read(4)
    return recs


def read_flow_node(d):
    path = os.path.join(d, "flow3d_node.dat")
    if not os.path.isfile(path):
        raise SystemExit("analyze_para: 缺少 %s（control.ec 里要 Iflag_flow_node=1）" % path)
    recs = _records(path)
    if len(recs) < 4:
        raise SystemExit("analyze_para: %s 记录数 %d 太少（需带头）" % (path, len(recs)))
    nb = struct.unpack("<i", recs[0])[0]
    hdr = struct.unpack("<%di" % (len(recs[1]) // 4), recs[1])
    if len(hdr) == 3 * nb and len(recs[2]) == 4:
        dims = [tuple(hdr[3 * k:3 * k + 3]) for k in range(nb)]
        nvar = struct.unpack("<i", recs[2])[0]
        data = recs[3:]
    elif len(hdr) == 4 * nb:            # 旧变体：记录2 = 逐块 (ni,nj,nk,nvar)
        dims = [(hdr[4 * k], hdr[4 * k + 1], hdr[4 * k + 2]) for k in range(nb)]
        nvar = hdr[3]
        data = recs[2:]
    else:
        raise SystemExit("analyze_para: %s 文件头无法识别（记录2 = %d 字节）"
                         % (path, len(recs[1])))
    if nvar != NVAR:
        raise SystemExit("analyze_para: %s 变量数 %d，本脚本需要 %d" % (path, nvar, NVAR))
    fields = {}
    for m, raw in enumerate(data, start=1):
        ni, nj, nk = dims[m - 1]
        n = ni * nj * nk
        vals = struct.unpack("<%dd" % (len(raw) // 8), raw)
        if len(vals) != NVAR * n:
            raise SystemExit("analyze_para: 块 %d 数据长度 %d != %d*%d"
                             % (m, len(vals), NVAR, n))
        for v in range(NVAR):
            fields[(m, VARNAMES[v])] = vals[v * n:(v + 1) * n]
    return dims, fields


# ------------------------------------------------------------------ vtk 压力
def _num(tok):
    """Fortran 可能打印缺 'E' 的科学计数（如 '1.234-108'）。"""
    try:
        return float(tok)
    except ValueError:
        for k in range(1, len(tok) - 1):
            if tok[k] in "+-" and tok[k - 1] not in "eEdD" and tok[k - 1].isdigit():
                return float(tok[:k] + "E" + tok[k:])
        raise


def read_vtk_pressure(path):
    """返回 pressure 值列表（单元中心，i 最快）；找不到返回 None。"""
    with open(path, "r", errors="replace") as f:
        lines = f.read().splitlines()
    for n, ln in enumerate(lines):
        sp = ln.split()
        if len(sp) > 1 and sp[0].upper() == "SCALARS" and sp[1].lower() == "pressure":
            vals = []
            m = n + 1
            while m < len(lines) and lines[m].strip().lower().startswith("lookup_table"):
                m += 1
            while m < len(lines):
                s = lines[m].strip()
                if not s:
                    m += 1
                    continue
                if s.split()[0].upper() in ("SCALARS", "VECTORS", "FIELD", "CELL_DATA",
                                            "POINT_DATA", "TENSORS", "LOOKUP_TABLE"):
                    break
                try:
                    vals.extend(_num(t) for t in s.split())
                except ValueError:
                    break
                m += 1
            return vals
    return None


# ------------------------------------------------------------------- 指标
def stat(ref, tst):
    """(max|Δ|, rms|Δ|, max|ref|, argmax 下标)。"""
    dmax = rmax = 0.0
    s2 = 0.0
    arg = -1
    for n in range(len(ref)):
        dd = ref[n] - tst[n]
        ad = abs(dd)
        if ad > dmax:
            dmax, arg = ad, n
        s2 += dd * dd
        if abs(ref[n]) > rmax:
            rmax = abs(ref[n])
    rms = (s2 / len(ref)) ** 0.5 if ref else 0.0
    return dmax, rms, rmax, arg


# --------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(prog="analyze_para.py",
                                 description="np=1 参考 <-> np>1 测试的逐点一致性断言")
    ap.add_argument("ref")
    ap.add_argument("test")
    ap.add_argument("--tol-vel", type=float, default=1e-5)
    ap.add_argument("--tol-rms", type=float, default=1e-6)
    ap.add_argument("--tol-p", type=float, default=1e-6)
    ap.add_argument("--tol-T", type=float, default=1e-6,
                    help="温度 T/Ts 的**相对**容差 max|dT|/max|T_ref|")
    ap.add_argument("--no-vtk", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    refd, tstd = a.ref.rstrip("/"), a.test.rstrip("/")

    dimsR, fR = read_flow_node(refd)
    dimsT, fT = read_flow_node(tstd)
    if dimsR != dimsT:
        raise SystemExit("analyze_para: 网格维数不同\n  ref =%s\n  test=%s" % (dimsR, dimsT))

    if not a.quiet:
        print("ref : %s" % refd)
        print("test: %s" % tstd)
        print("=" * 92)
        print("%-5s %-3s | %-12s %-12s %-12s | %-9s %s"
              % ("block", "var", "max|d|", "rms|d|", "max|ref|", "rel", "argmax"))
        print("-" * 92)

    worst = {"vel": 0.0, "rms": 0.0, "d": 0.0, "Trel": 0.0}
    fail = []
    for m in range(1, len(dimsR) + 1):
        for v in VARNAMES:
            r, t = fR[(m, v)], fT[(m, v)]
            d1, r1, big, arg = stat(r, t)
            rel = d1 / big if big > 0 else d1
            if not a.quiet:
                print("%-5d %-3s | %-12.4g %-12.4g %-12.4g | %-9.2g %d"
                      % (m, v, d1, r1, big, rel, arg + 1))
            if v in ("u", "v", "w"):
                worst["vel"] = max(worst["vel"], d1)
                worst["rms"] = max(worst["rms"], r1)
            elif v == "d":
                worst["d"] = max(worst["d"], d1)
            elif v in ("T", "Ts"):
                worst["Trel"] = max(worst["Trel"], rel)
    if worst["vel"] > a.tol_vel or worst["rms"] > a.tol_rms:
        fail.append("速度 max|d|=%.3g (tol %.3g) / rms=%.3g (tol %.3g)"
                    % (worst["vel"], a.tol_vel, worst["rms"], a.tol_rms))
    if worst["d"] > a.tol_vel:
        fail.append("密度 max|d|=%.3g (tol %.3g)" % (worst["d"], a.tol_vel))
    if worst["Trel"] > a.tol_T:
        fail.append("温度 相对 max|dT|/max|T| = %.3g > tol %.3g"
                    % (worst["Trel"], a.tol_T))

    # ---- 压力（vtk 单元中心）：去常数后比较
    pmsg = "压力场：跳过（--no-vtk）"
    if not a.no_vtk:
        files = sorted(glob.glob(os.path.join(refd, "flow3d_block_*.vtk")),
                       key=lambda s: int(re.sub(r"\D", "", os.path.basename(s)) or 0))
        if not files:
            pmsg = "压力场：跳过（没有 flow3d_block_*.vtk）"
        else:
            praw = pdev = poff = pmax = 0.0
            ncnt = 0
            for fr in files:
                bn = os.path.basename(fr)
                ft = os.path.join(tstd, bn)
                if not os.path.isfile(ft):
                    fail.append("缺少 %s" % ft)
                    continue
                pr = read_vtk_pressure(fr)
                pt = read_vtk_pressure(ft)
                if pr is None or pt is None:
                    continue
                if len(pr) != len(pt):
                    fail.append("%s 压力长度 %d != %d" % (bn, len(pr), len(pt)))
                    continue
                dmean = sum(y - x for x, y in zip(pr, pt)) / len(pr)
                for x, y in zip(pr, pt):
                    praw = max(praw, abs(y - x))
                    pdev = max(pdev, abs((y - x) - dmean))
                    pmax = max(pmax, abs(x))
                poff = max(poff, abs(dmean))
                ncnt += len(pr)
            if ncnt:
                relp = pdev / pmax if pmax > 0 else pdev
                pmsg = ("压力场：%d 单元  max|dp|=%.3g  去常数 max|dp-mean|=%.3g "
                        "(rel %.2g)  常数偏移 max|mean_dp|=%.3g"
                        % (ncnt, praw, pdev, relp, poff))
                if relp > a.tol_p:
                    fail.append("压力 去常数 max|dp-mean|/max|p| = %.3g > tol %.3g"
                                % (relp, a.tol_p))

    print("-" * 92)
    print("最坏值：速度 max|d|=%.3g  rms=%.3g   密度 max|d|=%.3g   温度 rel=%.3g"
          % (worst["vel"], worst["rms"], worst["d"], worst["Trel"]))
    print(pmsg)
    if fail:
        print("结论：FAIL")
        for s in fail:
            print("  - " + s)
        return 1
    print("结论：PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

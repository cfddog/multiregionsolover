#!/usr/bin/env python3
"""check_flow3d_node.py — 校验 flow3d_node.dat（节点值、SI 单位）与 Mesh3d.x 是否配套。

用法：python3 check_flow3d_node.py [算例目录]        （默认当前目录）

flow3d_node.dat 由 src/sub_IO.f90::output_flow_node 写出（control.ec 的
$flow_ec: Iflag_flow_node=1）：unformatted、块序与 Mesh3d.x 相同、
每条记录含 6 个节点场： d, u, v, w, T, Ts
  T  = 流体温度 [K]（固体块为固体温度 Ts，兼容 flow3d.vtk 约定）
  Ts = 固体/骨架温度 [K]（固体与多孔为真实值；流体/低速块无固相，Ts=T）

文件布局（2026-09-25 起带 PLOT3D 解(Q)/函数型文件头）：
  记录1 = Total_block                     块数
  记录2 = (ni,nj,nk) m=1..Total_block     各块节点维数
  记录3 = 变量数 (= 6)
  其后逐块一条数据记录： 6*ni*nj*nk 个 double（i 最快、j 次、k 最慢）
旧文件（无头）仍可识别：此时维数由 Mesh3d.x 提供。
"""
import struct
import sys
import os

d = (sys.argv[1] if len(sys.argv) > 1 else ".").rstrip("/")
mesh = os.path.join(d, "Mesh3d.x")
nodef = os.path.join(d, "flow3d_node.dat")
NVAR = 6


def read_plot3d_grid_dims(path):
    """返回 (nblock, [(ni,nj,nk), ...])（只读前两条记录）。"""
    with open(path, "rb") as f:
        rc = struct.unpack("<I", f.read(4))[0]
        nb = struct.unpack("<i", f.read(rc))[0]
        f.read(4)
        rc = struct.unpack("<I", f.read(4))[0]
        dims = struct.unpack("<%di" % (rc // 4), f.read(rc))
        f.read(4)
    return nb, [tuple(dims[3 * k:3 * k + 3]) for k in range(nb)]


def read_records(path):
    """返回每条 Fortran unformatted 记录的 payload 字节串列表。"""
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


nb, dims = read_plot3d_grid_dims(mesh)
print("Mesh3d.x    : nblock=%d  dims=%s" % (nb, dims))

recs = read_records(nodef)

# ---- 头记录识别：新格式（带 PLOT3D Q 头） vs 旧格式（无头） ----
nvars = NVAR
hdr_dims = None
hdr_nb = None
data = recs
fmt = "legacy (no header)"
if len(recs) >= 3 and len(recs[0]) == 4 and struct.unpack("<i", recs[0])[0] == nb:
    hdr = struct.unpack("<%di" % (len(recs[1]) // 4), recs[1])
    if len(hdr) == 4 * nb:
        # 过渡布局（2026-09-26 曾在工作树中出现）：记录2 = 逐块 (ni,nj,nk,nvar)，
        # 没有独立的记录3。保留兼容，便于校验当时写出的旧文件。
        hdr_nb = struct.unpack("<i", recs[0])[0]
        hdr_dims = [hdr[4 * k + c] for k in range(nb) for c in range(3)]
        nvars = hdr[3]
        data = recs[2:]
        fmt = ("PLOT3D-Q header, 2-record variant (nblock=%d, dims+nvars=%d)"
               % (hdr_nb, nvars))
    elif len(hdr) == 3 * nb and len(recs[2]) == 4:
        # 标准写法：记录1 块数 / 记录2 (ni,nj,nk) 逐块 / 记录3 变量数
        hdr_nb = struct.unpack("<i", recs[0])[0]
        hdr_dims = hdr
        nvars = struct.unpack("<i", recs[2])[0]
        data = recs[3:]
        fmt = "PLOT3D-Q header (nblock=%d, dims, nvars=%d)" % (hdr_nb, nvars)

print("flow3d_node : %d records, %d vars/record, format=%s"
      % (len(data), nvars, fmt))

ok = True
if hdr_nb is not None:
    if len(hdr_dims) != 3 * nb:
        ok = False
        print("  HEADER: dims record has %d ints, expect %d" % (len(hdr_dims), 3 * nb))
    else:
        hd = [tuple(hdr_dims[3 * k:3 * k + 3]) for k in range(nb)]
        if hd != dims:
            ok = False
            print("  HEADER: dims %s != Mesh3d.x %s" % (hd, dims))
    if nvars != NVAR:
        ok = False
        print("  HEADER: nvars=%d, this checker expects %d" % (nvars, NVAR))

if len(data) != nb:
    ok = False
    print("  RECORDS: got %d, Mesh3d.x has %d blocks" % (len(data), nb))

for m, raw in enumerate(data):
    if m >= len(dims):
        break
    ni, nj, nk = dims[m]
    n = ni * nj * nk
    vals = struct.unpack("<%dd" % (len(raw) // 8), raw)
    if len(vals) != nvars * n:
        ok = False
        print("  block %d: SIZE MISMATCH got %d expect %d" % (m + 1, len(vals), nvars * n))
        continue
    if nvars < 6:
        print("  block %d: npts=%-6d (=%dx%dx%d)  %d vars (no d/T/Ts report)"
              % (m + 1, n, ni, nj, nk, nvars))
        continue
    dd = vals[0 * n:1 * n]
    uu = vals[1 * n:2 * n]
    tt = vals[4 * n:5 * n]
    ts = vals[5 * n:6 * n]
    print("  block %d: npts=%-6d (=%dx%dx%d)  d[%.4g,%.4g] |u|[%.4g,%.4g]  "
          "T[%.4g,%.4g] K  Ts[%.4g,%.4g] K  max|T-Ts|=%.4g K"
          % (m + 1, n, ni, nj, nk, min(dd), max(dd),
             min(min(uu), min(vals[2 * n:3 * n]), min(vals[3 * n:4 * n])),
             max(max(uu), max(vals[2 * n:3 * n]), max(vals[3 * n:4 * n])),
             min(tt), max(tt), min(ts), max(ts),
             max(abs(a - b) for a, b in zip(tt, ts))))
print("STRUCTURE OK" if ok else "STRUCTURE MISMATCH")

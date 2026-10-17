#!/usr/bin/env python3
"""compare_flow3d_node.py — 逐点比较两个算例的 flow3d_node.dat（节点值、SI 单位）。

用法：
  python3 compare_flow3d_node.py <ref_dir> <test_dir> [选项]

选项：
  --vars   d,u,v,w,T,Ts       要比较的变量（默认全部 6 个）
  --blocks 1,2                要比较的块号（默认全部）
  --iface  1:j1,2:i101        各块界面节点行的显式位置（默认本工程 fluid_solid
                              拓扑：块1 界面在 j=1，块2 界面在 i=nx）
  --csv    out.csv            把逐变量指标另存为 CSV
  --atol / --rtol             判据阈值（默认 atol=1 K, rtol=1e-3），仅用于结论行

flow3d_node.dat 由 src/sub_IO.f90::output_flow_node 写出（$flow_ec: Iflag_flow_node=1）：
  记录1 = Total_block；记录2 = (ni,nj,nk) 逐块；记录3 = 变量数(=6)
  其后逐块一条数据记录：6*ni*nj*nk 个 double（i 最快、j 次、k 最慢）
  变量顺序：d, u, v, w, T, Ts      T=流体温度[K]（固体块为 Ts），Ts=骨架温度[K]

指标（逐块逐变量）：
  全部节点 / 排除界面节点行 / 排除界面侧第一层单元（2 行节点）
  max|Δ|, rms|Δ|, max|ref|, 相对量 max|Δ|/max|ref|, max|Δ| 位置
"""
import os
import struct
import sys

NVAR = 6
VARNAMES = ["d", "u", "v", "w", "T", "Ts"]


def read_records(path):
    """Fortran unformatted 文件的逐记录 payload（含 4 字节记录长度头尾）。"""
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
    """返回 (dims, fields)：dims=[(ni,nj,nk),...]，fields[(blk,var)]=list（i 最快）。

    文件头（src/sub_IO.f90::output_flow_node，PLOT3D 解(Q)/函数型）：
      记录1 = Total_block
      记录2 = 逐块 3 个整数 (ni,nj,nk)
      记录3 = nvar (=6)
      其后逐块一条数据记录：6*ni*nj*nk 个 double
    （兼容 2026-09-26 工作树曾短暂写出的变体：记录2 = 逐块 4 个整数
      (ni,nj,nk,nvar)、无记录3；当时用旧二进制写出的算例仍可读。）
    """
    recs = read_records(os.path.join(d, "flow3d_node.dat"))
    if len(recs) < 4:
        raise SystemExit("flow3d_node.dat %s: 记录数 %d 太少（需带头）" % (d, len(recs)))
    nb = struct.unpack("<i", recs[0])[0]
    hdr = struct.unpack("<%di" % (len(recs[1]) // 4), recs[1])
    nvars = NVAR
    if len(hdr) == 4 * nb:
        dims = [(hdr[4 * k], hdr[4 * k + 1], hdr[4 * k + 2]) for k in range(nb)]
        nvars = hdr[3]
        data = recs[2:]
    elif len(hdr) == 3 * nb and len(recs[2]) == 4:
        dims = [tuple(hdr[3 * k:3 * k + 3]) for k in range(nb)]
        nvars = struct.unpack("<i", recs[2])[0]
        data = recs[3:]
    else:
        raise SystemExit("flow3d_node.dat %s: 无法识别文件头（记录2 = %d 字节）"
                         % (d, len(recs[1])))
    if nvars != NVAR:
        raise SystemExit("flow3d_node.dat %s: nvars=%d，本脚本需要 6" % (d, nvars))
    if len(data) != nb:
        raise SystemExit("flow3d_node.dat %s: %d 条数据记录，头声明 %d 块" % (d, len(data), nb))
    fields = {}
    for m, raw in enumerate(data):
        ni, nj, nk = dims[m]
        n = ni * nj * nk
        vals = struct.unpack("<%dd" % (len(raw) // 8), raw)
        if len(vals) != NVAR * n:
            raise SystemExit("block %d: 数据长度 %d != 6*%d" % (m + 1, len(vals), n))
        for v in range(NVAR):
            fields[(m + 1, VARNAMES[v])] = vals[v * n:(v + 1) * n]
    return dims, fields


def read_mesh_xyz(d):
    """Mesh3d.x：逐块一条记录 (x,y,z)（节点坐标，SI 米）。返回 [ (x,y,z), ... ]。"""
    recs = read_records(os.path.join(d, "Mesh3d.x"))
    nb = struct.unpack("<i", recs[0])[0]
    dims_all = struct.unpack("<%di" % (len(recs[1]) // 4), recs[1])
    dims = [tuple(dims_all[3 * k:3 * k + 3]) for k in range(nb)]
    xyz = []
    for m in range(nb):
        ni, nj, nk = dims[m]
        n = ni * nj * nk
        a = struct.unpack("<%dd" % (3 * n), recs[2 + m])
        xyz.append((a[0:n], a[n:2 * n], a[2 * n:3 * n]))
    return dims, xyz


def metrics(ref, tst, dims3, mask=None):
    """返回 (max|Δ|, rms|Δ|, max|ref|, (i,j,k), imax)。mask=None 表示全部节点。"""
    ni, nj, nk = dims3
    dmax = 0.0
    s2 = 0.0
    rmax = 0.0
    idx = -1
    cnt = 0
    for n in range(len(ref)):
        if mask is not None and not mask(n):
            continue
        cnt += 1
        dd = abs(ref[n] - tst[n])
        if cnt == 1 or dd > dmax:
            dmax = dd
            idx = n
        s2 += dd * dd
        if abs(ref[n]) > rmax:
            rmax = abs(ref[n])
    rms = (s2 / cnt) ** 0.5 if cnt else 0.0
    if idx < 0:
        return 0.0, 0.0, rmax, (-1, -1, -1), 0
    k = idx // (ni * nj)
    j = (idx - k * ni * nj) // ni
    i = idx - k * ni * nj - j * ni
    return dmax, rms, rmax, (i + 1, j + 1, k + 1), idx


def mask_face_and_cell(dims3, block, iface):
    """返回 (mask_excl_face, mask_excl_cell, 说明文字)。iface: 'j1' / 'i101' 等。"""
    ni, nj, nk = dims3
    axis = iface[0].lower()
    row = int(iface[1:])
    if axis == "i":
        def _ind(i, j, k):
            return i
        nrow = ni
    elif axis == "j":
        def _ind(i, j, k):
            return j
        nrow = nj
    else:
        def _ind(i, j, k):
            return k
        nrow = nk

    def fn(n):
        k = n // (ni * nj)
        j = (n - k * ni * nj) // ni
        i = n - k * ni * nj - j * ni
        return _ind(i + 1, j + 1, k + 1)

    if row <= 1:
        face = lambda n: fn(n) > row
        cell = lambda n: fn(n) > row + 1
        txt = "%s=%d（低侧）" % (axis, row)
    else:
        face = lambda n: fn(n) < row
        cell = lambda n: fn(n) < row - 1
        txt = "%s=%d（高侧）" % (axis, row)
    return face, cell, txt


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = [a for a in sys.argv[1:] if a.startswith("--")]
    if len(args) < 2:
        raise SystemExit(__doc__)
    refd, tstd = args[0].rstrip("/"), args[1].rstrip("/")

    def optval(name, default):
        for i, a in enumerate(sys.argv):
            if a == name:
                return sys.argv[i + 1]
            if a.startswith(name + "="):
                return a.split("=", 1)[1]
        return default

    varlist = optval("--vars", ",".join(VARNAMES)).split(",")
    blocklist = [int(x) for x in optval("--blocks", "").split(",") if x]
    ifspec = optval("--iface", "1:j1,2:i101")
    ifmap = {}
    for item in ifspec.split(","):
        if not item.strip():
            continue
        b, s = item.split(":")
        ifmap[int(b)] = s.strip()
    csvpath = optval("--csv", None)
    atol = float(optval("--atol", "1.0"))
    rtol = float(optval("--rtol", "1e-3"))

    dimsR, fR = read_flow_node(refd)
    dimsT, fT = read_flow_node(tstd)
    if dimsR != dimsT:
        raise SystemExit("两块网格维数不同：%s vs %s" % (dimsR, dimsT))
    _, xyz = read_mesh_xyz(refd)

    print("ref : %s" % refd)
    print("test: %s" % tstd)
    print("=" * 108)
    hdr = ("%-5s %-3s | %-11s %-11s %-11s %-9s | %-17s | %-11s %-11s"
           % ("block", "var", "max|d|", "rms|d|", "max|ref|", "rel", "argmax (i,j,k)",
              "max|d|*", "rms|d|*"))
    print(hdr)
    print("-" * 108)
    nfail = 0
    rows = []
    for m in range(1, len(dimsR) + 1):
        if blocklist and m not in blocklist:
            continue
        ni, nj, nk = dimsR[m - 1]
        iface = ifmap.get(m)
        mface = mcell = None
        itxt = ""
        if iface:
            mface, mcell, itxt = mask_face_and_cell(dimsR[m - 1], m, iface)
        for v in varlist:
            r = fR[(m, v)]
            t = fT[(m, v)]
            d1, r1, big, loc, _ = metrics(r, t, dimsR[m - 1])
            rel = d1 / big if big > 0 else d1
            if mface is not None:
                d2, r2, _, loc2, _ = metrics(r, t, dimsR[m - 1], mface)
                d3, r3, _, _, _ = metrics(r, t, dimsR[m - 1], mcell)
            else:
                d2 = r2 = d3 = r3 = float("nan")
            x, y, z = xyz[m - 1]
            ymax = float("nan")
            if loc[0] > 0:
                ymax = y[(loc[2] - 1) * ni * nj + (loc[1] - 1) * ni + (loc[0] - 1)]
            print("%-5d %-3s | %-11.4g %-11.4g %-11.4g %-9.2g | (%4d,%4d,%4d) y=%-8.3g | %-11.4g %-11.4g"
                  % (m, v, d1, r1, big, rel, loc[0], loc[1], loc[2], ymax, d2, r2))
            rows.append((m, v, d1, r1, big, rel, d2, r2, d3, r3))
            if d1 > atol and rel > rtol and v in ("T", "Ts"):
                nfail += 1
    print("-" * 108)
    print("说明：'*' 列为排除界面节点行后的指标；'rel' = max|d| / max|ref|")
    for m in sorted(ifmap):
        print("      块%d 界面 %s" % (m, ifmap[m]))
    print("结论参考：atol=%.3g K(k/J/kg), rtol=%.3g -> 超限的 (块,变量) 数 = %d"
          % (atol, rtol, nfail))
    if csvpath:
        with open(csvpath, "w") as f:
            f.write("block,var,max_abs,rms_abs,max_ref,rel,excl_iface_max,excl_iface_rms,"
                    "excl_1cell_max,excl_1cell_rms\n")
            for r in rows:
                f.write("%d,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n" % r)
        print("CSV -> %s" % csvpath)


if __name__ == "__main__":
    main()

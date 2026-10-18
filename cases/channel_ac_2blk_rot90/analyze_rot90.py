#!/usr/bin/env python3
"""i<->j (90-degree index rotation) two-block channel check.

The physical grid is the same as cases/channel_ac_2blk (80+80 cells in x,
40 in y), but block 2 is written TRANSPOSED (i -> y, j -> x), so the
interface is  blk-1 i+ <-> blk-2 j-  and Convert_bc has to build a
descriptor with swapped index axes, L = (2, +1, -3) for the 'p' variant
(blk-2 i runs along +y, its k along -z so the triad stays right-handed) and
L = (2, -1, +3) for the 'n' variant (blk-2 i runs along -y).

This script therefore cannot concatenate blocks along an index axis: it
reads every block, works out which index runs along x / y from the VTK point
order, and scatters the cell data onto the physical (x, y) grid derived from
the cell centres.  On that grid it reports

  1. grid completeness (every one of the 160x40 cells covered exactly once),
  2. u_max / u_min and column sums of u (inlet, outlet) -> global mass balance,
  3. the jump of u, v, p across the interface plane x = Lx/2,
  4. if a reference run directory is given: RMS(u), max|du|, max|dp| and the
     mean pressure offset against it.

Usage:  python3 analyze_rot90.py <dir> [<ref_dir>]
"""
import sys
import glob
import numpy as np

LX, H = 8.0, 1.0


def _num(tok):
    """Fortran may print the exponent without the 'E' (e.g. '1.234-108')."""
    try:
        return float(tok)
    except ValueError:
        for k in range(1, len(tok) - 1):
            if tok[k] in '+-' and tok[k-1] not in 'eEdD' and tok[k-1].isdigit():
                return float(tok[:k] + 'E' + tok[k:])
        raise


def read_vtk(fn):
    """Cell-centred fields + cell-centre x/y + 'i runs along y' flag."""
    lines = open(fn).read().splitlines()

    def find(pre):
        for n, ln in enumerate(lines):
            if ln.startswith(pre):
                return n
        raise SystemExit('cannot find %r in %s' % (pre, fn))

    p0 = find('POINTS')
    npts = int(lines[p0].split()[1])
    XYZ = np.array([[float(x) for x in l.split()]
                    for l in lines[p0+1:p0+1+npts]])
    # the VTK point list has i fastest: whichever coordinate moves between
    # the first two points is the i direction -> decides nx / ny of the reshape
    d = XYZ[1] - XYZ[0]
    i_is_y = abs(d[1]) > abs(d[0])
    if not i_is_y:
        nx = len(np.unique(np.round(XYZ[:, 0], 10)))
        ny = len(np.unique(np.round(XYZ[:, 1], 10)))
    else:
        nx = len(np.unique(np.round(XYZ[:, 1], 10)))
        ny = len(np.unique(np.round(XYZ[:, 0], 10)))
    nz = npts // (nx * ny)
    P = XYZ.reshape(nz, ny, nx, 3)                      # (k, j, i, xyz)
    if not i_is_y:          # i -> x, j -> y  (normal block)
        xc = 0.5 * (P[0, 0, 1:, 0] + P[0, 0, :-1, 0])   # (nx-1,) along i
        yc = 0.5 * (P[0, 1:, 0, 1] + P[0, :-1, 0, 1])   # (ny-1,) along j
    else:                   # i -> y, j -> x  (transposed block)
        xc = 0.5 * (P[0, 1:, 0, 0] + P[0, :-1, 0, 0])   # (ny-1,) along j
        yc = 0.5 * (P[0, 0, 1:, 1] + P[0, 0, :-1, 1])   # (nx-1,) along i

    cd = find('CELL_DATA')
    ncel = int(lines[cd].split()[1])
    marks = [(n, ln.split()[0], ln.split()[1]) for n, ln in enumerate(lines)
             if ln.startswith(('SCALARS', 'VECTORS'))]
    D = {}
    for ii in range(len(marks)):
        n0, typ, nm = marks[ii]
        n1 = marks[ii+1][0] if ii+1 < len(marks) else len(lines)
        body = [_num(x) for l in lines[n0+1:n1] for x in l.split()
                if not l.startswith('LOOKUP')]
        if typ == 'SCALARS':
            D[nm] = np.array(body[:ncel]).reshape(nz-1, ny-1, nx-1)[0]
        else:
            D[nm] = np.array(body[:3*ncel]).reshape(nz-1, ny-1, nx-1, 3)[0]
    return xc, yc, D, i_is_y



def _axis_index(centres, allc):
    """Map cell-centre values to global columns (tolerant of 1e-15 noise)."""
    idx = []
    for c in centres:
        j = int(np.argmin(np.abs(allc - c)))
        if abs(allc[j] - c) > 1e-6:
            raise SystemExit('cell centre %.12g not on the global grid' % c)
        idx.append(j)
    return np.array(idx)


def load_run(dirn):
    files = sorted(glob.glob(dirn.rstrip('/') + '/flow3d_block_*.vtk'))
    if not files:
        raise SystemExit('no flow3d_block_*.vtk in ' + dirn)
    blocks = [read_vtk(f) for f in files]
    xs = np.unique(np.concatenate([np.round(b[0], 9) for b in blocks]))
    ys = np.unique(np.concatenate([np.round(b[1], 9) for b in blocks]))
    ny, nx = ys.size, xs.size
    U = np.full((ny, nx), np.nan)
    V = np.full((ny, nx), np.nan)
    W = np.full((ny, nx), np.nan)
    filled = np.zeros((ny, nx), dtype=int)
    sizes = []
    for xc, yc, D, i_is_y in blocks:
        sizes.append(xc.size)
        ix = _axis_index(xc, xs)
        jy = _axis_index(yc, ys)
        vel = D['velocity']
        # U(2:4) are CARTESIAN velocity components, independent of the block
        # index orientation: the VTK 'velocity' vector is (u,v,w) in both
        # layouts.  Only the array axes differ: for the transposed block the
        # stored array is (j,i) = (x,y), so transpose it to (y,x).
        uu = vel[:, :, 0].T if i_is_y else vel[:, :, 0]
        vv = vel[:, :, 1].T if i_is_y else vel[:, :, 1]
        U[np.ix_(jy, ix)] = uu
        V[np.ix_(jy, ix)] = vv
        if 'pressure' in D:
            W[np.ix_(jy, ix)] = D['pressure'].T if i_is_y else D['pressure']
        filled[np.ix_(jy, ix)] += 1
    return xs, ys, U, V, W, filled, sizes


def main():
    d = sys.argv[1]
    xs, ys, U, V, P, filled, sizes = load_run(d)
    print('%s : grid %dx%d cells, block cell counts (along their i) %s'
          % (d, xs.size, ys.size, sizes))
    print('   coverage: min/max writes per cell = %d / %d %s'
          % (filled.min(), filled.max(),
             'OK (each cell exactly once)'
             if filled.min() == 1 and filled.max() == 1
             else '!!! HOLE or OVERLAP'))
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]
    print('   dx=%.6g dy=%.6g ; x in [%.4f,%.4f] y in [%.4f,%.4f]'
          % (dx, dy, xs[0]-dx/2, xs[-1]+dx/2, ys[0]-dy/2, ys[-1]+dy/2))
    print('   u_max = %+.6f  u_min = %+.6f  v_max = %+.6f  v_min = %+.6f'
          % (U.max(), U.min(), V.max(), V.min()))
    print('   column sum u: inlet %.6f  outlet %.6f  (analytic 40 / 40)'
          % (U[:, 0].sum() * dy, U[:, -1].sum() * dy))
    print('   max |v| on the walls = %.3e / %.3e'
          % (np.abs(V[0, :]).max(), np.abs(V[-1, :]).max()))
    # interface plane x = LX/2 : the columns just left / right of it
    ii = int(np.argmin(np.abs(xs - (LX / 2 - dx / 2))))
    io = int(np.argmin(np.abs(xs - (LX / 2 + dx / 2))))
    print('   interface x=%.2f (cols %d|%d): max|du|=%.3e max|dv|=%.3e '
          'max|dp|=%.3e' % (LX / 2, ii, io,
                            np.max(np.abs(U[:, io] - U[:, ii])),
                            np.max(np.abs(V[:, io] - V[:, ii])),
                            np.nanmax(np.abs(P[:, io] - P[:, ii]))))
    print('   interface column sum: left %.6f  right %.6f  (diff %.3e)'
          % (U[:, ii].sum() * dy, U[:, io].sum() * dy,
             abs(U[:, ii].sum() - U[:, io].sum()) * dy))
    if len(sys.argv) > 2:
        x2, y2, U2, V2, P2, filled2, _ = load_run(sys.argv[2])
        if U.shape == U2.shape:
            print('   vs %s : RMS(u)=%.3e RMS(v)=%.3e max|du|=%.3e '
                  'max|dv|=%.3e max|dp|=%.3e mean(p-p2)=%+.3e'
                  % (sys.argv[2], np.sqrt(np.mean((U - U2) ** 2)),
                     np.sqrt(np.mean((V - V2) ** 2)),
                     np.max(np.abs(U - U2)), np.max(np.abs(V - V2)),
                     np.nanmax(np.abs(P - P2)), np.nanmean(P - P2)))
        else:
            print('   vs %s : grid mismatch %s vs %s'
                  % (sys.argv[2], U.shape, U2.shape))


if __name__ == '__main__':
    main()

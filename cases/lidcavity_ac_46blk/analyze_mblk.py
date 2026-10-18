#!/usr/bin/env python3
"""4/6-block lid-driven cavity check (AC solver, LS_Algorithm=3).

Every block of these layouts has i -> x and j -> y, so the block cell data can
be scattered onto the physical 42x42 cell grid (uniform: h = 1/42):

  m4 : 2 x 2 blocks of 21x21 cells     (block boundaries at x = y = 1/2)
  m6 : 3 x 2 blocks of 14x21 cells     (boundaries at x = 1/3, 2/3 and y = 1/2)

The layout is taken from the run directory name (m4*/m6*).  The script reports

  1. grid completeness (each of the 42x42 cells written exactly once),
  2. u(y) at x = 0.5 and v(x) at y = 0.5 against Ghia et al. (1982) Re=100,
  3. the velocity/pressure jump across every interior block boundary,
  4. if a reference directory is given: RMS(u), RMS(v), max|du|, max|dp| and
     the mean pressure offset -> this is the anchoring / numbering check
     (velocity must match, p may differ by a constant) and the np-independence
     check (everything must match).

Usage:  python3 analyze_mblk.py <dir> [<ref_dir>]
"""
import sys
import os
import re
import glob
import numpy as np

GH100 = {
 'y': [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,
       0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000],
 'u': [1.00000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,
       -0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717,0.0],
 'x': [0.0,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,
       0.8516,0.9531,0.9609,0.9688,0.9766,1.0],
 'v': [0.0,0.09233,0.10091,0.15358,0.17442,0.17386,0.13860,0.05918,
       -0.05918,-0.13860,-0.17386,-0.17442,-0.15358,-0.10091,-0.09233,0.0],
}


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
    """Cell-centred fields + cell-centre x/y (blocks here are always i->x)."""
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
    nx = len(np.unique(np.round(XYZ[:, 0], 10)))
    ny = len(np.unique(np.round(XYZ[:, 1], 10)))
    nz = npts // (nx * ny)
    P = XYZ.reshape(nz, ny, nx, 3)
    xc = 0.5 * (P[0, 0, 1:, 0] + P[0, 0, :-1, 0])
    yc = 0.5 * (P[0, 1:, 0, 1] + P[0, :-1, 0, 1])

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
    return xc, yc, D


def _axis_index(centres, allc):
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
    for xc, yc, D in blocks:
        ix = _axis_index(xc, xs)
        jy = _axis_index(yc, ys)
        U[np.ix_(jy, ix)] = D['velocity'][:, :, 0]
        V[np.ix_(jy, ix)] = D['velocity'][:, :, 1]
        if 'pressure' in D:
            W[np.ix_(jy, ix)] = D['pressure']
        filled[np.ix_(jy, ix)] += 1
    return xs, ys, U, V, W, filled


def main():
    d = sys.argv[1]
    layout = 'm2' if re.search(r'm2', os.path.basename(os.path.normpath(d))) \
        else ('m4' if re.search(r'm4', os.path.basename(os.path.normpath(d))) else 'm6')
    xs, ys, U, V, P, filled = load_run(d)
    print('%s : layout %s, grid %dx%d cells' % (d, layout, xs.size, ys.size))
    print('   coverage: min/max writes per cell = %d / %d %s'
          % (filled.min(), filled.max(),
             'OK' if filled.min() == 1 and filled.max() == 1
             else '!!! HOLE or OVERLAP'))
    g = GH100
    ic = int(np.argmin(np.abs(xs - 0.5)))
    jc = int(np.argmin(np.abs(ys - 0.5)))
    up = np.interp(g['y'], ys, U[:, ic])
    vp = np.interp(g['x'], xs, V[jc, :])
    print('   u(y) at x=0.5: RMS vs Ghia = %.5f  (min %+.5f, Ghia %+.5f)'
          % (np.sqrt(np.mean((up - g['u']) ** 2)), U[:, ic].min(), min(g['u'])))
    print('   v(x) at y=0.5: RMS vs Ghia = %.5f  (max %+.5f, Ghia %+.5f)'
          % (np.sqrt(np.mean((vp - g['v']) ** 2)), V[jc, :].max(), max(g['v'])))
    print('   u range [%+.5f, %+.5f]  v range [%+.5f, %+.5f]  p range [%+.4e, %+.4e]'
          % (U.min(), U.max(), V.min(), V.max(),
             np.nanmin(P), np.nanmax(P)))
    # interior block boundaries of the layout, in cell units
    if layout == 'm2':
        cuts_x, cuts_y = [21], []
    elif layout == 'm4':
        cuts_x, cuts_y = [21], [21]
    else:
        cuts_x, cuts_y = [14, 28], [21]
    for c in cuts_x:
        du = np.max(np.abs(U[:, c] - U[:, c-1]))
        dv = np.max(np.abs(V[:, c] - V[:, c-1]))
        dp = np.nanmax(np.abs(P[:, c] - P[:, c-1]))
        print('   cut x=%.4f (cols %d|%d): max|du|=%.3e max|dv|=%.3e max|dp|=%.3e'
              % (xs[c], c-1, c, du, dv, dp))
    for c in cuts_y:
        du = np.max(np.abs(U[c, :] - U[c-1, :]))
        dv = np.max(np.abs(V[c, :] - V[c-1, :]))
        dp = np.nanmax(np.abs(P[c, :] - P[c-1, :]))
        print('   cut y=%.4f (rows %d|%d): max|du|=%.3e max|dv|=%.3e max|dp|=%.3e'
              % (ys[c], c-1, c, du, dv, dp))
    if len(sys.argv) > 2:
        x2, y2, U2, V2, P2, filled2 = load_run(sys.argv[2])
        if U.shape == U2.shape:
            print('   vs %s : RMS(u)=%.3e RMS(v)=%.3e max|du|=%.3e max|dv|=%.3e '
                  'max|dp|=%.3e mean(p-p2)=%+.3e'
                  % (sys.argv[2], np.sqrt(np.mean((U - U2) ** 2)),
                     np.sqrt(np.mean((V - V2) ** 2)), np.max(np.abs(U - U2)),
                     np.max(np.abs(V - V2)), np.nanmax(np.abs(P - P2)),
                     np.nanmean(P - P2)))
        else:
            print('   vs %s : grid mismatch %s vs %s'
                  % (sys.argv[2], U.shape, U2.shape))


if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""3-block lid-driven cavity check (AC solver, LS_Algorithm=3).

The case is an equal x-split of a 42x42 cell cavity into 3 blocks of 14x42
cells, so the three flow3d_block_*.vtk files tile the domain in x and share
the same y grid.  This script

  1. loads all blocks of a run directory, concatenates them in +x order into a
     single 42x42 cell grid (so np=1 / np=2 / np=3 and the single-block
     reference case/lidcavity_ac/re100 are all compared on the same grid),
  2. prints u(x=0.5) and v(y=0.5) against Ghia et al. (1982) Table I/II,
  3. if a reference directory is given, prints the grid-to-grid RMS of u and v
     against it (layout equivalence check),
  4. prints the velocity/pressure jump across the two block interfaces.

Usage:  python3 analyze_3blk.py <dir> [<ref_dir>]
"""
import sys
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


def read_vtk(fn):
    """Return dict of cell-centred fields + cell-centre x/y of one block."""
    lines = open(fn).read().splitlines()
    def find(pre):
        for n, ln in enumerate(lines):
            if ln.startswith(pre):
                return n
        raise SystemExit('cannot find %r in %s' % (pre, fn))
    p0 = find('POINTS')
    npts = int(lines[p0].split()[1])
    XYZ = np.array([[float(x) for x in l.split()] for l in lines[p0+1:p0+1+npts]])
    nx = len(np.unique(np.round(XYZ[:, 0], 10)))
    ny = len(np.unique(np.round(XYZ[:, 1], 10)))
    nz = npts // (nx * ny)
    P = XYZ.reshape(nz, ny, nx, 3)                     # (k, j, i, xyz)
    xc = 0.5*(P[0, 0, 1:, 0] + P[0, 0, :-1, 0])        # (nx-1,)
    yc = 0.5*(P[0, 1:, 0, 1] + P[0, :-1, 0, 1])        # (ny-1,)
    cd = find('CELL_DATA')
    ncel = int(lines[cd].split()[1])
    marks = [(n, ln.split()[0], ln.split()[1]) for n, ln in enumerate(lines)
             if ln.startswith(('SCALARS', 'VECTORS'))]
    D = {}
    for ii in range(len(marks)):
        n0, typ, nm = marks[ii]
        n1 = marks[ii+1][0] if ii+1 < len(marks) else len(lines)
        body = [float(x) for l in lines[n0+1:n1] for x in l.split()
                if not l.startswith('LOOKUP')]
        if typ == 'SCALARS':
            D[nm] = np.array(body[:ncel]).reshape(nz-1, ny-1, nx-1)[0]
        else:
            D[nm] = np.array(body[:3*ncel]).reshape(nz-1, ny-1, nx-1, 3)[0]
    return xc, yc, D


def load_run(dirn):
    files = sorted(glob.glob(dirn.rstrip('/') + '/flow3d_block_*.vtk'))
    if not files:
        raise SystemExit('no flow3d_block_*.vtk in ' + dirn)
    blocks = []
    for fn in files:
        xc, yc, D = read_vtk(fn)
        blocks.append((xc, yc, D))
    # order the blocks by their minimum x centre
    blocks.sort(key=lambda b: b[0][0])
    xc = np.concatenate([b[0] for b in blocks])
    yc = blocks[0][1]
    order = np.argsort(xc)
    xc = xc[order]
    U = np.concatenate([b[2]['velocity'][:, :, 0] for b in blocks], axis=1)[:, order]
    V = np.concatenate([b[2]['velocity'][:, :, 1] for b in blocks], axis=1)[:, order]
    P = None
    for key in ('pressure', 'p', 'Pressure'):
        if key in blocks[0][2]:
            P = np.concatenate([b[2][key] for b in blocks], axis=1)[:, order]
            break
    ifc = []
    nb = 0
    for b in blocks[:-1]:
        nb += b[0].size
        ifc.append(nb)          # number of columns left of this interface
    return xc, yc, U, V, P, ifc


def main():
    d = sys.argv[1]
    xc, yc, U, V, P, ifc = load_run(d)
    g = GH100
    ic = int(np.argmin(np.abs(xc-0.5)))
    jc = int(np.argmin(np.abs(yc-0.5)))
    up = np.interp(g['y'], yc, U[:, ic])
    vp = np.interp(g['x'], xc, V[jc, :])
    eu = np.sqrt(np.mean((up-g['u'])**2))
    ev = np.sqrt(np.mean((vp-g['v'])**2))
    print('%s : grid %dx%d, blocks end at columns %s'
          % (d, U.shape[1], U.shape[0], ifc))
    print('   u(x=0.5): RMS vs Ghia = %.5f, min = %+.5f (Ghia %+.5f)'
          % (eu, U[:, ic].min(), min(g['u'])))
    print('   v(y=0.5): RMS vs Ghia = %.5f, max = %+.5f (Ghia %+.5f)'
          % (ev, V[jc, :].max(), max(g['v'])))
    print('   vertical u(x=0.5):   y / sim / Ghia')
    for k in [1, 3, 5, 7, 9, 10, 12, 14]:
        print('     %.4f  %+.5f %+.5f' % (g['y'][k], up[k], g['u'][k]))
    print('   horizontal v(y=0.5): x / sim / Ghia')
    for k in [1, 3, 5, 7, 8, 10, 12, 14]:
        print('     %.4f  %+.5f %+.5f' % (g['x'][k], vp[k], g['v'][k]))
    # interface jumps (a converged multi-block solution must be continuous)
    for icf in ifc:
        du = np.max(np.abs(U[:, icf] - U[:, icf-1]))
        dv = np.max(np.abs(V[:, icf] - V[:, icf-1]))
        msg = '   interface after column %2d (%d blocks left): max|du|=%.3e max|dv|=%.3e' \
              % (icf, icf//(ifc[0]), du, dv)
        if P is not None:
            dp = np.max(np.abs(P[:, icf] - P[:, icf-1]))
            msg += ' max|dp|=%.3e' % dp
        print(msg)
    if len(sys.argv) > 2:
        x2, y2, U2, V2, P2, _ = load_run(sys.argv[2])
        if U.shape == U2.shape:
            print('   vs %s : RMS(u)=%.3e  RMS(v)=%.3e  max|du|=%.3e'
                  % (sys.argv[2], np.sqrt(np.mean((U-U2)**2)),
                     np.sqrt(np.mean((V-V2)**2)), np.max(np.abs(U-U2))))
        else:
            print('   vs %s : grid mismatch %s vs %s'
                  % (sys.argv[2], U.shape, U2.shape))


if __name__ == '__main__':
    main()

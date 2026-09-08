#!/usr/bin/env python3
"""Compare lid-cavity centreline profiles vs Ghia et al. (1982) Re=1000.
Usage: python3 compare_profiles.py <run_dir> [<run_dir> ...]
Each run_dir must contain Mesh3d.x and flow3d.dat.
"""
import struct
import sys
import numpy as np


def rd(f):
    n = struct.unpack('I', f.read(4))[0]
    d = f.read(n)
    struct.unpack('I', f.read(4))
    return d


GHIA_U = np.array([
    [1.0000, 1.00000], [0.9766, 0.65928], [0.9688, 0.57492], [0.9609, 0.51117],
    [0.9531, 0.46604], [0.8516, 0.33304], [0.7344, 0.18719], [0.6172, 0.05702],
    [0.5000, -0.06080], [0.4531, -0.10648], [0.2813, -0.27805], [0.1719, -0.38289],
    [0.1016, -0.29730], [0.0703, -0.22220], [0.0625, -0.20196], [0.0547, -0.18109],
    [0.0000, 0.00000]])
GHIA_V = np.array([
    [0.0000, 0.00000], [0.0625, 0.09245], [0.0703, 0.10091], [0.1016, 0.15444],
    [0.1719, 0.17570], [0.2813, 0.17527], [0.4531, 0.13690], [0.5000, 0.06080],
    [0.6172, -0.05702], [0.7344, -0.18719], [0.8516, -0.33304], [0.9531, -0.46604],
    [0.9609, -0.51117], [0.9688, -0.57492], [0.9766, -0.65928], [1.0000, 0.00000]])


def load(run_dir):
    with open(f'{run_dir}/Mesh3d.x', 'rb') as f:
        rd(f)
        dims = np.frombuffer(rd(f), dtype=np.int32)
        nx, ny, nz = dims
        x = np.frombuffer(rd(f), dtype=np.float64).reshape(tuple(dims), order='F')
        y = np.frombuffer(rd(f), dtype=np.float64).reshape(tuple(dims), order='F')
    with open(f'{run_dir}/flow3d.dat', 'rb') as f:
        n = struct.unpack('I', f.read(4))[0]
        arr = np.frombuffer(f.read(n), dtype=np.float64)
    # flow3d stores (0:nx,0:ny,0:nz+1,5) for these grids (3rd dim = 3)
    per = arr.size // 5
    nz3 = per // ((nx + 1) * (ny + 1))
    U = arr.reshape((nx + 1, ny + 1, nz3, 5), order='F')
    u = U[:, :, 1, 1]
    v = U[:, :, 1, 2]
    # cell centres
    xc = 0.5 * (x[:-1, :-1, 0] + x[1:, 1:, 0])
    yc = 0.5 * (y[:-1, :-1, 0] + y[1:, 1:, 0])
    return xc, yc, u, v


def profiles(run_dir):
    xc, yc, u, v = load(run_dir)
    ic = np.argmin(np.abs(xc[:, 0] - 0.5))
    jc = np.argmin(np.abs(yc[0, :] - 0.5))
    us = np.interp(GHIA_U[:, 0], yc[ic, :], u[ic + 1, 1:-1])
    vs = np.interp(GHIA_V[:, 0], xc[:, jc], v[1:-1, jc + 1])
    return us, vs


def main():
    rows = []
    all_us, all_vs = {}, {}
    for d in sys.argv[1:]:
        us, vs = profiles(d)
        all_us[d], all_vs[d] = us, vs
        eu = max(abs(us - GHIA_U[:, 1]))
        ev = max(abs(vs - GHIA_V[:, 1]))
        ru = np.sqrt(np.mean((us - GHIA_U[:, 1]) ** 2))
        rv = np.sqrt(np.mean((vs - GHIA_V[:, 1]) ** 2))
        rows.append((d, eu, ev, ru, rv))
        print(f'{d:28s}  u: max={eu:.4f} rms={ru:.4f}   v: max={ev:.4f} rms={rv:.4f}')
    print('\nu(x=0.5,y):  sim vs Ghia')
    print('     y      Ghia     ' + '  '.join(f'{d.split("/")[-1]:>10s}' for d in sys.argv[1:]))
    for k in [1, 3, 5, 7, 8, 10, 11, 13, 14]:
        line = f'  {GHIA_U[k,0]:.4f}  {GHIA_U[k,1]:+.5f}  '
        line += '  '.join(f'{all_us[d][k]:+.5f}' for d in sys.argv[1:])
        print(line)
    print('\nv(x,y=0.5):  sim vs Ghia')
    print('     x      Ghia     ' + '  '.join(f'{d.split("/")[-1]:>10s}' for d in sys.argv[1:]))
    for k in [1, 3, 5, 7, 8, 10, 11, 13, 14]:
        line = f'  {GHIA_V[k,0]:.4f}  {GHIA_V[k,1]:+.5f}  '
        line += '  '.join(f'{all_vs[d][k]:+.5f}' for d in sys.argv[1:])
        print(line)


if __name__ == '__main__':
    main()

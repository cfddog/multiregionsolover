#!/usr/bin/env python3
# Compare cavity u-velocity along vertical centerline (x=0.5) vs Ghia et al. (1982), Re=100
import sys
import numpy as np

# Ghia, Ghia & Shin (1982), TABLE I: u-velocity along vertical line x=0.5, Re=100
ghia_y = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
          0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
ghia_u = [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641,
          -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000]


def read_vtk(path):
    lines = open(path).read().splitlines()
    fields, cur, buf = {}, None, []
    i = 0
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith(('SCALARS', 'VECTORS')):
            if cur:
                fields[cur] = np.array(buf)
            cur, buf = l.split()[1], []
            i += 1
            if i < len(lines) and lines[i].strip().startswith('LOOKUP_TABLE'):
                i += 1
            continue
        if cur is not None:
            try:
                buf += [float(t) for t in l.split()]
            except ValueError:
                fields[cur] = np.array(buf)
                cur, buf = None, []
        i += 1
    if cur:
        fields[cur] = np.array(buf)
    return fields


def u_profile(path, n=40):
    # structured grid, cells ordered i(x)-fastest, then j(y), then k
    f = read_vtk(path)
    vel = f['velocity'].reshape(-1, 3)  # (nx-1)*(ny-1)*(nz-1) cells
    u = vel[:, 0].reshape(n, n)         # k=1 layer only (nz-1=1)
    # vertical centerline x=0.5: average the two cells adjacent to x=0.5 (i=20,21 -> idx 19,20)
    return u[:, 19] * 0.5 + u[:, 20] * 0.5, np.arange(n) / n + 0.5 / n


def main(paths_labels):
    for path, label in paths_labels:
        prof, yc = u_profile(path)
        # interpolate Ghia onto our cell-center y (Ghia y is ascending in table? given descending)
        gy = np.array(ghia_y)[::-1]
        gu = np.array(ghia_u)[::-1]
        ui = np.interp(yc, gy, gu)
        err = prof - ui
        rms = np.sqrt(np.mean(err**2))
        # key locations
        umin = prof.min()
        print(f"--- {label} ---")
        print(f"  RMS error vs Ghia (40 pts): {rms:.5f}   max|err|: {np.abs(err).max():.5f}")
        print(f"  u_min (vortex backflow): {umin:+.5f}  (Ghia -0.21090)  err {umin+0.21090:+.5f}")
        j = np.argmin(prof)
        print(f"  u at y=0.5: {prof[19]:+.5f} (Ghia -0.20581)   u at y~0.95: {prof[37]:+.5f}")
    print()


if __name__ == '__main__':
    args = sys.argv[1:]
    pairs = [(args[i], args[i + 1]) for i in range(0, len(args), 2)]
    main(pairs)

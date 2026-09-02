#!/usr/bin/env python3
"""Generate a lid-driven cavity mesh (unit square, quasi-2D) with optional
tangent-hyperbolic boundary-layer clustering on walls.

Output: Mesh3d.x (Fortran unformatted, PLOT3D single block).
Default: 41 x 41 x 2 nodes  ->  40 x 40 x 1 cells.
Clustering controlled by CLI or constants:
  stretch_x/y: True  -> tangent-hyperbolic clustering toward both wall ends
  beta_x/y   : stretching parameter (0 = uniform; typical 0.6..1.4 for mild/strong)
"""
import struct
import array
import math
import argparse


def tangent_hyperbolic(n: int, L: float, beta: float) -> list[float]:
    """Return `n` node coordinates in [0, L] symmetrically clustered toward
    both ends (0 and L) using the tangent-hyperbolic stretching function of
    Vinokur (1983). beta=0 recovers uniform spacing; larger beta squeezes more
    points toward the walls.
    """
    if n <= 1 or beta <= 1e-12:
        return [L * i / max(n - 1, 1) for i in range(n)]
    coords = []
    for i in range(n):
        xi = i / (n - 1)                        # logical in [0,1]
        eta = 0.5 + 0.5 * math.tanh(beta * (xi - 0.5)) / math.tanh(beta / 2.0)
        coords.append(L * eta)
    return coords


def write_fortran_record(f, data, fmt):
    if isinstance(data, array.array):
        raw = data.tobytes()
    else:
        raw = struct.pack(fmt, *data)
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nx', type=int, default=41)
    ap.add_argument('--ny', type=int, default=41)
    ap.add_argument('--nz', type=int, default=2)
    ap.add_argument('--Lx', type=float, default=1.0)
    ap.add_argument('--Ly', type=float, default=1.0)
    ap.add_argument('--dz', type=float, default=1.0)
    ap.add_argument('--stretch_x', action='store_true')
    ap.add_argument('--stretch_y', action='store_true')
    ap.add_argument('--beta_x', type=float, default=1.2)
    ap.add_argument('--beta_y', type=float, default=1.2)
    args = ap.parse_args()

    xs = tangent_hyperbolic(args.nx, args.Lx, args.beta_x if args.stretch_x else 0.0)
    ys = tangent_hyperbolic(args.ny, args.Ly, args.beta_y if args.stretch_y else 0.0)

    x_coords, y_coords, z_coords = [], [], []
    for k in range(args.nz):
        zk = k * args.dz
        for j in range(args.ny):
            yj = ys[j]
            for i in range(args.nx):
                x_coords.append(xs[i])
                y_coords.append(yj)
                z_coords.append(zk)

    with open("Mesh3d.x", "wb") as f:
        write_fortran_record(f, [1], 'i')
        write_fortran_record(f, [args.nx, args.ny, args.nz], 'iii')
        write_fortran_record(f, array.array('d', x_coords), '')
        write_fortran_record(f, array.array('d', y_coords), '')
        write_fortran_record(f, array.array('d', z_coords), '')

    info = [f"Mesh3d.x generated: {args.nx}x{args.ny}x{args.nz} nodes"]
    for axis, flag, beta, nodes, L in [("x", args.stretch_x, args.beta_x, xs, args.Lx),
                                       ("y", args.stretch_y, args.beta_y, ys, args.Ly)]:
        if flag:
            import numpy as np
            nn = np.asarray(nodes)
            h_min = float(nn[1] - nn[0])
            h_max = float((nn[1:] - nn[:-1]).max())
            info.append(f"  {axis}-clustered (beta={beta}):  h_min={h_min:.4f}  h_max={h_max:.4f}  (L={L})")
        else:
            info.append(f"  {axis}-uniform: h={L/(args.nx if axis=='x' else args.ny-1):.4f}")
    print("\n".join(info))


if __name__ == "__main__":
    main()

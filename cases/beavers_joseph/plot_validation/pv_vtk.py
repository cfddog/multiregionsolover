#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generic loader for OpenCFD-EC structured VTK output (ASCII).

Reads POINTS (nodes, i-fastest ordering) and the cell-centred fields
density / velocity / temperature / pressure / Ts written by output_vtk.
Returns arrays shaped (nyc, nxc) for the k=1 cell plane (quasi-2D), plus
cell-centre coordinates computed from the node block.

Numbering: each block's node file has DIMENSIONS nx ny nz; cell data count
(nx-1)(ny-1)(nz-1), ordered k outer, j, i (i fastest).
"""
import numpy as np


def _read_section(lines, start):
    """Return list of floats following a header starting at 'start'."""
    vals = []
    i = start
    while i < len(lines):
        t = lines[i].split()
        if not t:
            i += 1
            continue
        try:
            vals.append(float(t[0]))
        except ValueError:
            break
        i += 1
    return vals, i


def load_block(path):
    lines = open(path).read().splitlines()
    # dimensions
    for l in lines:
        if l.startswith('DIMENSIONS'):
            _, nx, ny, nz = l.split()
            nx, ny, nz = int(nx), int(ny), int(nz)
            break
    ipts = 1
    for idx, l in enumerate(lines):
        if l.startswith('POINTS'):
            ipts = idx
            break
    npts = nx * ny * nz
    coords = []
    for k in range(ipts + 1, min(ipts + 1 + npts, len(lines))):
        t = lines[k].split()
        if len(t) >= 3:
            coords.append((float(t[0]), float(t[1]), float(t[2])))
        if len(coords) == npts:
            break
    xs = np.array([c[0] for c in coords]).reshape((nz, ny, nx), order='C')
    ys = np.array([c[1] for c in coords]).reshape((nz, ny, nx), order='C')
    zs = np.array([c[2] for c in coords]).reshape((nz, ny, nx), order='C')
    # cell centres (average of the 2x2x2 corners of each cell)
    xc = 0.125 * (xs[:-1, :-1, :-1] + xs[1:, :-1, :-1] + xs[:-1, 1:, :-1] +
                  xs[1:, 1:, :-1] + xs[:-1, :-1, 1:] + xs[1:, :-1, 1:] +
                  xs[:-1, 1:, 1:] + xs[1:, 1:, 1:])
    yc = 0.125 * (ys[:-1, :-1, :-1] + ys[1:, :-1, :-1] + ys[:-1, 1:, :-1] +
                  ys[1:, 1:, :-1] + ys[:-1, :-1, 1:] + ys[1:, :-1, 1:] +
                  ys[:-1, 1:, 1:] + ys[1:, 1:, 1:])
    nc = (nx - 1) * (ny - 1) * (nz - 1)

    def scalar(name):
        for i, l in enumerate(lines):
            if l.startswith('SCALARS %s double 1' % name):
                j = next(k for k in range(i + 1, len(lines))
                         if lines[k].startswith('LOOKUP_TABLE'))
                vals, _ = _read_section(lines, j + 1)
                v = np.array(vals[:nc]).reshape((nz - 1, ny - 1, nx - 1))
                return v[0]                 # quasi-2D plane k=0
        raise KeyError('scalar ' + name)

    def vector(name):
        for i, l in enumerate(lines):
            if l.startswith('VECTORS %s double' % name):
                vals = []
                for k in range(i + 1, len(lines)):
                    t = lines[k].split()
                    if not t:
                        continue
                    try:
                        vals += [float(x) for x in t]
                    except ValueError:
                        break
                a = np.array(vals[:3 * nc]).reshape((nc, 3))
                v = a.reshape((nz - 1, ny - 1, nx - 1, 3))
                return v[0]
        raise KeyError('vector ' + name)

    rho = scalar('density')
    T = scalar('temperature')
    p = scalar('pressure')
    vel = vector('velocity')
    return dict(nx=nx, ny=ny, nz=nz, xc=xc[0], yc=yc[0], rho=rho,
                p=p, T=T, u=vel[:, :, 0], v=vel[:, :, 1], w=vel[:, :, 2])

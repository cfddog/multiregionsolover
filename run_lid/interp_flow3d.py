#!/usr/bin/env python3
"""Interpolate coarse (41x41x2) flow3d.dat onto fine (81x81x2) mesh.

flow3d.dat stores U(0:nx, 0:ny, 0:nz, 5) where nx = PLOT3D ni.
PLOT3D mesh has ni nodes (1..ni), flow3d.dat adds ghost at index 0
so array size = (ni+1, nj+1, nk+1) per variable.

Usage:
  python3 interp_flow3d.py coarse_mesh3d.x coarse_flow3d.dat \
                          fine_mesh3d.x   fine_flow3d.dat
"""
import sys, struct
import numpy as np

def read_plot3d_mesh(path):
    """PLOT3D single-block unformatted: nblock, dims, x, y, z."""
    with open(path, 'rb') as f:
        # Record 1: nblock
        reclen = struct.unpack('<I', f.read(4))[0]
        nblock = struct.unpack('<i', f.read(reclen))[0]
        f.read(4)
        # Record 2: dims
        reclen = struct.unpack('<I', f.read(4))[0]
        dims = struct.unpack('<iii', f.read(reclen))
        f.read(4)
        ni, nj, nk = dims
        n = ni * nj * nk
        # Record 3: x
        reclen = struct.unpack('<I', f.read(4))[0]
        x = np.frombuffer(f.read(reclen), dtype='<f8').reshape(nk, nj, ni)
        f.read(4)
        # Record 4: y
        reclen = struct.unpack('<I', f.read(4))[0]
        y = np.frombuffer(f.read(reclen), dtype='<f8').reshape(nk, nj, ni)
        f.read(4)
        # Record 5: z
        reclen = struct.unpack('<I', f.read(4))[0]
        z = np.frombuffer(f.read(reclen), dtype='<f8').reshape(nk, nj, ni)
        f.read(4)
    return ni, nj, nk, x, y, z

def read_flow3d_dat(path, ni, nj, nk):
    """Read flow3d.dat: one record, U(0:nx, 0:ny, 0:nz, 5) where nx=ni.
    Array dims = (ni+1, nj+1, nk+1) per variable.
    Layout: ((((U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=1,5)
    """
    nx_dat = ni   # nx in Fortran = ni from PLOT3D
    ny_dat = nj
    nz_dat = nk
    # Total doubles: (nx+1)*(ny+1)*(nz+1)*5
    total = (nx_dat+1) * (ny_dat+1) * (nz_dat+1) * 5
    with open(path, 'rb') as f:
        reclen = struct.unpack('<I', f.read(4))[0]
        data = np.frombuffer(f.read(reclen), dtype='<f8')
        f.read(4)
    assert data.size == total, f"Expected {total} doubles, got {data.size}"
    # Reshape: m1 outermost (slowest), then k, j, i innermost (fastest)
    # C-order reshape(5, nz+1, ny+1, nx+1) gives i fastest. OK.
    arr = data.reshape(5, nz_dat+1, ny_dat+1, nx_dat+1)
    return arr  # arr[m1, k, j, i], i=0..nx

def write_flow3d_dat(path, arr, ni, nj, nk):
    """Write flow3d.dat: one record U(0:nx, 0:ny, 0:nz, 5).
    arr shape: (5, nz+1, ny+1, nx+1)
    """
    # Flatten in Fortran order: i fastest, j next, k next, m1 outermost (slowest)
    # arr shape is (5, nk+1, nj+1, ni+1) = (m1, k, j, i)
    # C-order ravel gives last axis (i) fastest, first axis (m1) slowest. Correct!
    flat = arr.ravel(order='C')
    flat_bytes = flat.astype('<f8').tobytes()
    reclen = len(flat_bytes)
    with open(path, 'wb') as f:
        f.write(struct.pack('<I', reclen))
        f.write(flat_bytes)
        f.write(struct.pack('<I', reclen))
    print(f"  Wrote {path}: {reclen+8} bytes")

def interp_2d(coarse_x, coarse_y, fine_x, fine_y, field):
    """Bilinear interp from coarse (nj_c, ni_c) to fine (nj_f, ni_f).
    Uses pure numpy (no scipy needed).
    """
    nj_c, ni_c = field.shape
    # Find x-index in coarse grid for each fine x
    ix = np.searchsorted(coarse_x, fine_x, side='right') - 1
    ix = np.clip(ix, 0, ni_c - 2)
    # Find y-index in coarse grid for each fine y
    jy = np.searchsorted(coarse_y, fine_y, side='right') - 1
    jy = np.clip(jy, 0, nj_c - 2)
    # Weights
    xL = coarse_x[ix]; xR = coarse_x[ix + 1]
    yL = coarse_y[jy]; yR = coarse_y[jy + 1]
    wx = np.where(xR - xL > 1e-14, (fine_x - xL) / (xR - xL), 0.5)
    wy = np.where(yR - yL > 1e-14, (fine_y - yL) / (yR - yL), 0.5)
    wx2 = wx[None, :]
    wy2 = wy[:, None]
    # 4 neighbors
    F00 = field[np.ix_(jy, ix)]
    F10 = field[np.ix_(jy, ix + 1)]
    F01 = field[np.ix_(jy + 1, ix)]
    F11 = field[np.ix_(jy + 1, ix + 1)]
    # Bilinear
    result = (1 - wy2) * ((1 - wx2) * F00 + wx2 * F10) + \
             wy2 * ((1 - wx2) * F01 + wx2 * F11)
    return result

def main():
    coarse_mesh = sys.argv[1]
    coarse_flow = sys.argv[2]
    fine_mesh   = sys.argv[3]
    fine_flow   = sys.argv[4]

    ni_c, nj_c, nk_c, xc, yc, zc = read_plot3d_mesh(coarse_mesh)
    print(f"Coarse mesh: {ni_c}x{nj_c}x{nk_c} PLOT3D nodes")

    ni_f, nj_f, nk_f, xf, yf, zf = read_plot3d_mesh(fine_mesh)
    print(f"Fine mesh: {ni_f}x{nj_f}x{nk_f} PLOT3D nodes")

    coarse = read_flow3d_dat(coarse_flow, ni_c, nj_c, nk_c)
    print(f"Coarse flow3d.dat shape: {coarse.shape}")

    # PLOT3D mesh coords: xc[k, j, i] with i=0..ni-1, j=0..nj-1, k=0..nk-1
    # flow3d.dat has ghost at index 0, real nodes at 1..ni
    # Extract real node coordinates (1D from k=0 layer)
    coarse_x_1d = xc[0, 0, :]  # (ni_c,)
    coarse_y_1d = yc[0, :, 0]  # (nj_c,)
    fine_x_1d   = xf[0, 0, :]  # (ni_f,)
    fine_y_1d   = yf[0, :, 0]  # (nj_f,)

    # flow3d.dat interior nodes: i=1..ni, j=1..nj (real mesh nodes)
    # Ghost at i=0, j=0 (extrapolated boundary)
    # Build fine array: (5, nk_f+1, nj_f+1, ni_f+1)
    fine = np.zeros((5, nk_f+1, nj_f+1, ni_f+1), dtype='<f8')

    for m1 in range(5):
        for k in range(nk_f+1):
            # For k=0 (ghost) or k=nk_f, use nearest real k
            k_use = min(max(k, 1), nk_c)  # map to coarse real k
            # Extract coarse real nodes: i=1..ni_c, j=1..nj_c
            field_c = coarse[m1, k_use, 1:nj_c+1, 1:ni_c+1]  # (nj_c, ni_c)
            # Interpolate to fine real nodes
            field_f = interp_2d(coarse_x_1d, coarse_y_1d,
                                 fine_x_1d, fine_y_1d, field_c)
            # Fill fine real nodes: i=1..ni_f, j=1..nj_f
            fine[m1, k, 1:nj_f+1, 1:ni_f+1] = field_f

    # Fill ghost at i=0, j=0 by copying nearest real node
    fine[:, :, 0, :] = fine[:, :, 1, :]       # i=0 ghost = i=1
    fine[:, :, :, 0] = fine[:, :, :, 1]       # j=0 ghost = j=1

    # Apply hard boundary conditions at real boundary nodes
    # rho = 1 everywhere
    fine[0, :, :, :] = 1.0
    # u: bottom(j=1)=0, left(i=1)=0, right(i=ni_f)=0, top(j=nj_f)=1(lid)
    fine[1, :, 1, :] = 0.0       # bottom wall
    fine[1, :, :, 1] = 0.0       # left wall
    fine[1, :, :, ni_f] = 0.0    # right wall
    fine[1, :, nj_f, :] = 1.0    # top wall (lid)
    # v, w: all walls = 0
    for v in [2, 3]:
        fine[v, :, 1, :] = 0.0
        fine[v, :, :, 1] = 0.0
        fine[v, :, :, ni_f] = 0.0
        fine[v, :, nj_f, :] = 0.0
    # T: 300 K everywhere
    fine[4, :, :, :] = 300.0

    write_flow3d_dat(fine_flow, fine, ni_f, nj_f, nk_f)
    print("Done.")

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Generate O-grid mesh for a 2D half-annulus (semicircle ring).
i: radial direction (inner to outer)
j: circumferential direction (0 to pi, flat ends)
k: 2 layers (quasi-2D)

Mesh: (nr+1) x (ntheta+1) x 2 nodes
Output: Mesh3d.x (Fortran unformatted binary, PLOT3D)
"""
import struct
import array
import math

# Parameters
R_inner = 0.1       # inner radius (m)
R_outer = 0.2       # outer radius (m)
nr = 40             # radial cells
ntheta = 80         # circumferential cells (0 to pi)
nk = 2              # 2 nodes in z (quasi-2D)

nx = nr + 1
ny = ntheta + 1
nz = nk

# Generate coordinates in PLOT3D order (i-fastest, j-medium, k-slowest)
x_coords = []
y_coords = []
z_coords = []
dz = 0.001
for k in range(nz):
    zk = k * dz
    for j in range(ny):
        theta = math.pi * j / ntheta    # 0 to pi (half-annulus)
        for i in range(nx):
            r = R_inner + (R_outer - R_inner) * i / nr
            x_coords.append(r * math.cos(theta))
            y_coords.append(r * math.sin(theta))
            z_coords.append(zk)

def write_fortran_record(f, data, fmt):
    if isinstance(data, (list, tuple)):
        raw = struct.pack(fmt, *data)
    elif isinstance(data, array.array):
        raw = data.tobytes()
    else:
        raw = data
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))

with open("Mesh3d.x", "wb") as f:
    # Record 1: NB (number of blocks)
    write_fortran_record(f, [1], 'i')
    # Record 2: NI, NJ, NK for each block
    write_fortran_record(f, [nx, ny, nz], 'iii')
    # Record 3: x coordinates
    write_fortran_record(f, array.array('d', x_coords), '')
    # Record 4: y coordinates
    write_fortran_record(f, array.array('d', y_coords), '')
    # Record 5: z coordinates
    write_fortran_record(f, array.array('d', z_coords), '')

print(f"Mesh3d.x generated: {nx} x {ny} x {nz} nodes (half-annulus)")
print(f"  R_inner={R_inner}m, R_outer={R_outer}m")
print(f"  nr={nr} radial cells, ntheta={ntheta} circumferential cells")
print(f"  Total cells: {nr * ntheta * (nk-1)}")
print(f"  Theta range: 0 to pi (half-annulus)")

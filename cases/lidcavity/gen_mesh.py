#!/usr/bin/env python3
"""Generate a lid-driven cavity mesh (unit square, quasi-2D).
nx=41 x ny=41 x nz=2 nodes  ->  40 x 40 x 1 cells.
Output: Mesh3d.x (Fortran unformatted, PLOT3D single block).
"""
import struct
import array

nx, ny, nz = 41, 41, 2
Lx, Ly = 1.0, 1.0
dz = 1.0

x_coords, y_coords, z_coords = [], [], []
for k in range(nz):
    zk = k * dz
    for j in range(ny):
        yj = Ly * j / (ny - 1)
        for i in range(nx):
            x_coords.append(Lx * i / (nx - 1))
            y_coords.append(yj)
            z_coords.append(zk)

def write_fortran_record(f, data, fmt):
    if isinstance(data, array.array):
        raw = data.tobytes()
    else:
        raw = struct.pack(fmt, *data)
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))

with open("Mesh3d.x", "wb") as f:
    write_fortran_record(f, [1], 'i')
    write_fortran_record(f, [nx, ny, nz], 'iii')
    write_fortran_record(f, array.array('d', x_coords), '')
    write_fortran_record(f, array.array('d', y_coords), '')
    write_fortran_record(f, array.array('d', z_coords), '')

print(f"Mesh3d.x generated: {nx} x {ny} x {nz} nodes (lid-driven cavity)")

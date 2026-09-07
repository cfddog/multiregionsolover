#!/usr/bin/env python3
"""Generate the 1D transpiration-cooling LTNE test mesh (porous block).

1D bar along x (thickness of the porous wall): nx nodes in x, ny=nz=9 nodes
in the (adiabatic) cross-section.  Length Lx = wall thickness.
Physical model (see README):
  x=0 (i- face, "cold end"): coolant inlet  Tf=300 K, G=2 kg/m2/s; frame
      convective BC h_c=31.4 W/m2/K to 300 K (solid_bc.inp face 1).
  x=L (i+ face, "hot end"): frame heat flux q=2e6 W/m2 (solid_bc face 4),
      fluid outlet (zero gradient).
  j- / j+ / k- / k+ faces: symmetry (adiabatic, impermeable).

Outputs (same style as cases/solid_1d/gen_mesh.py):
  Mesh3d.dat  - formatted PLOT3D
  Mesh3d.x    - binary (unformatted Fortran, gfortran record markers)
  bc3d.inp    - boundary file: code 5 inflow (x-), 6 outflow (x+),
                3 symmetry (other four faces)

Usage: python3 gen_case.py [nx] [ny] [nz]   (defaults 401 9 9)
"""
import struct
import array
import sys

nx = int(sys.argv[1]) if len(sys.argv) > 1 else 401
ny = int(sys.argv[2]) if len(sys.argv) > 2 else 9
nz = int(sys.argv[3]) if len(sys.argv) > 3 else 9
Lx = 0.010           # wall thickness 10 mm
Ly = 0.1
Lz = 0.1
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
dz = Lz / (nz - 1)

# Coordinates in PLOT3D order (i fastest, j medium, k slowest)
x_coords, y_coords, z_coords = [], [], []
for k in range(nz):
    zk = k * dz
    for j in range(ny):
        yj = j * dy
        for i in range(nx):
            x_coords.append(i * dx)
            y_coords.append(yj)
            z_coords.append(zk)

with open("Mesh3d.dat", "w") as f:
    f.write("1\n")
    f.write(f"{nx} {ny} {nz}\n")
    for v in x_coords:
        f.write(f"{v:.15e} ")
    f.write("\n")
    for v in y_coords:
        f.write(f"{v:.15e} ")
    f.write("\n")
    for v in z_coords:
        f.write(f"{v:.15e} ")
    f.write("\n")

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
    write_fortran_record(f, [1], 'i')
    write_fortran_record(f, [nx, ny, nz], 'iii')
    write_fortran_record(f, array.array('d', x_coords), '')
    write_fortran_record(f, array.array('d', y_coords), '')
    write_fortran_record(f, array.array('d', z_coords), '')

with open("bc3d.inp", "w") as f:
    f.write("! Gridgen boundary condition file\n")
    f.write("     1\n")
    f.write(f"   {nx:5d}   {ny:5d}   {nz:5d}\n")
    f.write("blk-1\n")
    f.write("6\n")
    # k- (k=1), k+ (k=nz): symmetry
    f.write(f"       1    {nx:5d}       1    {ny:5d}       1       1       3\n")
    f.write(f"       1    {nx:5d}       1    {ny:5d}    {nz:5d}    {nz:5d}       3\n")
    # j- (j=1), j+ (j=ny): symmetry
    f.write(f"       1    {nx:5d}       1       1       1    {nz:5d}       3\n")
    f.write(f"       1    {nx:5d}    {ny:5d}    {ny:5d}       1    {nz:5d}       3\n")
    # i- (x=0): low-speed inlet;  i+ (x=L): outlet
    f.write(f"       1       1       1    {ny:5d}       1    {nz:5d}       5\n")
    f.write(f"    {nx:5d}    {nx:5d}       1    {ny:5d}       1    {nz:5d}       6\n")

print(f"Mesh3d.dat/.x generated: nx={nx}, ny={ny}, nz={nz}, Lx={Lx} m")
print(f"dx = {dx:.3e} m   ({nx-1} cells along wall thickness)")

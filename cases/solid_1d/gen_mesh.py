#!/usr/bin/env python3
"""Generate a 1D solid heat conduction test case mesh.

1D bar: nx=51 nodes (50 cells) in x, ny=9, nz=9 nodes in y/z
Length L=1.0m, cross-section 0.1m x 0.1m

Outputs:
  Mesh3d.dat  - formatted PLOT3D format
  Mesh3d.x    - binary (unformatted Fortran) format

Boundary conditions:
  Left end  (i=1):  isothermal wall, T=300K
  Right end (i=nx): wall (defaults to adiabatic since no solid_bc entry)
  Other faces:      symmetry (adiabatic)

Steady-state analytical solution: T(x) = 300K (uniform)
"""
import struct
import array

nx, ny, nz = 51, 9, 9
Lx = 1.0
Ly = 0.1
Lz = 0.1
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
dz = Lz / (nz - 1)

# Generate coordinates in PLOT3D order (i-fastest, j-medium, k-slowest)
# This matches Fortran read: read(99) (((Ux(i,j,k,1),i=1,nx),j=1,ny),k=1,nz)
x_coords = []
y_coords = []
z_coords = []
for k in range(nz):
    zk = k * dz
    for j in range(ny):
        yj = j * dy
        for i in range(nx):
            x = i * dx
            x_coords.append(x)
            y_coords.append(yj)
            z_coords.append(zk)

# Write formatted Mesh3d.dat
with open("Mesh3d.dat", "w") as f:
    f.write(f"{1}\n")
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

print(f"Mesh3d.dat generated: nx={nx}, ny={ny}, nz={nz}")

# Write binary Mesh3d.x (Fortran unformatted sequential access)
# gfortran record format: 4-byte length, data, 4-byte length
def write_fortran_record(f, data, fmt):
    """Write a Fortran unformatted record with gfortran-style record markers."""
    if isinstance(data, (list, tuple)):
        raw = struct.pack(fmt, *data)
    elif isinstance(data, array.array):
        # array.tobytes() gives the raw binary bytes
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
    data_x = array.array('d', x_coords)
    write_fortran_record(f, data_x, '')
    # Record 4: y coordinates
    data_y = array.array('d', y_coords)
    write_fortran_record(f, data_y, '')
    # Record 5: z coordinates
    data_z = array.array('d', z_coords)
    write_fortran_record(f, data_z, '')

print(f"Mesh3d.x generated: {len(x_coords)} nodes per coordinate")
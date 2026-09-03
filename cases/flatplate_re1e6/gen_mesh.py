"""Regenerate the flat-plate laminar boundary-layer case (low-speed solver).

Domain:  x in [0,1] (120 cells, uniform), y in [0,0.05] (60 cells, geometric,
first-cell ~4.3e-6, r=1.13), z one cell.  BC types (explicit codes):
  i- = 21 low-speed pressure inlet (LS_Inlet_Type=3, LS_P_in ~ 1/2 rho U^2)
  i+ = 22 low-speed pressure outlet
  j- = farfield (4) on x in [0, 6*dx] (blunt leading edge), then wall (2)
  j+ = farfield (4)
  k- / k+ = symmetry (3)
Run from this directory:  python3 gen_mesh.py
"""
import struct
import numpy as np

nx, ny, nz = 121, 61, 2
x = np.linspace(0.0, 1.0, nx)
r = 1.13
d1 = 0.05 * (r - 1.0) / (r**60 - 1.0)          # first cell height ~4.3e-6
y = np.zeros(ny)
for j in range(1, ny):
    y[j] = y[j - 1] + d1 * r ** (j - 1)
y[ny - 1] = 0.05
z = np.array([0.0, 1.0])
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')


def rec_bytes(a):
    b = np.ascontiguousarray(a, dtype='<f8').tobytes()
    return struct.pack('<I', len(b)) + b + struct.pack('<I', len(b))


def rec_ints(v):
    b = np.ascontiguousarray(np.array(v, dtype='<i4')).tobytes()
    return struct.pack('<I', len(b)) + b + struct.pack('<I', len(b))


coords = [X.ravel(order='F'), Y.ravel(order='F'), Z.ravel(order='F')]
# Mesh3d.x: three coordinate records (Mesh_File_Format=2); Mesh3d.dat: one
for fn, nrec in (('Mesh3d.dat', 1), ('Mesh3d.x', 3)):
    with open(fn, 'wb') as f:
        f.write(rec_ints([1]))
        f.write(rec_ints([nx, ny, nz]))
        if nrec == 1:
            f.write(rec_bytes(np.concatenate(coords)))
        else:
            for c in coords:
                f.write(rec_bytes(c))


def F7(v):
    return ''.join('%8d' % t for t in v)


lines = ['     1', '     1', '%8d%8d%8d' % (nx, ny, nz), 'plate', '7']
faces = [
    (1, nx, 1, ny, 1, 1, 3),      # k- symmetry
    (1, nx, 1, ny, 2, 2, 3),      # k+ symmetry
    (1, 1, 1, ny, 1, 2, 21),      # i- pressure inlet
    (nx, nx, 1, ny, 1, 2, 22),    # i+ pressure outlet
    (1, 6, 1, 1, 1, 2, 4),        # j- leading-edge segment: farfield (slip)
    (7, nx, 1, 1, 1, 2, 2),       # j- plate: adiabatic no-slip wall
    (1, nx, ny, ny, 1, 2, 4),     # j+ farfield
]
for fv in faces:
    lines.append(F7(fv))
with open('bc3d.inp', 'w') as f:
    f.write('\n'.join(lines) + '\n')
print('Mesh3d.dat / Mesh3d.x / bc3d.inp written (dy1=%.2e)' % d1)

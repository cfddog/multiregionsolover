#!/usr/bin/env python3
"""Half-cylinder (upper half, symmetry plane y=0) polar O-grid for the 2D
circular cylinder at Re=40.

Mesh topology (ONE block, structured, body-fitted):

      i : radial,      i=1 -> r=R (cylinder wall), i=nx -> r=R_out
      j : azimuthal,   j=1 -> theta=0 (+x axis),   j=ny -> theta=pi (-x axis)
      k : z,           k=1..nz (nz=2 -> one cell thick, quasi-2D)

  only the UPPER half (theta = 0..pi) of the cylinder is meshed; the two
  theta=const. faces (j=1 and j=ny) both lie on the symmetry plane y=0
  (downstream and upstream ray) and carry BC_Symmetry, which halves the
  grid size and (for this symmetric case) also enforces Cl == 0 exactly.

Boundary conditions written to bc3d.inp:

  i=1  face            : BC_Wall      (2)  no-slip cylinder, upper surface
  i=nx face, j<  jsplit: BC_Outflow   (6)  downstream far field, p = LS_P_out
  i=nx face, j>=jsplit : BC_Inflow    (5)  upstream   far field, u = LS_U_in
  j=1 , j=ny faces     : BC_Symmetry  (3)  symmetry plane y = 0
  k=1 , k=nz  faces    : BC_Symmetry  (3)  quasi-2D (w = 0)

  The far-field split is placed exactly at theta=pi/2, where the free stream
  is parallel to the boundary (no mass flux through the split point).

Usage:
  python3 gen_case.py [nx] [ny] [R_out] [beta]

  nx,ny   node numbers (default 121 x 81 -> 120x80 cells)
  R_out   far-field radius in diameters (default 25)
  beta    tanh wall clustering parameter (default 2.5)

Always writes Mesh3d.x (PLOT3D unformatted, "separate x/y/z records"
layout used by the rest of this repository), bc3d.inp and material.in
in the current directory.
"""
import struct, array, math, sys

D     = 1.0                 # cylinder diameter  (length unit = D)
R_in  = 0.5 * D             # cylinder radius
nx    = int(sys.argv[1]) if len(sys.argv) > 1 else 121
ny    = int(sys.argv[2]) if len(sys.argv) > 2 else 81
R_out = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0 * D
beta  = float(sys.argv[4]) if len(sys.argv) > 4 else 2.5
nz    = 2                   # two k-planes -> one cell in z (quasi 2D)

if ny % 2 == 0:
    sys.exit('ny must be odd so that theta=pi/2 is a grid line')

def r_of(i):                # tanh clustering towards the wall
    xi = i / (nx - 1.0)
    return R_in + (R_out - R_in) * (1.0 - math.tanh(beta * (1.0 - xi)) / math.tanh(beta))

# ---- node coordinates, ordering i fastest, then j, then k ------------------
x = array.array('d'); y = array.array('d'); z = array.array('d')
for k in range(nz):
    for j in range(ny):
        th = math.pi * j / (ny - 1.0)
        for i in range(nx):
            r = r_of(i)
            x.append(r * math.cos(th)); y.append(r * math.sin(th)); z.append(float(k))

def wr(f, raw):
    f.write(struct.pack('I', len(raw))); f.write(raw); f.write(struct.pack('I', len(raw)))

with open('Mesh3d.x', 'wb') as f:
    wr(f, struct.pack('i', 1))                    # total blocks
    wr(f, struct.pack('3i', nx, ny, nz))          # block 1 node numbers
    for a in (x, y, z):                           # three records per block
        wr(f, a.tobytes())

# ---- boundary conditions --------------------------------------------------
js = (ny + 1) // 2            # node index of theta = pi/2
jsplit = js                   # outlet : j = 1 .. js ; inlet : j = js .. ny
faces = [
    (  1,   1,     1,   ny, 1,   nz,  2),   # i=1  : cylinder wall (no-slip)
    ( nx,  nx,     1,   js, 1,   nz,  6),   # i=nx : downstream far field (p_out)
    ( nx,  nx,   js,   ny, 1,   nz,  5),   # i=nx : upstream  far field (u_in)
    (  1,  nx,     1,    1, 1,   nz,  3),   # j=1  : symmetry plane y=0 (downstream ray)
    (  1,  nx,   ny,   ny, 1,   nz,  3),   # j=ny : symmetry plane y=0 (upstream ray)
    (  1,  nx,     1,   ny, 1,    1,  3),   # k=1  : symmetry (quasi 2D)
    (  1,  nx,     1,   ny, nz, nz,  3),   # k=nz : symmetry (quasi 2D)
]
lines = ['! half-cylinder polar grid, upper half with symmetry plane y=0',
         '%6d' % 1,
         '%6d%6d%6d' % (nx, ny, nz),
         'blk-1',
         '%6d' % len(faces)]
for fc in faces:
    lines.append(''.join('%8d' % v for v in fc))
open('bc3d.inp', 'w').write('\n'.join(lines) + '\n')

# ---- material : one BLOCK_LOWSPEED block ---------------------------------
open('material.in', 'w').write('1\n2\n1.0 1.0 0.0\n')

dr1 = r_of(1) - R_in
print('half cylinder %dx%dx%d  (cells %dx%dx%d)  R_in=%.3f R_out=%.1f beta=%.2f'
      % (nx, ny, nz, nx - 1, ny - 1, nz - 1, R_in, R_out, beta))
print('  first cell height dr1=%.5f D (d+ = %.2f at Re=40)  dtheta=%.3f deg, arc=%.5f D'
      % (dr1, 0.5 * dr1 * 40.0, 180.0 / (ny - 1), R_in * math.pi / (ny - 1)))

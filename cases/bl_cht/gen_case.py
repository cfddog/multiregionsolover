#!/usr/bin/env python3
"""Generate the incompressible boundary-layer + solid-plate conjugate case.

Geometry (SI, L = 1 m plate):
  blk-1 (LOWSPEED fluid):  x in [0,4] (plate on x in [0,1], wake region up to 4)
                           y in [0,2], free stream U=1 m/s, T=300 K
  blk-2 (SOLID plate):     x in [0,1], y in [-0.1,0], bottom y=-0.1 at 800 K
Conjugate interface: blk-1 j- face over x in [0,1]  <->  blk-2 j+ face.
Re_L = rho U L / mu = 1000  (rho=1, U=1, mu=1e-3);  Pr = mu Cp / k = 0.7.
Writes Mesh3d.x (multi-block Fortran unformatted) plus bc3d.inp,
bc3d_interface.inp, material.in, solid_bc.inp.
"""
import struct, array, math, os

def record(f, data, fmt):
    if isinstance(data, array.array):
        raw = data.tobytes()
    else:
        raw = struct.pack(fmt, *data)
    f.write(struct.pack('I', len(raw))); f.write(raw); f.write(struct.pack('I', len(raw)))

def tangent_hyperbolic(n, L, beta):
    if n <= 1 or beta <= 1e-12:
        return [L*i/max(n-1,1) for i in range(n)]
    c = []
    for i in range(n):
        xi = i/(n-1)
        c.append(L*(0.5+0.5*math.tanh(beta*(xi-0.5))/math.tanh(beta/2.0)))
    return c

# ---------- coordinates ----------
# fluid x: plate part [0,1] dx=0.01 (101 nodes, 1:1 with solid), wake [1,4]
# stretched, 50 cells
dxp = 0.01
xs_f = [dxp*i for i in range(101)]
rest_n = 50
x = 1.0
for k in range(rest_n):
    x += 0.02 + k*0.001633
    xs_f.append(round(x, 12))
nx1 = len(xs_f)              # 151
nz1 = 2

# y in [0,2]: geometric clustering from the wall (dy1, ratio r) up to ~0.93,
# then near-uniform spacing to the top.  ~38 cells inside the thermal BL at
# x=1 (delta_T ~ 0.17 m), first cell 1e-3 m.
dy1, r, n1 = 1.0e-3, 1.07, 62
ys_f = [0.0]
for k in range(n1):
    ys_f.append(ys_f[-1] + dy1*r**k)
L2 = 2.0 - ys_f[-1]
n2 = math.ceil(L2/0.05)
dy2 = L2/n2
for k in range(n2):
    ys_f.append(ys_f[-1] + dy2)
ny1 = len(ys_f)
print(f"fluid y: n={ny1}, h_min={ys_f[1]-ys_f[0]:.2e}, h_max={ys_f[-1]-ys_f[-2]:.3f}")

# solid: same x as plate part, uniform y in [-0.1, 0]
xs_s = [dxp*i for i in range(101)]
nx2 = len(xs_s)
ny2 = 11
nz2 = 2
ys_s = [-0.1 + 0.01*j for j in range(ny2)]   # 10 cells of 0.01

def coords(xs, ys, nz, dz=1.0):
    xc, yc, zc = [], [], []
    for k in range(nz):
        for j in range(len(ys)):
            for i in range(len(xs)):
                xc.append(xs[i]); yc.append(ys[j]); zc.append(k*dz)
    return xc, yc, zc

with open("Mesh3d.x", "wb") as f:
    record(f, [2], 'i')
    record(f, [nx1, ny1, nz1, nx2, ny2, nz2], 'iiiiii')
    for xs, ys, nz in [(xs_f, ys_f, nz1), (xs_s, ys_s, nz2)]:
        xc, yc, zc = coords(xs, ys, nz)
        record(f, array.array('d', xc), '')
        record(f, array.array('d', yc), '')
        record(f, array.array('d', zc), '')
print(f"Mesh3d.x: blk1 {nx1}x{ny1}x{nz1}, blk2 {nx2}x{ny2}x{nz2}")

# ---------- bc3d.inp (physical faces; interface segment listed as wall(2),
# overwritten at runtime by the coupling) ----------
def F(i1,i2,j1,j2,k1,k2,bc): return f"{i1:8d}{i2:8d}{j1:8d}{j2:8d}{k1:8d}{k2:8d}{bc:8d}"
lines = ["", "     2"]
lines += [f"{nx1:8d}{ny1:8d}{nz1:8d}", "blk-1", "7"]
lines += [F(1,1,1,ny1,1,2,5),          # i- inflow
          F(nx1,nx1,1,ny1,1,2,6),      # i+ outflow
          F(1,nx2,1,1,1,2,2),          # j- plate part (interface, wall placeholder)
          F(nx2+1,nx1,1,1,1,2,3),      # j- wake part, symmetry
          F(1,nx1,ny1,ny1,1,2,3),      # j+ symmetry
          F(1,nx1,1,ny1,1,1,3),        # k-
          F(1,nx1,1,ny1,2,2,3)]        # k+
lines += [f"{nx2:8d}{ny2:8d}{nz2:8d}", "blk-2", "6"]
lines += [F(1,1,1,ny2,1,2,2),          # i- end face (adiabatic)
          F(nx2,nx2,1,ny2,1,2,2),      # i+ end face (adiabatic)
          F(1,nx2,1,1,1,2,2),          # j- bottom 800 K
          F(1,nx2,ny2,ny2,1,2,2),      # j+ interface (wall placeholder)
          F(1,nx2,1,ny2,1,1,2),        # k-
          F(1,nx2,1,ny2,2,2,2)]        # k+
open("bc3d.inp", "w").write("\n".join(lines) + "\n")

# ---------- bc3d_interface.inp (same face layout per block; the interface
# faces get negative bc, followed by a pairing line (range on the OTHER block
# + its block number). Order of faces must match bc3d.inp per block. ----------
lines = ["", "     2"]
lines += [f"{nx1:8d}{ny1:8d}{nz1:8d}", "blk-1", "7"]
lines += [F(1,1,1,ny1,1,2,5),
          F(nx1,nx1,1,ny1,1,2,6),
          F(1,nx2,1,1,1,2,-1),         # interface j- over plate
          F(1,nx2,ny2,ny2,1,2,2),      # pairing line: solid j+ range, block 2
          F(nx2+1,nx1,1,1,1,2,3),
          F(1,nx1,ny1,ny1,1,2,3),
          F(1,nx1,1,ny1,1,1,3),
          F(1,nx1,1,ny1,2,2,3)]
lines += [f"{nx2:8d}{ny2:8d}{nz2:8d}", "blk-2", "6"]
lines += [F(1,1,1,ny2,1,2,2),
          F(nx2,nx2,1,ny2,1,2,2),
          F(1,nx2,1,1,1,2,2),
          F(1,nx2,ny2,ny2,1,2,-1),     # interface j+
          F(1,nx2,1,1,1,2,1),          # pairing line: fluid j- range, block 1
          F(1,nx2,1,ny2,1,1,2),
          F(1,nx2,1,ny2,2,2,2)]
open("bc3d_interface.inp", "w").write("\n".join(lines) + "\n")

# ---------- material.in ----------
open("material.in", "w").write("2\n2 1\n1.0 1000.0 1.4286\n7800.0 460.0 50.0\n")

# ---------- solid_bc.inp : solid block 2, face 2 (j-) at 800 K ----------
open("solid_bc.inp", "w").write("1\n2\n1\n2 800.0 0.0\n")

print("input files written")

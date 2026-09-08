#!/usr/bin/env python3
"""Validate solid_annulus steady radial conduction vs analytic log profile.

Geometry: half-annulus O-grid, R_in = 0.1 m, R_out = 0.2 m, quasi-2D.
BCs (solid_bc.inp): i- (inner, face 1) heat flux q=800 W/m2 into domain;
                    i+ (outer, face 4) isothermal 300 K.
Material (material.in): k = 15.1 W/mK.
Steady radial analytic:  T(r) = T_out + (q*R_in/k) * ln(R_out/r)
Reads Tecplot output Ts_block_1.dat (I=nr=40 J=ntheta=80 K=nz-1=1, POINT).
Data order: i fastest -> row = i + nr*j (k fixed at 0).
"""
import numpy as np

Q, K = 800.0, 15.1
RI, RO, T_OUT = 0.1, 0.2, 300.0
NR, NTH = 40, 80

d = np.loadtxt("Ts_block_1.dat", skiprows=3)
x, y, T = d[:, 0], d[:, 1], d[:, 3]
r = np.hypot(x, y)

i_idx = np.arange(d.shape[0]) % NR          # radial layer (i), i fastest
T_an = T_OUT + (Q * RI / K) * np.log(RO / np.clip(r, 1e-12, None))

# radial-layer means over theta (per i)
r_i = np.array([r[i_idx == i].mean() for i in range(NR)])
T_i = np.array([T[i_idx == i].mean() for i in range(NR)])
T_an_i = T_OUT + (Q * RI / K) * np.log(RO / np.clip(r_i, 1e-12, None))

print(f"cells: {r.size}   r in [{r.min():.5f}, {r.max():.5f}]")
print(f"layer1 (inner) T = {T_i[0]:.6f}  analytic {T_an_i[0]:.6f}  (r={r_i[0]:.5f})")
print(f"layer20 (mid) T = {T_i[19]:.6f}  analytic {T_an_i[19]:.6f}  (r={r_i[19]:.5f})")
print(f"layer40 (outer)T = {T_i[-1]:.6f}  analytic {T_an_i[-1]:.6f}  (r={r_i[-1]:.5f})")
print(f"radial RMS err (40 layer means) = {np.sqrt(np.mean((T_i-T_an_i)**2)):.3e}")
print(f"per-cell RMS err = {np.sqrt(np.mean((T-T_an)**2)):.3e}   max|err| = {np.max(np.abs(T-T_an)):.3e}")

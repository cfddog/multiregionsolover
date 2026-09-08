#!/usr/bin/env python3
"""Validate solid_1d steady conduction vs analytic linear profile.

BCs (solid_bc.inp): i- (x=0) isothermal 300 K; i+ (x=L=1 m) heat flux q=800 W/m2.
Material (material.in): k = 15.1 W/mK.
Steady 1D analytic:  T(x) = 300 + (q/k)*x
Reads Tecplot output Ts_block_1.dat (I=nx-1 J=ny-1 K=nz-1, POINT datapacking).
Data order: i fastest (I per row), then j, then k.
"""
import numpy as np

Q, K = 800.0, 15.1          # W/m2, W/mK
T_LEFT = 300.0

d = np.loadtxt("Ts_block_1.dat", skiprows=3)
x, T = d[:, 0], d[:, 3]

# group over (j,k) duplicates: unique x -> mean T per x (1D profile)
xu, inv = np.unique(np.round(x, 12), return_inverse=True)
T_x = np.zeros_like(xu)
for k in range(xu.size):
    T_x[k] = T[inv == k].mean()

T_an = T_LEFT + (Q / K) * xu
err = T_x - T_an
rms = np.sqrt(np.mean(err ** 2))
linf = np.max(np.abs(err))
iL, iR, iM = 0, xu.size - 1, np.argmin(np.abs(xu - 0.5))
print(f"profile pts (distinct x): {xu.size}   x in [{xu.min():.4f}, {xu.max():.4f}]")
print(f"T(x=0.01) = {T_x[iL]:.6f}  analytic {T_an[iL]:.6f}")
print(f"T(x=0.50) = {T_x[iM]:.6f}  analytic {T_an[iM]:.6f}")
print(f"T(x=0.99) = {T_x[iR]:.6f}  analytic {T_an[iR]:.6f}")
print(f"RMS err    = {rms:.3e}   max|err| = {linf:.3e}")
print(f"rel. max   = {linf/(T_an[iR]-T_an[iL]):.2e}")

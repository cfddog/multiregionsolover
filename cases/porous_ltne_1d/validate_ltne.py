#!/usr/bin/env python3
"""Validate porous_ltne_1d (1D transpiration cooling, LTNE two-temperature)
against the closed-form analytical solution of the two-temperature model.

Model (steady, 1D, y in [0,L], cold inlet at y=0, hot face at y=L):
  fluid:  G*cp * dTf/dy       = hv*(Ts - Tf),        Tf(0)=Tin
  solid:  (1-eps)*ks * d2Ts/dy2 + hv*(Tf - Ts) = 0
          cold face: kse*dTs/dy(0)  = hc*(Ts(0)-Tinf)
          hot  face: kse*dTs/dy(L)  = q        (per total wall area)

Closed form: Ts''' + A Ts'' - B Ts' = 0,
  A = hv/(G*cp),  B = hv/kse,  kse=(1-eps)*ks,  beta=kse/(G*cp)
  Ts  = C0 + C1*exp(l1*y) + C2*exp(l2*y)
  Tf  = Ts - (kse/hv)*Ts'' = C0 + beta*(C1*l1*exp(l1*y)+C2*l2*exp(l2*y))
  l1,2 = (-A +/- sqrt(A^2+4B))/2
  BCs: Tf(0)=Tin; kse*(C1*l1+C2*l2)=hc*(C0+C1+C2-Tinf); kse*Ts'(L)=q

Numerical profiles: Ts from Ts_block_1.dat (cell-centred), Tf from
flow3d_block_1.vtk scalar 'temperature' (cell-centred); both share the same
cell-centre x coordinates (i fastest; 1 cell in j,k in the quasi-1D grid).
"""
import sys
import os
import numpy as np

# ------------------------- parameters (case) --------------------------------
L   = 0.010
eps = 0.30
ks  = 13.4
G   = 2.0            # kg/m2/s
cp  = 1007.0
hv  = 2.0e6          # W/m3/K (Wakao-Kaguei, dp=1mm, air 300K -> ~2.0e6)
q   = 2.0e6          # W/m2  hot-face heat flux (into the solid frame)
hc  = 31.4           # W/m2/K cold-face convective coefficient
Tinf = 300.0
Tin  = 300.0
kse = (1.0 - eps) * ks

# ------------------------- analytical solution ------------------------------
def analytic(y):
    A = hv / (G * cp)
    B = hv / kse
    beta = kse / (G * cp)
    D = np.sqrt(A * A + 4.0 * B)
    l1 = (-A + D) / 2.0
    l2 = (-A - D) / 2.0
    e1, e2 = np.exp(l1 * y), np.exp(l2 * y)
    e1L, e2L = np.exp(l1 * L), np.exp(l2 * L)
    s1, s2 = l1, l2
    # linear system for [C0, C1, C2]
    M = np.array([
        [1.0,          beta * l1,     beta * l2],
        [-hc,          kse * l1 - hc, kse * l2 - hc],
        [0.0,          kse * l1 * e1L, kse * l2 * e2L],
    ])
    rhs = np.array([Tin, -hc * Tinf, q])
    C0, C1, C2 = np.linalg.solve(M, rhs)
    Ts = C0 + C1 * e1 + C2 * e2
    Tf = C0 + beta * (C1 * l1 * e1 + C2 * l2 * e2)
    return Ts, Tf

# ------------------------- read numerical fields ----------------------------
def read_Ts(path="Ts_block_1.dat"):
    d = np.loadtxt(path, skiprows=3)
    return d[:, 0], d[:, 3]          # cell-centre x, Ts

def read_flow_t(path="flow3d_block_1.vtk"):
    lines = open(path).read().splitlines()
    vals = []
    i = 0
    # jump to CELL_DATA scalar temperature
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("SCALARS") and lines[i].split()[1] == "temperature":
            i += 2
            while i < len(lines) and not lines[i].strip().startswith(
                    ("SCALARS", "VECTORS", "POINT_DATA", "CELL_DATA", "FIELD")):
                vals += lines[i].split()
                i += 1
            break
        i += 1
    return np.array(vals, float)

base = os.path.dirname(os.path.abspath(__file__))
xts, Ts_num = read_Ts(os.path.join(base, "Ts_block_1.dat"))
Tf_num = read_flow_t(os.path.join(base, "flow3d_block_1.vtk"))
# quasi-1D grid: cell centres coincide for Ts and Tf (same ordering)
x = xts
Tf_num = Tf_num[:x.size]

Ts_an, Tf_an = analytic(x)

# ------------------------- report -------------------------------------------
print(f"case: L={L} m, eps={eps}, ks={ks}, hv={hv:.3e} W/m3K, G={G}, cp={cp}")
print(f"kse=(1-eps)ks = {kse:.4f} W/mK   G*cp = {G*cp:.1f} W/m2K")
eTf = Tf_num - Tf_an
eTs = Ts_num - Ts_an
i0 = np.argmin(np.abs(x - 0.0)); iL = np.argmin(np.abs(x - L))
ic = np.argmin(np.abs(x - L / 2.0))
print(f"{'station':>10} {'Tf_num':>12} {'Tf_an':>12} {'Ts_num':>12} {'Ts_an':>12} {'Ts-Tf_num':>10}")
for nm, ii in [("y=0", i0), ("y=L/4", np.argmin(np.abs(x-L/4))), ("y=L/2", ic),
               ("y=3L/4", np.argmin(np.abs(x-3*L/4))), ("y=L", iL)]:
    print(f"{nm:>10} {Tf_num[ii]:12.4f} {Tf_an[ii]:12.4f} {Ts_num[ii]:12.4f} "
          f"{Ts_an[ii]:12.4f} {Ts_num[ii]-Tf_num[ii]:10.4f}")
print(f"Tf RMS err  = {np.sqrt(np.mean(eTf**2)):.4e} K   max|err| = {np.max(np.abs(eTf)):.4e} K")
print(f"Ts RMS err  = {np.sqrt(np.mean(eTs**2)):.4e} K   max|err| = {np.max(np.abs(eTs)):.4e} K")
dTf_an = np.max(Tf_an) - np.min(Tf_an)
print(f"rel err (Tf, vs dTf) = {np.max(np.abs(eTf))/dTf_an:.3e}")

# energy balance of the numerical field (per unit wall area):
#  in  = q ; out = G*cp*(Tf(y=L)-Tin) + hc*(Ts(y=0)-Tinf)  (fluid axial cond.=0)
TfL_num = Tf_num[iL]; Ts0_num = Ts_num[i0]
out_num = G * cp * (TfL_num - Tin) + hc * (Ts0_num - Tinf)
print(f"energy balance: q={q:.3e}; out=Gcp dTf+ hc dTs ={out_num:.3e}; rel gap={(q-out_num)/q:.3e}")

# ------------------------- plot (optional) ----------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(x * 1e3, Ts_an, "-", color="tab:red", lw=1.6, label="Ts analytic")
    ax.plot(x * 1e3, Tf_an, "-", color="tab:blue", lw=1.6, label="Tf analytic")
    ax.plot(x * 1e3, Ts_num, "o", ms=2.6, color="tab:red", alpha=0.55, label="Ts CFD")
    ax.plot(x * 1e3, Tf_num, "s", ms=2.6, color="tab:blue", alpha=0.55, label="Tf CFD")
    ax.set_xlabel("y  (mm)")
    ax.set_ylabel("T  (K)")
    ax.set_title("porous_ltne_1d: transpiration cooling LTNE vs analytical")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(base, "validation_ltne.png"), dpi=150)
    print("saved validation_ltne.png")
except Exception as exc:
    print("plot not written:", exc)

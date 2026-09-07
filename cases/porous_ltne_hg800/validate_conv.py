#!/usr/bin/env python3
"""Porous LTNE 1D with a convective (Newton) hot face - analytical validation.

Same two-temperature model as validate_ltne.py, but the hot-end heat load is
provided by the high-speed flow through a convective BC:
    kse*dTs/dy(L) = hg*(T_rec - Ts(L))
where hg and T_rec are calibrated from cases/high_low_fluid_800K (Step-1).
Cold end unchanged: Tf(0)=Tin, kse*dTs/dy(0)=hc*(Ts(0)-Tinf).
"""
import os
import numpy as np

L    = 0.010
eps  = 0.30
ks   = 13.4
G    = 2.0
cp   = 1007.0
hv   = 2.0e6
hc   = 31.4
Tinf = 300.0
Tin  = 300.0
HG   = 1.48e3       # W/m2/K  (from high_low_fluid_800K interface)
T_REC = 2240.0      # K        (Ma=3, T_inf=800 K total temperature)
kse = (1.0 - eps) * ks

def analytic(y):
    A = hv / (G * cp); B = hv / kse; beta = kse / (G * cp)
    D = np.sqrt(A * A + 4.0 * B)
    l1 = (-A + D) / 2.0; l2 = (-A - D) / 2.0
    e1, e2 = np.exp(l1 * y), np.exp(l2 * y)
    e1L, e2L = np.exp(l1 * L), np.exp(l2 * L)
    M = np.array([
        [1.0, beta * l1, beta * l2],
        [-hc, kse * l1 - hc, kse * l2 - hc],
        [HG,  (HG + kse * l1) * e1L, (HG + kse * l2) * e2L],
    ])
    rhs = np.array([Tin, -hc * Tinf, HG * T_REC])
    C0, C1, C2 = np.linalg.solve(M, rhs)
    Ts = C0 + C1 * e1 + C2 * e2
    Tf = C0 + beta * (C1 * l1 * e1 + C2 * l2 * e2)
    return Ts, Tf

base = os.path.dirname(os.path.abspath(__file__))
d = np.loadtxt(os.path.join(base, "Ts_block_1.dat"), skiprows=3)
x, Ts_num = d[:, 0], d[:, 3]

lines = open(os.path.join(base, "flow3d_block_1.vtk")).read().splitlines()
vals = []; i = 0
while i < len(lines):
    s = lines[i].strip()
    if s.startswith("SCALARS") and lines[i].split()[1] == "temperature":
        i += 2
        while i < len(lines) and not lines[i].strip().startswith(
                ("SCALARS", "VECTORS", "POINT_DATA", "CELL_DATA", "FIELD")):
            vals += lines[i].split(); i += 1
        break
    i += 1
Tf_num = np.array(vals, float)[:x.size]

Ts_an, Tf_an = analytic(x)
eTf, eTs = Tf_num - Tf_an, Ts_num - Ts_an
i0 = np.argmin(np.abs(x - 0)); iL = np.argmin(np.abs(x - L))
print(f"hot-face convective BC: hg={HG:.0f} W/m2K, T_rec={T_REC:.0f} K  (Ma=3, Tinf=800 K)")
print(f"{'station':>8} {'Tf_num':>11} {'Tf_an':>11} {'Ts_num':>11} {'Ts_an':>11} {'Ts-Tf':>9}")
for nm, ii in [("y=0", i0), ("y=L/2", np.argmin(np.abs(x - L/2))), ("y=L", iL)]:
    print(f"{nm:>8} {Tf_num[ii]:11.3f} {Tf_an[ii]:11.3f} {Ts_num[ii]:11.3f} "
          f"{Ts_an[ii]:11.3f} {Ts_num[ii]-Tf_num[ii]:9.3f}")
print(f"Tf RMS {np.sqrt(np.mean(eTf**2)):.4e}  max {np.max(np.abs(eTf)):.4e}  "
      f"Ts RMS {np.sqrt(np.mean(eTs**2)):.4e}  max {np.max(np.abs(eTs)):.4e}")
# heat fluxes: q_hot = hg*(T_rec - Ts(L))  vs  coolant + cold-face loss
q_hot_num = HG * (T_REC - Ts_num[iL])
out = G * cp * (Tf_num[iL] - Tin) + hc * (Ts_num[i0] - Tinf)
print(f"q_hot(hg,Tw) = {q_hot_num:.3e} W/m2 ; coolant+cold loss = {out:.3e} ; rel gap {(q_hot_num-out)/q_hot_num:.3e}")
print(f"wall face Ts(L)={Ts_num[iL]:.1f} K, coolant exit Tf(L)={Tf_num[iL]:.1f} K")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(x*1e3, Ts_an, "-", color="tab:red", lw=1.6, label="Ts analytic")
    ax.plot(x*1e3, Tf_an, "-", color="tab:blue", lw=1.6, label="Tf analytic")
    ax.plot(x*1e3, Ts_num, "o", ms=2.4, color="tab:red", alpha=.55, label="Ts CFD")
    ax.plot(x*1e3, Tf_num, "s", ms=2.4, color="tab:blue", alpha=.55, label="Tf CFD")
    ax.set_xlabel("y (mm)"); ax.set_ylabel("T (K)")
    ax.set_title("porous LTNE wall heated by Ma=3, 800 K flow (convective hot face)")
    ax.legend(); ax.grid(alpha=.3); fig.tight_layout()
    fig.savefig(os.path.join(base, "validation_ltne_hg800.png"), dpi=150)
    print("saved validation_ltne_hg800.png")
except Exception as e:
    print("plot skip:", e)

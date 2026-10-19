#!/usr/bin/env python3
"""Validation of the code-14 SOLID<->SOLID interface (cases/solid_2blk_chain).

Reads Ts_block_1.dat / Ts_block_2.dat written by the solver and compares the
cell-centre profile with the analytic series-conduction solution

    R1 = L1/k1, R2 = L2/k2,  q  = dTtot/(R1+R2)
    T(x) = 400 - q x/k1          0   <= x <= L1
         = 350 + q (L1+L2-x)/k2  L1  <= x <= L1+L2

It also reports the CELL-CENTRE jump across the interface

    dT_cc = T(blk-1 last cell) - T(blk-2 first cell)

whose expected value is q*(dx1/(2 k1) + dx2/(2 k2)) -- NOT zero, because the
interface temperature itself sits half a cell away from either centre.

Usage:  python3 analyze_solid_2blk_chain.py [run_dir]   (default: np1)
Exit code 0 = pass, 1 = fail.
"""
import sys, os, numpy as np

L1, L2 = 0.05, 0.05
K1, K2 = 16.3, 200.0
TL, TR = 400.0, 350.0
NX, DX = 41, L1 / 40.0

q_ex   = (TL - TR) / (L1 / K1 + L2 / K2)
Ti_ex  = TL - q_ex * L1 / K1
dT_ex  = q_ex * (DX / (2 * K1) + DX / (2 * K2))

run = sys.argv[1] if len(sys.argv) > 1 else 'np1'

def profile(fname):
    """Return (x_unique, T_mean_over_yz, T_all) for one block."""
    d = np.loadtxt(fname, skiprows=3)
    x, T = d[:, 0], d[:, 3]
    xu, inv = np.unique(np.round(x, 12), return_inverse=True)
    Tx = np.array([T[inv == k].mean() for k in range(xu.size)])
    return xu, Tx, T

def T_exact(x):
    return np.where(x <= L1 + 1e-12,
                    TL - q_ex * x / K1,
                    TR + q_ex * (L1 + L2 - x) / K2)

ok = True
print('analytic: q = %.6g W/m2   T_iface = %.6f K   dT_cc = %.6f K'
      % (q_ex, Ti_ex, dT_ex))

prof = {}
emin = None
for m in (1, 2):
    f = os.path.join(run, 'Ts_block_%d.dat' % m)
    if not os.path.exists(f):
        print('MISSING', f)
        sys.exit(1)
    xu, Tx, Tall = profile(f)
    prof[m] = (xu, Tx)
    err = Tx - T_exact(xu)
    e = np.abs(err).max()
    emin = e if emin is None else max(emin, e)
    print('blk-%d: %2d columns, max|T-T_exact| = %.3e K (at x=%.6f m)'
          % (m, xu.size, e, xu[np.abs(err).argmax()]))
    if e > 1e-6:
        ok = False
    # physical sanity: monotone decreasing in x inside each block
    if np.any(np.diff(Tx) > 0.0):
        print('  WARNING: non-monotone T(x) in blk-%d' % m)
        ok = False

T1 = prof[1][1][-1]
T2 = prof[2][1][0]
dT = T1 - T2
print('interface: T(blk-1 last centre) = %.6f K, T(blk-2 first centre) = %.6f K'
      % (T1, T2))
print('           dT_cc = %.6f K   (analytic %.6f K, |diff| = %.3e K)'
      % (dT, dT_ex, abs(dT - dT_ex)))
if abs(dT - dT_ex) > 1e-3:
    ok = False

print('profile near the interface (last 6 of blk-1 | first 6 of blk-2):')
for i in range(-6, 0):
    print('   x=%.6f  T=%.6f  T_ex=%.6f  err=% .3e'
          % (prof[1][0][i], prof[1][1][i], T_exact(prof[1][0][i]),
             prof[1][1][i] - T_exact(prof[1][0][i])))
for i in range(0, 6):
    print('   x=%.6f  T=%.6f  T_ex=%.6f  err=% .3e'
          % (prof[2][0][i], prof[2][1][i], T_exact(prof[2][0][i]),
             prof[2][1][i] - T_exact(prof[2][0][i])))

print('RESULT:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)

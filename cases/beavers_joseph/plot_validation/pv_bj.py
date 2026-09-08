#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Beavers-Joseph validation: compare the two-block CFD solution with the
analytical reference (finite-bed Darcy-Brinkman composite, which reduces to
the classical BJ slip profile for a deep bed).

Usage (run from cases/beavers_joseph):
   python3 plot_validation/pv_bj.py [--case .] [--out validation_bj.png]
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import pv_vtk
import bj_analytical as bj

LX = 12.0
H_FL = 1.0
HP = 1.2
EPS, DP, MU = 0.8, 0.5, 0.02


def station(d, x):
    return int(np.argmin(np.abs(d['xc'][0, :] - x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='.', help='case directory')
    ap.add_argument('--out', default='validation_bj.png')
    args = ap.parse_args()

    b1 = pv_vtk.load_block(os.path.join(args.case, 'flow3d_block_1.vtk'))
    b2 = pv_vtk.load_block(os.path.join(args.case, 'flow3d_block_2.vtk'))
    # block1 (fluid): y in [0,h]; block2 (porous): y in [-Hp,0]
    p = bj.params(EPS, DP, MU)
    x = b1['xc'][0]
    dx = x[1] - x[0]
    L = x[-1] - x[0] + dx          # physical streamwise length

    # ------------------------------------------------------------------
    # 1) pressure fit in the developed region -> dp/dx
    xwin = (x > 0.25 * L) & (x < 0.75 * L)
    jf = b1['p'].shape[0] // 2          # mid-fluid row
    jp = b2['p'].shape[0] // 2          # mid-porous row
    coef_f = np.polyfit(x[xwin], b1['p'][jf, xwin], 1)
    coef_p = np.polyfit(x[xwin], b2['p'][jp, xwin], 1)
    P = 0.5 * (coef_f[0] + coef_p[0])   # dp/dx (<0)
    G = -P
    print('pressure fit: dp/dx fluid=%.6e  porous=%.6e  mean=%.6e Pa/m'
          % (coef_f[0], coef_p[0], P))

    # ------------------------------------------------------------------
    # 2) analytical solution with the fitted G
    sol = bj.solve_composite(G, EPS, DP, MU, h=H_FL, Hp=HP)
    sol_d = bj.solve_deepbed(G, EPS, DP, MU, h=H_FL)
    print('uD(an)=%.6e  u_int(an)=%.6e  lambda=%.5f  kappa*Hp=%.2f'
          % (sol['uD'], sol['u_int'], sol['lam'], sol['kappa'] * HP))
    print('deep-bed BJ(alpha=1): u_int=%.6e ; finite-bed vs deep-bed diff=%.4f%%'
          % (sol_d['u_int'], 100 * (sol['u_int'] - sol_d['u_int']) /
             sol_d['u_int']))

    # ------------------------------------------------------------------
    # 3) CFD profile at the central station
    i0 = station(b1, 0.5 * L)
    uf = b1['u'][:, i0]          # fluid cell centres, y ascending from 0
    up = b2['u'][:, i0]          # porous cell centres, y ascending from -Hp
    yf = b1['yc'][:, i0]
    yp = b2['yc'][:, i0]
    yc_all = np.concatenate([yp, yf])
    uc_all = np.concatenate([up, uf])
    ua_all = bj.profile(yc_all, sol)
    umax_a = np.max(np.abs(bj.profile(np.linspace(-HP, H_FL, 2001), sol)))
    n_p = yp.size
    err_f = np.abs(uf - ua_all[n_p:]) / umax_a
    err_p = np.abs(up - ua_all[:n_p]) / umax_a
    err_all = np.abs(uc_all - ua_all) / umax_a
    print('profile station x=%.3f : L2(fluid)=%.3e  L2(porous)=%.3e  '
          'L2(total)=%.3e  Linf(total)=%.3e'
          % (0.5 * L, np.sqrt(np.mean(err_f ** 2)),
             np.sqrt(np.mean(err_p ** 2)), np.sqrt(np.mean(err_all ** 2)),
             np.max(err_all)))

    # ------------------------------------------------------------------
    # 4) interface quantities from the CFD (linear extrapolation, fluid side)
    dy = yf[1] - yf[0]
    u_i = uf[0] - (uf[1] - uf[0]) * 0.5        # extrapolate to y=0
    dui = (uf[1] - uf[0]) / dy                 # du/dy at interface
    mD = (yp > -0.9 * HP) & (yp < -0.25 * HP)
    uD_cfd = np.mean(up[mD])
    alpha_eff = dui * sol['lam'] / (u_i - uD_cfd)
    print('interface: u_i(cfd)=%.6e (an=%.6e, err=%.3f%%)' %
          (u_i, sol['u_int'], 100 * (u_i - sol['u_int']) / sol['u_int']))
    print('Darcy plateau uD: cfd=%.6e  an=%.6e' % (uD_cfd, sol['uD']))
    print('du/dy|0 = %.6e 1/s ;  BJ alpha_eff = lam*du/dy/(u_i-uD) = %.4f '
          '(theory 1.0)' % (dui, alpha_eff))

    # ------------------------------------------------------------------
    # 5) developed-flow check
    stations = [0.3 * L, 0.5 * L, 0.7 * L]
    ya = np.linspace(-HP, H_FL, 1201)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    ax[0, 0].plot(x, b1['p'][jf, :], '.', ms=3, label='p fluid (mid row)')
    ax[0, 0].plot(x, b2['p'][jp, :], '.', ms=3, label='p porous (mid row)')
    ax[0, 0].plot(x, np.polyval(coef_f, x), '--', color='k', lw=1,
                  label='fit dp/dx=%.4g' % P)
    ax[0, 0].axvspan(0.25 * L, 0.75 * L, alpha=0.08)
    ax[0, 0].set_xlabel('x (m)'); ax[0, 0].set_ylabel('p (Pa)')
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].set_title('(a) pressure drop')

    ax[0, 1].plot(uc_all, yc_all, 'o', ms=3, mfc='none', color='tab:blue',
                  label='CFD (x=Lx/2)')
    ax[0, 1].plot(bj.profile(ya, sol), ya, '-', color='tab:red', lw=1.5,
                  label='analytic (BJ composite)')
    y_ns = np.linspace(0, H_FL, 200)
    ax[0, 1].plot(bj.u_poiseuille_noslip(y_ns, G, MU, H_FL), y_ns, ':',
                  color='gray', lw=1.5, label='no-slip Poiseuille (ref)')
    ax[0, 1].axhline(0, color='k', lw=0.5)
    ax[0, 1].set_xlabel('u (m/s)'); ax[0, 1].set_ylabel('y (m)')
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].set_title('(b) streamwise velocity profile')

    dy2 = 0.6 * sol['lam']
    ax[1, 0].plot(uc_all, yc_all, 'o', ms=3.5, mfc='none', color='tab:blue')
    ax[1, 0].plot(bj.profile(ya, sol), ya, '-', color='tab:red', lw=1.5)
    ax[1, 0].axhline(0, color='k', lw=0.6)
    ax[1, 0].axvline(sol['u_int'], ls=':', color='gray')
    ax[1, 0].axvline(sol['uD'], ls=':', color='green')
    ax[1, 0].set_xlim(0, 1.2 * np.max(uc_all))
    ax[1, 0].set_ylim(-dy2, dy2)
    ax[1, 0].set_xlabel('u (m/s)'); ax[1, 0].set_ylabel('y (m)')
    ax[1, 0].set_title('(c) interface region zoom  (0.6 lambda)')

    for xs0, c in zip(stations, ['tab:orange', 'tab:red', 'tab:green']):
        ii = station(b1, xs0)
        uc = np.concatenate([b2['u'][:, ii], b1['u'][:, ii]])
        ax[1, 1].plot(uc,
                      np.concatenate([b2['yc'][:, ii], b1['yc'][:, ii]]),
                      '.', ms=2.5, color=c, label='x=%.1f' % xs0)
    ax[1, 1].plot(bj.profile(ya, sol), ya, '-', color='k', lw=1)
    ax[1, 1].set_xlabel('u (m/s)'); ax[1, 1].set_ylabel('y (m)')
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].set_title('(d) developed-flow check (3 stations)')

    fig.suptitle('Beavers-Joseph: two-block low-speed/porous validation '
                 '(eps=%.2g, dp=%.2g m, mu=%.3g Pa s)' % (EPS, DP, MU))
    fig.tight_layout()
    fig.savefig(os.path.join(args.case, args.out), dpi=130)
    print('saved', os.path.join(args.case, args.out))


if __name__ == '__main__':
    main()


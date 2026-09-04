import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pv_io import load_vtk, load_ts, cell_centers, brinkman_profile, perme_ergun

def plot_darcy(args, p, u, xs, ys):
    eps, dp, mu = args.eps, args.dp, args.mu
    K = perme_ergun(eps, dp)
    jc = NYC = u.shape[0]
    imid = np.arange(10, 34)
    coef = np.polyfit(xs[imid], p[jc-1, imid], 1)
    dpdx = coef[0]
    i0 = 20
    uprof = u[:, i0]
    ubar = u.mean(axis=0)
    ya = np.linspace(0, args.LY, 201)
    ub, uD, _ = brinkman_profile(ya, dpdx, eps, dp, mu)
    alpha = 0.5*args.LY/(np.sqrt(K)/eps)
    ubarB = uD*(1.0 - np.tanh(alpha)/alpha)
    ubar_mid = ubar[10:34].mean()
    print('eps=%.3g dp=%.3g K=%.3e m2' % (eps, dp, K))
    print('fitted dp/dx = %.6e Pa/m' % dpdx)
    print('u_max CFD = %.6e , ubar CFD = %.6e , ubar Brinkman = %.6e (err %+.2f%%)'
          % (uprof.max(), ubar_mid, ubarB, 100*(ubar_mid-ubarB)/ubarB))
    if args.u_in:
        print('prescribed u_in = %g -> ubar/u_in = %.3f' % (args.u_in, ubar_mid/args.u_in))

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(xs, p[jc-1, :], 'o-', ms=4, label='p centerline')
    ax[0, 0].plot(xs, np.polyval(coef, xs), '--', color='tab:red',
                  label='lin fit dp/dx=%.4g' % dpdx)
    ax[0, 0].set_xlabel('x (m)'); ax[0, 0].set_ylabel('p'); ax[0, 0].legend()
    ax[0, 0].set_title('(a) centerline pressure (pressure-driven plug)')

    ax[0, 1].plot(uprof, ys, 'o-', ms=4, label='CFD u(y) x=L/2')
    ax[0, 1].plot(ub, ya, '--', color='tab:red', label='Brinkman(fit dp/dx)')
    ax[0, 1].axvline(ubar_mid, color='tab:green', ls=':', label='ubar CFD')
    ax[0, 1].axvline(ubarB, color='tab:purple', ls=':', label='ubar Brink')
    ax[0, 1].set_xlabel('u'); ax[0, 1].set_ylabel('y'); ax[0, 1].legend(fontsize=8)
    ax[0, 1].set_title('(b) u(y) vs Darcy-Brinkman solution')

    ax[1, 0].plot(xs, ubar, 'o-', ms=3, label='col-avg u')
    ax[1, 0].axhline(ubarB, color='tab:red', ls='--', label='ubar_Brinkman')
    if args.u_in:
        ax[1, 0].axhline(args.u_in, color='tab:green', ls=':', label='u_in')
    ax[1, 0].set_xlabel('x'); ax[1, 0].set_ylabel('col-avg u'); ax[1, 0].legend()
    ax[1, 0].set_title('(c) streamwise mean u (mass-flow check)')

    ax[1, 1].plot(ys, u[:, i0], 'o-', ms=3)
    ax[1, 1].set_xlabel('y'); ax[1, 1].set_ylabel('u')
    ax[1, 1].set_title('(d) u(y) at x=L/2 (all rows)')
    fig.suptitle('Porous plug validation - Darcy/Brinkman (eps=%.2g, dp=%.2g m)' % (eps, dp))
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=120)
    print('saved', args.out_png)

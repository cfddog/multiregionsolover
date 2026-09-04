import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pv_io import load_vtk, load_ts, cell_centers

def plot_ltne(args, p, T, Ts, u, xs, ys):
    jc = 10; i0 = 20
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(T[:, i0], ys, 'o-', ms=3, label='Tf')
    ax[0, 0].plot(Ts[:, i0], ys, 's--', ms=3, label='Ts')
    ax[0, 0].set_xlabel('T (K)'); ax[0, 0].set_ylabel('y'); ax[0, 0].legend()
    ax[0, 0].set_title('(a) Tf & Ts across channel at x=L/2')

    ax[0, 1].plot(xs, T[jc, :], 'o-', ms=3, label='Tf')
    ax[0, 1].plot(xs, Ts[jc, :], 's--', ms=3, label='Ts')
    ax[0, 1].set_xlabel('x'); ax[0, 1].set_ylabel('T'); ax[0, 1].legend()
    ax[0, 1].set_title('(b) Tf & Ts along centerline')

    dT = T - Ts
    pm = ax[1, 0].pcolormesh(xs, ys, dT, shading='auto')
    fig.colorbar(pm, ax=ax[1, 0])
    ax[1, 0].set_xlabel('x'); ax[1, 0].set_ylabel('y')
    ax[1, 0].set_title('(c) Tf - Ts (K)')

    ax[1, 1].plot(xs, T[0, :], 'o-', ms=3, label='Tf j=0')
    ax[1, 1].plot(xs, Ts[0, :], 's--', ms=3, label='Ts j=0')
    ax[1, 1].plot(xs, T[-1, :], 'o-', ms=3, alpha=.6, label='Tf j=19')
    ax[1, 1].plot(xs, Ts[-1, :], 's--', ms=3, alpha=.6, label='Ts j=19')
    ax[1, 1].set_xlabel('x'); ax[1, 1].set_ylabel('T'); ax[1, 1].legend(fontsize=8)
    ax[1, 1].set_title('(d) near-wall temperatures')

    print('Tf range %.2f..%.2f K, Ts range %.2f..%.2f K'
          % (T.min(), T.max(), Ts.min(), Ts.max()))
    print('center Tf,Ts at inlet/mid/outlet x: ',
          T[jc, 0], Ts[jc, 0], T[jc, 20], Ts[jc, 20], T[jc, 39], Ts[jc, 39])
    fig.suptitle('Porous plug LTNE validation (Tf / Ts, hv>0, wall Tw=400 K)')
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=120)
    print('saved', args.out_png)

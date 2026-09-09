#!/usr/bin/env python3
"""Side-by-side plug/darcy quantitative comparison: SIMPLE reference vs AC.

Reads the two `flow3d_block_1.vtk` outputs and reports, on IDENTICAL
post-processing conventions:
    u_bar : column-averaged streamwise velocity averaged over a fixed window
            (default: cell-centres with x in [0.25*Lx, 0.85*Lx), reproducing
            pv_darcy.py cols 10..33 on the 40-cell grid);
    dp/dx : linear pressure fit over the same window, for the all-row-mean
            pressure and for the last row (pv_darcy.py centreline convention);
    diagnostics: end-cell residuals, interior-line extrapolation to the faces,
            L_eff = dP/|dpdx|.

Usage:
    python3 pv_plug_compare.py --simple <SIMPLE_dir> --ac <AC_dir> \
        [--out png] [--pin 0.25 --pout 0.0] [--xw0 0.25 --xw1 0.85]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_vtk_cell_fields(path):
    """Parse a structured VTK: cell arrays dict + x cell centres.

    Cell data assumed written i-fastest, then j, then k (OpenCFD-EC writer).
    """
    lines = open(path).read().splitlines()
    di = [l for l in lines if l.startswith('DIMENSIONS')][0].split()
    nxp, nyp, nzp = int(di[1]), int(di[2]), int(di[3])
    nxc, nyc, nzc = nxp - 1, nyp - 1, nzp - 1
    nc = nxc * nyc * nzc

    pi = [i for i, l in enumerate(lines) if l.startswith('POINTS')][0]
    node = []
    for s in lines[pi + 1:pi + 1 + nxp * nyp * nzp]:
        t = s.split()
        if len(t) >= 3:
            node.append((float(t[0]), float(t[1]), float(t[2])))
    node = np.asarray(node)
    x0 = node[0:nxp, 0]
    xc = 0.5 * (x0[:-1] + x0[1:])

    def read_scalar(name):
        i = lines.index('SCALARS %s double 1' % name)
        j = lines.index('LOOKUP_TABLE default', i)
        a = np.array([float(s.split()[0]) for s in lines[j + 1:j + 1 + nc]])
        return a.reshape((nzc, nyc, nxc))

    def read_vectors(name):
        i = lines.index('VECTORS %s double' % name)
        v = np.zeros((nc, 3))
        for n, s in enumerate(lines[i + 1:i + 1 + nc]):
            t = s.split()
            if len(t) >= 3:
                v[n] = float(t[0]), float(t[1]), float(t[2])
        return v.reshape((nzc, nyc, nxc, 3))

    fields = {'p': read_scalar('pressure'), 'T': read_scalar('temperature'),
              'vel': read_vectors('velocity')}
    return fields, xc


def analyze(fields, xc, xw0, xw1, pin, pout):
    p = fields['p']
    p = p[0, :, :] if p.ndim == 3 else p
    u = fields['vel'][..., 0]
    u = u[0, :, :] if u.ndim == 3 else u
    nyc, nxc = p.shape
    colmean = u.mean(axis=0)
    Lx = float(xc[-1] + 0.5 * (xc[1] - xc[0]))
    cols = np.where((xc >= xw0 * Lx) & (xc < xw1 * Lx))[0]
    if cols.size == 0:
        raise RuntimeError('empty fit window')
    pmean = p.mean(axis=0)
    co_all = np.polyfit(xc[cols], pmean[cols], 1)
    co_row = np.polyfit(xc[cols], p[-1, cols], 1)
    fit_all = np.polyval(co_all, xc)
    r = {
        'nx': int(nxc), 'ny': int(nyc),
        'xw_first': int(cols[0]), 'xw_last': int(cols[-1]),
        'ubar': float(colmean[cols].mean()),
        'dpdx_all': float(co_all[0]), 'dpdx_row': float(co_row[0]),
        'res0': float(pmean[0] - fit_all[0]), 'resN': float(pmean[-1] - fit_all[-1]),
        'pf_in': float(np.polyval(co_all, 0.0)), 'pf_out': float(np.polyval(co_all, Lx)),
        'Leff': float((pin - pout) / max(abs(co_all[0]), 1e-30)),
    }
    return r, pmean, colmean, xc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--simple', required=True)
    ap.add_argument('--ac', required=True)
    ap.add_argument('--out', default=None)
    ap.add_argument('--pin', type=float, default=0.25)
    ap.add_argument('--pout', type=float, default=0.0)
    ap.add_argument('--xw0', type=float, default=0.25)
    ap.add_argument('--xw1', type=float, default=0.85)
    ap.add_argument('--label', default='')
    args = ap.parse_args()

    rows = []
    for tag, d in [('SIMPLE', args.simple), ('AC', args.ac)]:
        vtk = os.path.join(d, 'flow3d_block_1.vtk')
        f, xc = read_vtk_cell_fields(vtk)
        r, pmean, colmean, xc = analyze(f, xc, args.xw0, args.xw1,
                                        args.pin, args.pout)
        rows.append((tag, r, pmean, colmean, xc))
        print('[%s] nx=%d ny=%d window cols %d..%d'
              % (tag, r['nx'], r['ny'], r['xw_first'], r['xw_last']))
        print('    u_bar            = % .6e' % r['ubar'])
        print('    dp/dx (all-rows) = % .6e' % r['dpdx_all'])
        print('    dp/dx (last row) = % .6e' % r['dpdx_row'])
        print('    p(x=0),p(x=L)    = % .6f , % .6f  (interior line)' % (r['pf_in'], r['pf_out']))
        print('    end-cell resid   = % .6f , % .6f' % (r['res0'], r['resN']))
        print('    L_eff            = %.4f m' % r['Leff'])

    _, rs, _, _, _ = rows[0]
    _, ra, _, _, _ = rows[1]
    print('\n== AC vs SIMPLE (relative, %% of SIMPLE) ==')
    for key, name in [('ubar', 'u_bar'), ('dpdx_all', 'dp/dx(all-rows)'),
                      ('dpdx_row', 'dp/dx(last-row)')]:
        rel = 100.0 * (ra[key] - rs[key]) / abs(rs[key])
        print('%-19s SIM % .6e   AC % .6e   rel %+.3f %%' % (name, rs[key], ra[key], rel))
    ok_u = abs(100.0 * (ra['ubar'] - rs['ubar']) / abs(rs['ubar'])) < 1.0
    ok_d = abs(100.0 * (ra['dpdx_all'] - rs['dpdx_all']) / abs(rs['dpdx_all'])) < 1.0
    print('AC within 1%% of SIMPLE?  u_bar: %s   dp/dx: %s' % ('YES' if ok_u else 'NO',
                                                               'YES' if ok_d else 'NO'))

    if args.out:
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        for tag, r, pmean, colmean, xc in rows:
            ax[0, 0].plot(xc, pmean, 'o-', ms=3, label=tag)
            ax[1, 0].plot(xc, colmean, '.-', ms=4, label=tag)
        ax[0, 0].legend(); ax[0, 0].set_title('(a) p along x (mean over rows)')
        ax[0, 0].set_xlabel('x'); ax[0, 0].set_ylabel('p')
        ax[1, 0].legend(); ax[1, 0].set_title('(b) column-averaged u')
        ax[1, 0].set_xlabel('x'); ax[1, 0].set_ylabel('u')
        _, _, _, cs, _ = rows[0]
        _, _, _, ca, _ = rows[1]
        ax[1, 1].plot(xc, 100.0 * (ca / cs - 1.0), 'o-', ms=4)
        ax[1, 1].axhline(0, color='k', lw=0.6)
        ax[1, 1].axhspan(-1, 1, color='g', alpha=0.15)
        ax[1, 1].set_title('(c) AC/SIMPLE col-mean u deviation (%)')
        ax[1, 1].set_xlabel('x')
        ax[0, 1].axis('off')
        fig.suptitle('plug/darcy SIMPLE vs AC ' + args.label)
        fig.tight_layout()
        fig.savefig(args.out, dpi=120)
        print('saved', args.out)


if __name__ == '__main__':
    main()


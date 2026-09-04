import argparse, os
from pv_io import load_vtk, load_ts, cell_centers
from pv_darcy import plot_darcy
from pv_ltne import plot_ltne

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case_dir')
    ap.add_argument('out_png')
    ap.add_argument('--mode', default='darcy', choices=['darcy', 'ltne'])
    ap.add_argument('--eps', type=float, default=0.9)
    ap.add_argument('--dp', type=float, default=0.10)
    ap.add_argument('--mu', type=float, default=0.02)
    ap.add_argument('--u_in', type=float, default=None)
    ap.add_argument('--LY', type=float, default=1.0)
    args = ap.parse_args()

    p, T, uvw = load_vtk(os.path.join(args.case_dir, 'flow3d_block_1.vtk'))
    Ts = load_ts(os.path.join(args.case_dir, 'Ts_block_1.dat'))
    u = uvw[:, :, 0]
    xs, ys = cell_centers()

    if args.mode == 'darcy':
        plot_darcy(args, p, u, xs, ys)
    else:
        plot_ltne(args, p, T, Ts, u, xs, ys)

if __name__ == '__main__':
    main()

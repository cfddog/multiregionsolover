#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the two-block Beavers-Joseph validation case for OpenCFD-EC.

Domain (physical, SI, Lscale=1 m):
   porous bed   : 0 <= x <= Lx , -Hp <= y <= 0   (block 2, BLOCK_POROUS)
   free channel : 0 <= x <= Lx ,   0 <= y <= h   (block 1, BLOCK_LOWSPEED)
   quasi-2D     : z in [0, dz] (2 nodes, symmetric k- and k+ faces)

Block 1 (fluid) j- face (node plane j=1, y=0) touches block 2 (porous)
j+ face (node plane j=nyp, y=0).  Node ranges in x and z are identical on
both blocks (conformal, index-aligned interface) and the cell size normal
to the interface is identical on the two sides (dyf = dyp).

Outputs written in this directory:
   Mesh3d.x            Fortran sequential unformatted, 2 blocks
   bc3d.inp            boundary file (physical + interface code 16)
   bc3d_interface.inp  same records (interface pairing for bc_msg2)
   material.in         block types: 2=fluid(low-speed), 3=porous
   porous.inp          porous block eps / dp / hv
   control.ec          solver control namelist
"""
import argparse
import struct
import math

# ----------------------------------------------------------------------
# physical case constants (m, SI)
LX   = 8.0       # streamwise length
H_FL = 1.0       # free-fluid channel height (y in 0..h)
HP   = 1.2       # porous-bed thickness       (y in -Hp..0)
DZ   = 1.0       # z thickness (one cell, symmetric z faces)
EPS  = 0.8       # porosity of the bed
DP   = 0.5       # particle diameter (m) -> Ergun permeability
MU   = 0.02      # dynamic viscosity (Pa s), LS_mu
P_IN = 0.005     # pressure-inlet value (Pa), LS_P_in
P_OUT = 0.0      # outlet pressure (Pa), LS_P_out


def ergun_permeability(eps, dp):
    return dp**2 * eps**3 / (150.0 * (1.0 - eps)**2)


def write_fortran_record(f, data, fmt):
    if isinstance(data, bytes):
        raw = data
    elif fmt == '':
        raw = data.tobytes()
    else:
        raw = struct.pack(fmt, *data)
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def gen_mesh3d(path, nx, nyf, nyp, nz, lx=LX):
    """Write Mesh3d.x: 2 blocks (fluid block 1 first, then porous block 2)."""
    dx = lx / (nx - 1)
    dyf = H_FL / (nyf - 1)
    dyp = HP / (nyp - 1)
    xs = [i * dx for i in range(nx)]
    zs = [0.0, DZ]
    with open(path, 'wb') as f:
        write_fortran_record(f, [2], 'i')
        write_fortran_record(f, [nx, nyf, nz, nx, nyp, nz], '6i')
        # block 1 (fluid): y in [0, H_FL], y(1)=0 (interface), y(nyf)=H_FL
        xc, yc, zc = [], [], []
        for k in range(nz):
            for j in range(nyf):
                yj = j * dyf
                for i in range(nx):
                    xc.append(xs[i]); yc.append(yj); zc.append(zs[k])
        write_fortran_record(f, struct.pack('<%dd' % len(xc), *xc), '')
        write_fortran_record(f, struct.pack('<%dd' % len(yc), *yc), '')
        write_fortran_record(f, struct.pack('<%dd' % len(zc), *zc), '')
        # block 2 (porous): y in [-Hp, 0], y(1)=-Hp, y(nyp)=0 (interface)
        xc, yc, zc = [], [], []
        for k in range(nz):
            for j in range(nyp):
                yj = -HP + j * dyp
                for i in range(nx):
                    xc.append(xs[i]); yc.append(yj); zc.append(zs[k])
        write_fortran_record(f, struct.pack('<%dd' % len(xc), *xc), '')
        write_fortran_record(f, struct.pack('<%dd' % len(yc), *yc), '')
        write_fortran_record(f, struct.pack('<%dd' % len(zc), *zc), '')
    return dx, dyf, dyp


def gen_bc(path, nx, nyf, nyp, nz, bc_if):
    """Write bc3d.inp / bc3d_interface.inp (identical records)."""
    with open(path, 'w') as f:
        f.write('! OpenCFD-EC boundary file: Beavers-Joseph two blocks\n')
        f.write('%d\n' % 2)
        # ------------- block 1 : free-fluid channel (low-speed) ---------
        f.write('%d %d %d\n' % (nx, nyf, nz))
        f.write('blk-1\n')
        f.write('%d\n' % 6)
        faces = [
            (1, 1, 1, nyf, 1, nz, 5),          # i-  pressure inlet
            (nx, nx, 1, nyf, 1, nz, 6),         # i+  outlet
            (1, nx, 1, 1, 1, nz, bc_if),        # j-  interface (porous)
            (1, nx, nyf, nyf, 1, nz, 2),        # j+  top wall
            (1, nx, 1, nyf, 1, 1, 3),           # k-  symmetry
            (1, nx, 1, nyf, nz, nz, 3),         # k+  symmetry
        ]
        for ib, ie, jb, je, kb, ke, bc in faces:
            f.write('%4d %4d %4d %4d %4d %4d %4d\n' % (ib, ie, jb, je, kb, ke, bc))
            if bc == bc_if:
                # neighbour = block 2, its j+ face: i 1..nx, j=nyp..nyp, k 1..nz
                f.write('%4d %4d %4d %4d %4d %4d %4d\n' % (1, nx, nyp, nyp, 1, nz, 2))
        # ------------- block 2 : porous bed ------------------------------
        f.write('%d %d %d\n' % (nx, nyp, nz))
        f.write('blk-2\n')
        f.write('%d\n' % 6)
        faces = [
            (1, 1, 1, nyp, 1, nz, 5),          # i-  pressure inlet
            (nx, nx, 1, nyp, 1, nz, 6),         # i+  outlet
            (1, nx, 1, 1, 1, nz, 2),            # j-  bottom wall
            (1, nx, nyp, nyp, 1, nz, bc_if),    # j+  interface (fluid)
            (1, nx, 1, nyp, 1, 1, 3),           # k-  symmetry
            (1, nx, 1, nyp, nz, nz, 3),         # k+  symmetry
        ]
        for ib, ie, jb, je, kb, ke, bc in faces:
            f.write('%4d %4d %4d %4d %4d %4d %4d\n' % (ib, ie, jb, je, kb, ke, bc))
            if bc == bc_if:
                # neighbour = block 1, its j- face: i 1..nx, j=1..1, k 1..nz
                f.write('%4d %4d %4d %4d %4d %4d %4d\n' % (1, nx, 1, 1, 1, nz, 1))


def gen_material(path):
    with open(path, 'w') as f:
        f.write('2\n')
        f.write('2 3\n')                 # block1 = LOWSPEED, block2 = POROUS
        f.write('0.0 0.0 0.0\n')         # block1 solid props (unused)
        f.write('2500.0 800.0 300.0\n')  # block2 skeleton rho,Cp,k


def gen_porous(path):
    with open(path, 'w') as f:
        f.write('1\n')
        f.write('2  %.4f  %.4f  0.0\n' % (EPS, DP))   # hv=0 -> frozen frame


def gen_control(path, t_end, max_iter, tol, ksave):
    with open(path, 'w') as f:
        f.write('$control_ec\n')
        f.write('  Iflag_init=0\n')
        f.write('  Ma=0.1d0\n')
        f.write('  Re=1000.d0\n')
        f.write('  AoA=0.d0\n')
        f.write('  If_viscous=1\n')
        f.write('  Iflag_turbulence_model=0\n')
        f.write('  Kstep_save=%d\n' % ksave)
        f.write('  t_end=%.0f.d0\n' % t_end)
        f.write('  CFL=1.0d0\n')
        f.write('  Time_Method=0\n')
        f.write('  T_inf=288.15d0\n')
        f.write('  Twall=-1.d0\n')
        f.write('  Iflag_Scheme=5\n')
        f.write('  Iflag_Flux=5\n')
        f.write('  IFlag_Reconstruction=0\n')
        f.write('  Ref_S=1.d0\n')
        f.write('  Ref_L=1.0d0\n')
        f.write('  Lscale=1.0d0\n')
        f.write('  Cood_Y_UP=1\n')
        f.write('  NUM_THREADS=1\n')
        f.write('  Mesh_File_Format=2\n')
        f.write('  Kstep_show=1\n')
        f.write('  Kstep_average=0\n')
        f.write('  Kstep_smooth=-1\n')
        f.write('  Kstep_init_smooth=0\n')
        f.write('  Iflag_local_dt=1\n')
        f.write('  dt_global=1.0d0\n')
        f.write('  dtmax=1000.d0\n')
        f.write('  dtmin=1.d-9\n')
        f.write('  If_Residual_smoothing=0\n')
        f.write('  w_LU=1.0d0\n')
        f.write('  If_dtime_mesh=1\n')
        f.write('  Num_Mesh=1\n')
        f.write('  IF_Debug=0\n')
        f.write('  Step_Inner_Limit=20\n')
        f.write('  Res_Inner_Limit=1.d-10\n')
        f.write('  MUT_MAX=-1.d0\n')
        f.write('  Bound_Scheme=-1\n')
        f.write('  Pre_Step_Mesh=0,0,0\n')
        f.write('  Centroid=0.d0,0.d0,0.d0\n')
        f.write('  IFLAG_LIMIT_FLOW=0\n')
        f.write('  Pdebug=1,1,1,1\n')
        f.write('  Ldmin=1.d-6\n')
        f.write('  Ldmax=1000.d0\n')
        f.write('  Lpmin=1.d-6\n')
        f.write('  Lpmax=1000.d0\n')
        f.write('  Lumax=1000.d0\n')
        f.write('  LSAmax=1000.d0\n')
        f.write('  CP1_NSA=0.2d0\n')
        f.write('  CP2_NSA=100.d0\n')
        f.write('  Periodic_dX=0.d0\n')
        f.write('  Periodic_dY=0.d0\n')
        f.write('  Periodic_dZ=0.d0\n')
        f.write('  IF_TurboMachinary=0\n')
        f.write('  Ref_medium_usrdef=0\n')
        f.write('  IF_Scheme_Positivity=1\n')
        f.write('  Turbo_P0=101330.d0\n')
        f.write('  Turbo_T0=288.15d0\n')
        f.write('  Turbo_L0=1.d0\n')
        f.write('  Turbo_w=0.d0\n')
        f.write('  Turbo_Periodic_seta=0.d0\n')
        f.write('  gamma=1.4d0\n')
        f.write('  PrL=0.7d0\n')
        f.write('  PrT=0.9d0\n')
        f.write('  IF_Innerflow=0\n')
        f.write('  Iflag_savefile=0\n')
        f.write('  LS_rho=1.0d0\n')
        f.write('  LS_mu=%.6fd0\n' % MU)
        f.write('  LS_k=0.0d0\n')
        f.write('  LS_Cp=1.0d0\n')
        f.write('  LS_T_ref=300.d0\n')
        f.write('  LS_Inlet_Type=3\n')
        f.write('  LS_U_in=0.0d0\n')
        f.write('  LS_V_in=0.d0\n')
        f.write('  LS_W_in=0.d0\n')
        f.write('  LS_Mdot_in=0.d0\n')
        f.write('  LS_P_in=%.6fd0\n' % P_IN)
        f.write('  LS_P_out=%.6fd0\n' % P_OUT)
        f.write('  LS_T_wall=-1.d0\n')
        f.write('  LS_U_lid=0.0d0\n')
        f.write('  LS_alpha_p=0.2d0\n')
        f.write('  LS_alpha_u=0.5d0\n')
        f.write('  LS_alpha_T=0.7d0\n')
        f.write('  LS_Max_Iter=%d\n' % max_iter)
        _tol = '1.d%d' % int(round(math.log10(tol)))
        f.write('  LS_Tol=%s\n' % _tol)
        f.write('  LS_Scheme=1\n')
        f.write('  LS_Algorithm=1\n')
        f.write('  Porous_alpha_Ts=0.7d0\n')
        f.write('  Porous_Max_Iter=%d\n' % max_iter)
        f.write('  Porous_Tol=%s\n' % _tol)
        f.write('  Porous_T_ref=300.d0\n')
        f.write('$end\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nx', type=int, default=81, help='x node count')
    ap.add_argument('--nyf', type=int, default=41, help='fluid block y node count')
    ap.add_argument('--nyp', type=int, default=49, help='porous block y node count')
    ap.add_argument('--lx', type=float, default=LX, help='streamwise length')
    ap.add_argument('--tend', type=float, default=150.0, help='outer steps')
    ap.add_argument('--maxiter', type=int, default=80, help='SIMPLE inner sweeps per block call')
    ap.add_argument('--tol', type=float, default=1.e-9)
    ap.add_argument('--ksave', type=int, default=100)
    args = ap.parse_args()

    dx, dyf, dyp = gen_mesh3d('Mesh3d.x', args.nx, args.nyf, args.nyp, 2, args.lx)
    gen_bc('bc3d.inp', args.nx, args.nyf, args.nyp, 2, 16)
    gen_bc('bc3d_interface.inp', args.nx, args.nyf, args.nyp, 2, 16)
    gen_material('material.in')
    gen_porous('porous.inp')
    gen_control('control.ec', args.tend, args.maxiter, args.tol, args.ksave)

    K = ergun_permeability(EPS, DP)
    Keff = K / EPS**2
    lam = math.sqrt(Keff)
    kappa = 1.0 / lam
    G = P_IN / args.lx
    uD = G * Keff / MU
    print('Mesh3d.x    : block1(fluid) %dx%dx%d , block2(porous) %dx%dx%d' %
          (args.nx, args.nyf, 2, args.nx, args.nyp, 2))
    print('cell sizes  : dx=%.5f  dy_fluid=%.5f  dy_porous=%.5f (match=%s)'
          % (dx, dyf, dyp, abs(dyf - dyp) < 1.e-12))
    print('eps=%.3f dp=%.3f m -> Ergun K=%.6e m2, K_eff=K/eps^2=%.6e m2'
          % (EPS, DP, K, Keff))
    print('lambda=sqrt(K_eff)=%.5f m , kappa=1/lambda=%.4f 1/m' % (lam, kappa))
    print('bed depth Hp=%.3f m -> kappa*Hp=%.2f (>=6 -> deep-bed BJ limit OK)' % (HP, kappa * HP))
    print('fluid height h=%.3f m -> kappa*h=%.2f' % (H_FL, kappa * H_FL))
    print('nominal dp/dx=%.4e Pa/m, u_Darcy=%.4e m/s' % (-G, uD))
    print('cells per lambda near interface: %.2f (fluid), %.2f (porous)'
          % (lam / dyf, lam / dyp))


if __name__ == '__main__':
    main()


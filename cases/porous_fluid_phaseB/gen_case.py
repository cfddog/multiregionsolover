#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-block 'transpiration wall under a hot supersonic channel' mesh (Phase B,
interface code 19, BLOCK_FLUID <-> BLOCK_POROUS blowing coupling).

Domain (SI, Lscale=1 m):
   compressible channel : 0<=x<=LX ,  0<=y<=HC    (block 1, BLOCK_FLUID  =0)
   porous wall          : 0<=x<=LX , -HP<=y<=0    (block 2, BLOCK_POROUS =3)
   quasi-2D             : z in [0,DZ] (2 nodes, symmetric z faces)

Interface at y=0: compressible block j- face  <->  porous block j+ face
(conformal, aligned, same dy).  Interface code = 19 written in bc3d.inp and
bc3d_interface.inp (pairing records like cases/beavers_joseph).

Flow BCs:
  compressible: x- Ma=3/800K inlet(5), x+ outlet(6), y+ symmetry(3), z 3;
                y- = interface(19)
  porous:       y- bottom = coolant velocity inlet (5, LS_Inlet_Type=1,
                normal +y = LS_V_in = G/LS_rho, 300 K); x walls(2);
                z symmetry(3); y+ = interface(19)
"""
import argparse
import struct

LX  = 0.5
HC  = 0.05
HP  = 0.010
DZ  = 1.0
EPS = 0.3
DP  = 1.0e-3
HV  = 2.0e6
G   = 2.0                 # coolant mass flux kg/m2/s
LS_RHO = 0.3172           # air at 300 K, p~27.3 kPa
KS  = 13.4
BC_IF = 19


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


def gen_mesh3d(path, nx, nyc, nyp, nz, lx=LX):
    dx = lx / (nx - 1)
    dyc = HC / (nyc - 1)
    dyp = HP / (nyp - 1)
    xs = [i * dx for i in range(nx)]
    zs = [0.0, DZ]
    assert abs(dyc - dyp) < 1.e-12, (dyc, dyp)
    with open(path, 'wb') as f:
        write_fortran_record(f, [2], 'i')
        write_fortran_record(f, [nx, nyc, nz, nx, nyp, nz], '6i')
        # block 1 compressible: y in [0,HC]; y(1)=0 interface
        for lo, hi, ny in ((0.0, HC, nyc), (-HP, 0.0, nyp)):
            xc, yc, zc = [], [], []
            for k in range(nz):
                for j in range(ny):
                    yj = lo + j * (hi - lo) / (ny - 1)
                    for i in range(nx):
                        xc.append(xs[i]); yc.append(yj); zc.append(zs[k])
            write_fortran_record(f, struct.pack('<%dd' % len(xc), *xc), '')
            write_fortran_record(f, struct.pack('<%dd' % len(yc), *yc), '')
            write_fortran_record(f, struct.pack('<%dd' % len(zc), *zc), '')
    return dx, dyc, dyp


def gen_bc(path, nx, nyc, nyp, nz):
    def face(ib, ie, jb, je, kb, ke, bc):
        return '%4d %4d %4d %4d %4d %4d %4d\n' % (ib, ie, jb, je, kb, ke, bc)

    with open(path, 'w') as f:
        f.write('! OpenCFD-EC boundary: transpiration wall under Ma=3 800K channel (interface 19)\n')
        f.write('2\n')
        # ---------------- block 1 : compressible channel -------------------
        f.write('%d %d %d\n' % (nx, nyc, nz))
        f.write('blk-1\n')
        f.write('6\n')
        f.write(face(1, 1, 1, nyc, 1, nz, 5))          # i-  supersonic inlet
        f.write(face(nx, nx, 1, nyc, 1, nz, 6))        # i+  outlet
        f.write(face(1, nx, 1, 1, 1, nz, BC_IF))       # j-  interface (porous)
        f.write(face(1, nx, nyp, nyp, 1, nz, 2))      # partner: block2 j+ face
        f.write(face(1, nx, nyc, nyc, 1, nz, 3))       # j+  top symmetry
        f.write(face(1, nx, 1, nyc, 1, 1, 3))          # k-  symmetry
        f.write(face(1, nx, 1, nyc, nz, nz, 3))        # k+  symmetry
        # ---------------- block 2 : porous wall ----------------------------
        f.write('%d %d %d\n' % (nx, nyp, nz))
        f.write('blk-2\n')
        f.write('6\n')
        f.write(face(1, 1, 1, nyp, 1, nz, 2))          # i-  wall
        f.write(face(nx, nx, 1, nyp, 1, nz, 2))        # i+  wall
        f.write(face(1, nx, 1, 1, 1, nz, 5))           # j-  coolant inlet
        f.write(face(1, nx, nyp, nyp, 1, nz, BC_IF))   # j+  interface (fluid)
        f.write(face(1, nx, 1, 1, 1, nz, 1))           # partner: block1 j- face
        f.write(face(1, nx, 1, nyp, 1, 1, 3))          # k-  symmetry
        f.write(face(1, nx, 1, nyp, nz, nz, 3))        # k+  symmetry


def gen_material(path):
    with open(path, 'w') as f:
        f.write('2\n')
        f.write('0 3\n')                      # block1 FLUID, block2 POROUS
        f.write('0.0 0.0 0.0\n')              # block1 solid props (unused)
        f.write('8000.0 500.0 %.4f\n' % KS)   # block2 skeleton rho Cp k


def gen_porous(path):
    with open(path, 'w') as f:
        f.write('1\n')
        f.write('2  %.4f  %.4e  %.3e\n' % (EPS, DP, HV))


def gen_control(path, tend, maxiter, ksave, gflag):
    v = G / LS_RHO          # coolant normal velocity in the porous block
    with open(path, 'w') as f:
        f.write('$control_ec\n')
        f.write('  Iflag_init=0\n')
        f.write('  Ma=3.0d0\n')
        f.write('  Re=5.584d6\n')
        f.write('  AoA=0.d0\n')
        f.write('  If_viscous=0\n')           # Phase B first pass (inviscid hot channel)
        f.write('  Iflag_turbulence_model=0\n')
        f.write('  Kstep_save=%d\n' % ksave)
        f.write('  t_end=%.0f.d0\n' % tend)
        f.write('  CFL=0.8d0\n')
        f.write('  Time_Method=0\n')
        f.write('  T_inf=800.d0\n')
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
        f.write('  Iflag_vtk_onefile=0\n')
        f.write('  Iflag_vtk_SI=1\n')
        f.write('  LS_rho=%.5fd0\n' % LS_RHO)
        f.write('  LS_mu=3.17d-3\n')
        f.write('  LS_k=2.62d-2\n')
        f.write('  LS_Cp=1005.d0\n')
        f.write('  LS_T_ref=300.d0\n')
        f.write('  LS_Inlet_Type=1\n')
        f.write('  LS_U_in=0.d0\n')
        f.write('  LS_V_in=%.4fd0\n' % (v if gflag else 0.0))
        f.write('  LS_W_in=0.d0\n')
        f.write('  LS_Mdot_in=0.d0\n')
        f.write('  LS_P_in=2.7314d4\n')
        f.write('  LS_P_out=2.7314d4\n')
        f.write('  LS_T_wall=-1.d0\n')
        f.write('  LS_U_lid=0.d0\n')
        f.write('  LS_alpha_p=0.2d0\n')
        f.write('  LS_alpha_u=0.5d0\n')
        f.write('  LS_alpha_T=1.0d0\n')
        f.write('  LS_Max_Iter=%d\n' % maxiter)
        f.write('  LS_Tol=1.d-7\n')
        f.write('  LS_Scheme=1\n')
        f.write('  LS_Algorithm=1\n')
        f.write('  Porous_alpha_Ts=1.0d0\n')
        f.write('  Porous_Max_Iter=%d\n' % maxiter)
        f.write('  Porous_Tol=1.d-7\n')
        f.write('  Porous_T_ref=300.d0\n')
        f.write('$end\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nx', type=int, default=121)
    ap.add_argument('--nyc', type=int, default=81)
    ap.add_argument('--nyp', type=int, default=17)
    ap.add_argument('--tend', type=float, default=60.0)
    ap.add_argument('--maxiter', type=int, default=1500)
    ap.add_argument('--ksave', type=int, default=15)
    ap.add_argument('--blow', type=int, default=1, help='1 = coolant blowing (LS_V_in>0)')
    args = ap.parse_args()
    dx, dyc, dyp = gen_mesh3d('Mesh3d.x', args.nx, args.nyc, args.nyp, 2)
    gen_bc('bc3d.inp', args.nx, args.nyc, args.nyp, 2)
    gen_bc('bc3d_interface.inp', args.nx, args.nyc, args.nyp, 2)
    gen_material('material.in')
    gen_porous('porous.inp')
    gen_control('control.ec', args.tend, args.maxiter, args.ksave, args.blow)
    print('block1(compressible) %dx%dx%d ; block2(porous) %dx%dx%d' %
          (args.nx, args.nyc, 2, args.nx, args.nyp, 2))
    print('dy=%.3e m, dx=%.4e m ; coolant V_in=%.3f m/s (G=%.2f kg/m2/s, rho=%.4f)'
          % (dyp, dx, G / LS_RHO if args.blow else 0.0, G, LS_RHO))


if __name__ == '__main__':
    main()

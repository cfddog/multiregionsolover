!===============================================================================
! Artificial-compressibility (AC) solver for BLOCK_POROUS blocks.
! Selected when  LS_Algorithm = 3  (SIMPLE/SIMPLEC porous path unchanged).
!
! Same cell-centred FV / AC-LU-SGS engine as the low-speed AC solver
! (src/sub_lowspeed_ac.f90: module lows_ac_work + ac_* routines), plus:
!   - Darcy + Forchheimer drag in the momentum equation (Ergun K, F; the same
!     volume-averaged coefficients as sub_porous.f90:841-853)
!       du/dtau + div(uu+qI) = nu Lap(u) - c_drag u
!       c_drag = (eps^2*mu/K + eps^3*F*rho*|u|/sqrt(K)) / rho
!     -> residual: RAC(m) -= c_drag*Vol*u_m ; LU-SGS diagonal: + c_drag*Vol
!   - LTNE (Tf in B%U(5), skeleton Ts) scalar pair solved after the flow
!     converges (placeholder; enabled when hv > 0).
! Variables/pressure/output identical to the SIMPLE porous path (B%p, U(2..5)).
!===============================================================================

  subroutine porous_ac_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_var
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock
   real(PRE_EC) :: Sfac, Sfac1
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, iter, i,j,k, isdual
   real(PRE_EC) :: uin_x, uin_y, uin_z, Uref, beta, U2, L
   real(PRE_EC) :: rq, rm
   logical :: have_pbc

   B => Mesh(nMesh)%Block(mBlock)
   nx = B%nx; ny = B%ny; nz = B%nz
   L  = max(Lscale, 1.d-30)

   if(.not. allocated(XW) .or. ac_nxw/=nx .or. ac_nyw/=ny .or. ac_nzw/=nz &
      .or. ac_lap /= LAP) then
     if(allocated(XW)) deallocate(XW,RAC,DU4,DTAC,SIG,dragc)
     allocate( XW(4, 1-LAP:nx+LAP-1, 1-LAP:ny+LAP-1, 1-LAP:nz+LAP-1) )
     allocate( RAC(4, nx, ny, nz), DU4(4, nx, ny, nz) )
     allocate( DTAC(nx, ny, nz), SIG(nx, ny, nz), dragc(nx, ny, nz) )
     ac_nxw=nx; ac_nyw=ny; ac_nzw=nz; ac_lap=LAP
   endif

   call lowspeed_inlet_velocity(B, uin_x, uin_y, uin_z)
   Uref = sqrt(uin_x*uin_x + uin_y*uin_y + uin_z*uin_z)
   if(Uref < 1.d-12) Uref = max(abs(LS_U_lid), 1.d-3)
   if(Uref < 1.d-12) Uref = sqrt(abs(LS_P_in-LS_P_out)/max(LS_rho,1.d-20))
   if(Uref < 1.d-12) Uref = 1.d0
   beta = max(AC_beta, 1.d-8) * Uref * Uref
   U2   = max(Uref*Uref, 1.d-20)

   isdual = 0
   if(Time_Method .eq. Time_Dual_LU_SGS) isdual = 1
   call ac_have_pressure_bc(nMesh, mBlock, have_pbc)

   if(my_id .eq. 0) print*, ' AC porous solver (LS_Algorithm=3) block', mBlock, &
      ' eps=', B%porous_eps, ' dp=', B%porous_dp, ' hv=', B%porous_hv, &
      ' beta=', beta, ' CFL=', AC_CFL

   do iter = 1, AC_Max_Iter
     call ac_fill_ghost(nMesh, mBlock, uin_x, uin_y, uin_z)
     call ac_load_state(nMesh, mBlock)
     call ac_boundary_flux(nMesh, mBlock, beta, uin_x, uin_y, uin_z)
     call ac_internal_flux(nMesh, mBlock, beta)
     if(If_viscous .eq. 1 .and. LS_mu .gt. 0.d0) call ac_viscous_res(nMesh, mBlock)
     call ac_porous_drag_res(nMesh, mBlock)              ! sets dragc + RAC(2:4)
     if(isdual .eq. 1) call ac_dual_res(nMesh, mBlock, Sfac)
     call ac_compute_dt(nMesh, mBlock, beta)
     call ac_res_norm(nMesh, mBlock, beta, U2, L, rq, rm)
     if(iter .eq. 1) then
       if(my_id .eq. 0) print*, '  AC-porous iter', iter, ' res_q=', rq, ' res_m=', rm
     endif
     call ac_lusgs_sweep(nMesh, mBlock, beta, isdual, Sfac1)
     do k = 1, nz-1
     do j = 1, ny-1
     do i = 1, nx-1
       B%p(i,j,k)    = LS_rho*XW(1,i,j,k)
       B%U(2,i,j,k)  = XW(2,i,j,k)
       B%U(3,i,j,k)  = XW(3,i,j,k)
       B%U(4,i,j,k)  = XW(4,i,j,k)
       B%U(1,i,j,k)  = LS_rho
     enddo; enddo; enddo
     if(.not. have_pbc) call ac_zero_mean_p(nMesh, mBlock)
     if(iter .eq. 1 .or. mod(iter, AC_Print) .eq. 0) then
       if(my_id .eq. 0) print*, '  AC-porous iter', iter, ' res_q=', rq, ' res_m=', rm
     endif
     if(iter .ge. 20 .and. rq .lt. AC_Tol .and. rm .lt. AC_Tol) exit
   enddo
   if(my_id .eq. 0) print*, ' AC porous solver block', mBlock, ': iterations =', iter, &
       ' final res_q=', rq, ' res_m=', rm

!  LTNE two-temperature scalar pair (Tf / Ts) after the flow has converged.
   if(B%porous_hv .gt. 0.d0 .or. LS_k .gt. 0.d0) then
     call ac_porous_ltne(nMesh, mBlock)
   endif
  end subroutine porous_ac_solver_one_block

!===============================================================================
! Darcy + Forchheimer drag (volume-averaged, kinematic form used by the AC
! momentum equation).  Sets the implicit coefficient dragc and adds
! -c_drag*Vol*u_i to RAC(2:4).
!===============================================================================
  subroutine ac_porous_drag_res(nMesh, mBlock)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k
   real(PRE_EC) :: epsp, Kp, Fc, vm, cdrag, mu, rho
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   rho = max(LS_rho, 1.d-20); mu = LS_mu
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     epsp = B%porous_eps
     Kp = B%porous_dp**2 * epsp**3 / max(150.d0*(1.d0-epsp)**2, 1.d-30)
     Fc = 1.75d0 / max(sqrt(150.d0)*epsp**1.5d0, 1.d-30)
     vm = sqrt(B%U(2,i,j,k)**2 + B%U(3,i,j,k)**2 + B%U(4,i,j,k)**2)
     cdrag = (epsp*epsp*mu/max(Kp,1.d-30) &
            + epsp**3 * Fc * rho * vm / max(sqrt(Kp),1.d-30)) / rho
     dragc(i,j,k) = cdrag
     RAC(2,i,j,k) = RAC(2,i,j,k) - cdrag*B%Vol(i,j,k)*B%U(2,i,j,k)
     RAC(3,i,j,k) = RAC(3,i,j,k) - cdrag*B%Vol(i,j,k)*B%U(3,i,j,k)
     RAC(4,i,j,k) = RAC(4,i,j,k) - cdrag*B%Vol(i,j,k)*B%U(4,i,j,k)
   enddo; enddo; enddo
  end subroutine ac_porous_drag_res

!===============================================================================
! LTNE scalar pair (placeholder; implemented when validating porous_ltne_1d)
!===============================================================================
  subroutine ac_porous_ltne(nMesh, mBlock)
   use Global_var
   implicit none
   integer :: nMesh, mBlock
   if(my_id .eq. 0) print*, ' AC-porous LTNE scalar solver: to be implemented'
  end subroutine ac_porous_ltne

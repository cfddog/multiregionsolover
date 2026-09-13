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
!  dedicated porous coolant inlet (overrides the low-speed LS_U_in if set)
   if(Porous_U_in /= 0.d0 .or. Porous_V_in /= 0.d0 .or. Porous_W_in /= 0.d0) then
     uin_x = Porous_U_in; uin_y = Porous_V_in; uin_z = Porous_W_in
   endif
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
      ' beta=', beta, ' CFL=', AC_CFL, ' flux=', AC_Flux, ' recon=', AC_Recon, &
      ' lim=', AC_Limiter

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
! AC-porous LTNE two-temperature scalar pair (Tf / Ts), solved after the AC
! flow has converged.
!
! The face mass fluxes are reconstructed from the converged AC velocity field
! and the coupled steady system is iterated with the SAME finite-volume
! discretisation as the validated SIMPLE porous LTNE path (sub_porous.f90):
!   fluid :  rho*cp u.grad(Tf)  = div(k_f,eff grad Tf) + hv*(Ts - Tf)
!   frame :  div(k_s,eff grad Ts) + hv*(Tf - Ts) = 0
!   k_f,eff = LS_k*eps   ;  k_s,eff = (1-eps)*solid_k
! including the porous boundary ghost conventions (Tf inlet/outlet/wall ghosts
! and the solid_bc.inp frame BCs: isothermal / heat-flux / Robin extrapolation).
!===============================================================================
  subroutine ac_porous_ltne(nMesh, mBlock)
   use Global_var
   use const_var
   use porous_work
   implicit none
   integer :: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, i,j,k, iter
   real(PRE_EC) :: rho, dTf, dTs, Tfmin, Tfmax, Tsmin, Tsmax
   real(PRE_EC), allocatable :: Tfold(:,:,:), Tsold(:,:,:)

   interface
     subroutine porous_boundary_block(nMesh, mBlock)
       use precision_EC
       implicit none
       integer, intent(in) :: nMesh, mBlock
     end subroutine porous_boundary_block
     subroutine porous_energy(nMesh, mBlock)
       use precision_EC
       implicit none
       integer, intent(in) :: nMesh, mBlock
     end subroutine porous_energy
     subroutine porous_energy_Ts(nMesh, mBlock)
       use precision_EC
       implicit none
       integer, intent(in) :: nMesh, mBlock
     end subroutine porous_energy_Ts
   end interface

   B => Mesh(nMesh)%Block(mBlock)
   nx = B%nx; ny = B%ny; nz = B%nz
   rho = max(LS_rho, 1.d-20)

!  boundary ghosts + boundary-face mass fluxes (recomputed inside the sweep
!  loop too, so the frame thermal BCs track the evolving Ts)
   call porous_boundary_block(nMesh, mBlock)

!  interior face mass fluxes (kg/s) from the converged AC velocity field.
!  Fi/Fj/Fk carry the flux in the +coordinate direction, exactly as in
!  porous_face_flux: Fi(i,j,k) = rho*Si(i,j,k)*u on the face between cells
!  i-1 and i.  (For a converged divergence-free AC field the plain average of
!  the two adjacent cell velocities is the consistent face value.)
   do k = 1, nz-1
   do j = 1, ny-1
   do i = 2, nx-1
     Fi(i,j,k) = rho*B%Si(i,j,k)*0.5d0*(B%U(2,i-1,j,k)+B%U(2,i,j,k))
   enddo; enddo; enddo
   do k = 1, nz-1
   do j = 2, ny-1
   do i = 1, nx-1
     Fj(i,j,k) = rho*B%Sj(i,j,k)*0.5d0*(B%U(3,i,j-1,k)+B%U(3,i,j,k))
   enddo; enddo; enddo
   do k = 2, nz-1
   do j = 1, ny-1
   do i = 1, nx-1
     Fk(i,j,k) = rho*B%Sk(i,j,k)*0.5d0*(B%U(4,i,j,k-1)+B%U(4,i,j,k))
   enddo; enddo; enddo

!  coupled Gauss-Seidel sweeps (INNER sweeps inside porous_energy / _Ts);
!  stop on the change of both temperatures over a 100-sweep window
   allocate(Tfold(nx,ny,nz), Tsold(nx,ny,nz))
   Tfold = 0.d0; Tsold = 0.d0
   do iter = 1, max(Porous_Max_Iter, 2000)
     if(mod(iter,100) .eq. 1) then
       do k = 1, nz-1; do j = 1, ny-1; do i = 1, nx-1
         Tfold(i,j,k) = B%U(5,i,j,k); Tsold(i,j,k) = B%Ts(i,j,k)
       enddo; enddo; enddo
     endif
     call porous_boundary_block(nMesh, mBlock)
     call porous_energy(nMesh, mBlock)
     call porous_energy_Ts(nMesh, mBlock)
     if(mod(iter,100) .eq. 0) then
       dTf = 0.d0; dTs = 0.d0
       do k = 1, nz-1; do j = 1, ny-1; do i = 1, nx-1
         dTf = max(dTf, abs(B%U(5,i,j,k)-Tfold(i,j,k)))
         dTs = max(dTs, abs(B%Ts(i,j,k)-Tsold(i,j,k)))
       enddo; enddo; enddo
       if(my_id .eq. 0 .and. mod(iter,2000) .eq. 0) then
         Tfmin = minval(B%U(5,1:nx-1,1:ny-1,1:nz-1))
         Tfmax = maxval(B%U(5,1:nx-1,1:ny-1,1:nz-1))
         Tsmin = minval(B%Ts(1:nx-1,1:ny-1,1:nz-1))
         Tsmax = maxval(B%Ts(1:nx-1,1:ny-1,1:nz-1))
         print*, '  AC-porous LTNE iter', iter, ' dTf=', dTf, ' dTs=', dTs, &
                 ' Tf[', Tfmin, ',', Tfmax, ']  Ts[', Tsmin, ',', Tsmax, ']'
       endif
       if(dTf .lt. 1.d-3 .and. dTs .lt. 1.d-3) exit
     endif
   enddo
   deallocate(Tfold, Tsold)
   if(my_id .eq. 0) print*, ' AC-porous LTNE block', mBlock, ': sweeps =', iter, &
       ' final dTf=', dTf, ' dTs=', dTs
  end subroutine ac_porous_ltne


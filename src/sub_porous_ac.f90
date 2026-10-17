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
   use region_stat_mod
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock
   real(PRE_EC) :: Sfac, Sfac1
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, iter, i,j,k, isdual
   integer :: imax_s, jmax_s, kmax_s
   real(PRE_EC) :: uin_x, uin_y, uin_z, Uref, beta, U2, L
   real(PRE_EC) :: rq, rm
   real(PRE_EC) :: rmax_s, rrms_s, rq_s, rm_s
   logical :: have_pbc

   B => Mesh(nMesh)%Block(mBlock)
   nx = B%nx; ny = B%ny; nz = B%nz
   L  = max(Lscale, 1.d-30)

   if(.not. allocated(XW) .or. ac_nxw/=nx .or. ac_nyw/=ny .or. ac_nzw/=nz &
      .or. ac_lap /= LAP) then
     if(allocated(XW)) deallocate(XW,RAC,DU4,DTAC,SIG,dragc,ac_wallface,QWF)
     allocate( XW(4, 1-LAP:nx+LAP-1, 1-LAP:ny+LAP-1, 1-LAP:nz+LAP-1) )
     allocate( RAC(4, nx, ny, nz), DU4(4, nx, ny, nz) )
     allocate( DTAC(nx, ny, nz), SIG(nx, ny, nz), dragc(nx, ny, nz) )
     allocate( ac_wallface(nx, ny, nz, 3), QWF(nx, ny, nz) )
     ac_nxw=nx; ac_nyw=ny; ac_nzw=nz; ac_lap=LAP
     ac_wallface = 0
   endif
   dragc = 0.d0
   call ac_check_controls()

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
     call ac_wall_pressure(nMesh, mBlock)   ! wall-face pressure (needs fresh XW)
     call ac_boundary_flux(nMesh, mBlock, beta, uin_x, uin_y, uin_z)
     call ac_internal_flux(nMesh, mBlock, beta)
     if(If_viscous .eq. 1 .and. LS_mu .gt. 0.d0) call ac_viscous_res(nMesh, mBlock)
     call ac_porous_drag_res(nMesh, mBlock)              ! sets dragc + RAC(2:4)
     if(isdual .eq. 1) call ac_dual_res(nMesh, mBlock, Sfac)
     call ac_compute_dt(nMesh, mBlock, beta)
     call ac_res_norm(nMesh, mBlock, beta, U2, L, rq, rm)
!    (the separate iter==1 print was removed 2026-09-26; the line below covers it)
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
!    eps <= EPS_SOLID_POROUS: the impermeable-solid limit -- freeze the state at
!    u = 0 (interior + ghost rows).  Together with dragc = 1e30 from
!    ac_porous_drag_res this makes the block a pure conduction region, i.e. the
!    flow field is independent of the drag model and equivalent to a
!    BLOCK_SOLID for the LTNE frame temperature.
     if(B%porous_eps .le. EPS_SOLID_POROUS) then
       do k = 1-LAP, nz+LAP-1
       do j = 1-LAP, ny+LAP-1
       do i = 1-LAP, nx+LAP-1
         B%U(2,i,j,k) = 0.d0; B%U(3,i,j,k) = 0.d0; B%U(4,i,j,k) = 0.d0
       enddo; enddo; enddo
     endif

     if(.not. have_pbc) call ac_zero_mean_p(nMesh, mBlock)
     if(iter .eq. 1 .or. mod(iter, AC_Print) .eq. 0) then
       if(my_id .eq. 0) print*, '  AC-porous iter', iter, ' res_q=', rq, ' res_m=', rm
     endif
     if(iter .ge. 20 .and. rq .lt. AC_Tol .and. rm .lt. AC_Tol) exit
   enddo
!  register the AC residual of this block for the staggered outer-loop stop
!  test (Iflag_Couple_Stop=1); see sub_lowspeed_ac.f90 for the convention.
   Couple_Res_AC_q = max(Couple_Res_AC_q, rq)
   Couple_Res_AC_m = max(Couple_Res_AC_m, rm)
!  区域汇总行：max / 平均 rms / max 位置（-iter: final=总伪时间步数）
   if(my_id .eq. 0) then
     call ac_res_stat(nMesh, mBlock, beta, U2, L, rmax_s, imax_s, jmax_s, kmax_s, &
                      rrms_s, rq_s, rm_s)
     call region_res_stat(nMesh, mBlock, 3, -iter, rmax_s, imax_s, jmax_s, kmax_s, &
                          rrms_s, rq, rm)
   endif

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
!    eps -> 0 (impermeable solid): the Ergun Darcy term eps^2*mu/K would
!    formally give 0/1e-30 = 0, i.e. NO drag, and the block would behave like a
!    free fluid.  Renormalise eps to EPS_SOLID_POROUS (K -> 1e-30, c_drag ->
!    1.8e9 1/s) and mark the diagonal as drag-dominated; the caller freezes the
!    state at u = 0, so eps = 0 is the exact solid limit and stays continuous
!    with eps ~ 1e-8 (both give c_drag ~ 1e9 1/s).
     if(epsp .le. EPS_SOLID_POROUS) then
       dragc(i,j,k) = 1.d30
       cycle
     endif
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
   use region_stat_mod
   use const_var
   use porous_work
   implicit none
   integer :: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, i,j,k, iter
   integer :: imax_f, jmax_f, kmax_f, imax_g, jmax_g, kmax_g, ncell_t
   real(PRE_EC) :: rho, dTf, dTs, Tfmin, Tfmax, Tsmin, Tsmax
   real(PRE_EC) :: rss_f, rss_g, rcell_f, rcell_g, rrms_f, rrms_g
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
   dTf = 0.d0; dTs = 0.d0
   imax_f = 1; jmax_f = 1; kmax_f = 1
   imax_g = 1; jmax_g = 1; kmax_g = 1
   rrms_f = 0.d0; rrms_g = 0.d0
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
       rss_f = 0.d0; rss_g = 0.d0
       ncell_t = 0
       do k = 1, nz-1; do j = 1, ny-1; do i = 1, nx-1
         rcell_f = abs(B%U(5,i,j,k)-Tfold(i,j,k))
         rcell_g = abs(B%Ts(i,j,k)-Tsold(i,j,k))
         if(rcell_f > dTf) then      ! 流体/骨架各自的最大变化位置
           dTf = rcell_f; imax_f = i; jmax_f = j; kmax_f = k
         endif
         if(rcell_g > dTs) then
           dTs = rcell_g; imax_g = i; jmax_g = j; kmax_g = k
         endif
         rss_f = rss_f + rcell_f*rcell_f
         rss_g = rss_g + rcell_g*rcell_g
         ncell_t = ncell_t + 1
       enddo; enddo; enddo
       rrms_f = sqrt(rss_f/dble(max(ncell_t,1)))
       rrms_g = sqrt(rss_g/dble(max(ncell_t,1)))
!      (the periodic ' AC-porous LTNE iter ...' line was removed 2026-09-26)
       if(dTf .lt. 1.d-3 .and. dTs .lt. 1.d-3) exit
     endif
   enddo
   deallocate(Tfold, Tsold)
!  区域汇总行：100 步窗口变化的 max / 平均 rms / max 位置（-iter: final=总扫掠数）
   if(my_id .eq. 0) then
     call region_res_stat(nMesh, mBlock, 7, -iter, dTf, imax_f, jmax_f, kmax_f, &
                          rrms_f, dTf, 0.d0)
     call region_res_stat(nMesh, mBlock, 8, -iter, dTs, imax_g, jmax_g, kmax_g, &
                          rrms_g, dTs, 0.d0)
   endif
  end subroutine ac_porous_ltne


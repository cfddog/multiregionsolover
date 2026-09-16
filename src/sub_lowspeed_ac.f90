!===============================================================================
! Artificial-compressibility (AC) low-speed solver for BLOCK_LOWSPEED blocks.
! Selected by  LS_Algorithm = 3  (1=SIMPLE, 2=SIMPLEC keep the old paths).
!   state: q = p/rho (kinematic), u,v,w ;  T in B%U(5) passive (later).
!   pseudo-time:  dq/dtau + beta*div(u)=0 ;  du/dtau + div(uu+qI)= nu Lap(u)
!   beta = AC_beta * U_ref^2 (U_ref local velocity scale).
! Cell-centred FV on the SAME structured multi-block mesh & geometry as the
! compressible solver; internal faces: Rusanov-type AC flux (v1 1st order);
! physical boundary faces: explicit BC flux; viscous: orthogonal-grid
! Laplacian through face gradients (nu=LS_mu/rho); time: pseudo-time implicit
! LU-SGS (steady) or dual-time BDF2 LU-SGS (Time_Dual_LU_SGS).
!===============================================================================

  module lows_ac_work
   use precision_EC
   implicit none
   real(PRE_EC), allocatable, dimension(:,:,:,:) :: XW    ! (4,[q,u,v,w], ghost)
   real(PRE_EC), allocatable, dimension(:,:,:,:) :: RAC   ! (4, cells)
   real(PRE_EC), allocatable, dimension(:,:,:,:) :: DU4   ! delta per sweep
   real(PRE_EC), allocatable, dimension(:,:,:)   :: DTAC  ! pseudo time step
   real(PRE_EC), allocatable, dimension(:,:,:)   :: SIG   ! spectral sum
   real(PRE_EC), allocatable, dimension(:,:,:)   :: dragc ! porous drag coeff (1/s, per cell)
   integer :: ac_nxw=0, ac_nyw=0, ac_nzw=0, ac_lap=0
!  ---- same-class (LOWSPEED<->LOWSPEED) block-interface flux registry --------
!  The OWNER of an interface (the block with the SMALLER block number) computes
!  the interface flux ONCE per pseudo-step and stores it in FSH indexed by its
!  own cell adjacent to the interface; every other block on that interface only
!  READS it and applies the opposite contribution to its own cell.  The flux is
!  therefore single valued and the interface is conservative to round-off.
!  Evaluating the flux on both sides from their own ghost buffers is NOT
!  conservative: the numerical dissipation term 0.5*lam*(qR-qL) has lam even in
!  the face normal, so the two one-sided fluxes differ by lam*(qB-qA) instead of
!  cancelling -- a spurious interface mass/momentum source.
!  If a pair cannot be resolved unambiguously (id = 0) the historical two-sided
!  flux is used, so the change is backward compatible.
   integer, parameter :: AC_IFC_MAX = 64
   integer :: ac_ifc_n = 0
   integer :: ac_ifc_mesh = -1
!  NOTE (known limitation, 2026-09-15): the registry is built for ONE mesh
!  (ac_ifc_mesh) and is rebuilt when the AC solver is called for another mesh,
!  which also zeroes FSH.  Same-class (BC_Inner) interfaces always live inside a
!  single mesh, so this is exact as long as the blocks on either side of an
!  interface are solved without an intervening switch to another mesh (all
!  current cases: 2-block / 4-block channel, BJ, bl_cht, fluid_solid).  If a
!  future run interleaves same-class AC interfaces across meshes, build the
!  registry once for all meshes (loop nMesh=1..Num_Mesh, store the mesh number
!  per interface and look the owner up in Mesh(ac_ifc_msh(id))).
   integer :: ac_ifc_own(AC_IFC_MAX) = 0
   integer :: ac_ifc_ks (AC_IFC_MAX) = 0
   integer :: ac_ifc_nb (AC_IFC_MAX) = 0
   real(PRE_EC), allocatable, dimension(:,:,:,:,:) :: FSH  ! (4,nx,ny,nz,nifc)
!  ---- wall-face map for the near-wall reconstruction (AC_WallRecon) ---------
!  ac_wallface(i,j,k,dir), dir = 1,2,3 for the i,j,k index direction:
!    1 -> the face on the LOW-index side of cell (i,j,k) is a wall (BC_Wall) or
!         a symmetry plane (BC_Symmetry),
!    2 -> the face on the HIGH-index side is,
!    0 -> neither.
!  Filled from the block subfaces by ac_fill_ghost, so partial boundary faces are
!  resolved cell by cell.  Used by ac_face_flux to de-bias the tangential states
!  of the first interior face next to a wall/symmetry plane.
   integer, allocatable, dimension(:,:,:,:) :: ac_wallface   ! (nx,ny,nz,3)
!  ---- wall-face pressure (AC_WallP) ----------------------------------------
!  QWF(i,j,k) = pressure q = p/rho to be used on the wall/symmetry face of cell
!  (i,j,k): 1.5*q_1 - 0.5*q_2 (linear reconstruction through the face) when
!  AC_WallP=1, else the cell value q_1.  Filled by ac_wall_pressure.
   real(PRE_EC), allocatable, dimension(:,:,:) :: QWF
!  ---- AC controls: validated/clamped once per run ---------------------------
   logical :: ac_ctl_checked = .false.
!  ---- optional diagnostic for the interface sharing (default OFF) -----------
!  ac_dbg=1 prints, for one face of every same-class interface (the first 3
!  pseudo-steps and then every 500), the owner cell / ghost state, the stored
!  flux and the neighbour's mapped read slot, plus the (mBlock,ksub) ->
!  (id,is_owner) resolution of the first few subfaces.  Used to verify that the
!  neighbour reads exactly the slot the owner wrote (2-block channel 2026-09-15).
   integer :: ac_dbg = 0
   integer :: ac_dbg_own = 0, ac_dbg_nb = 0, ac_dbg_find = 0
  end module lows_ac_work

!===============================================================================
! AC driver
!===============================================================================
  subroutine lowspeed_ac_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_var
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock
   real(PRE_EC) :: Sfac, Sfac1
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, iter, i,j,k, isdual
   real(PRE_EC) :: uin_x, uin_y, uin_z, Uref, beta, U2, L
   real(PRE_EC) :: rq, rm, rq0, rm0
   logical :: have_pbc
   interface
     subroutine update_buffer_onemesh(nMesh)
       use precision_EC
       integer:: nMesh
     end subroutine update_buffer_onemesh
   end interface

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
     dragc = 0.d0
     ac_wallface = 0
   endif
   dragc = 0.d0   ! clear residual drag coefficient (porous AC sets it per iteration)
   call ac_check_controls()

   call lowspeed_inlet_velocity(B, uin_x, uin_y, uin_z)
   Uref = sqrt(uin_x*uin_x + uin_y*uin_y + uin_z*uin_z)
   if(Uref < 1.d-12) Uref = max(abs(LS_U_lid), 1.d-3)      ! lid-driven cavity etc.
   if(Uref < 1.d-12) Uref = sqrt(abs(LS_P_in-LS_P_out)/max(LS_rho,1.d-20))
   if(Uref < 1.d-12) Uref = 1.d0
   beta = max(AC_beta, 1.d-8) * Uref * Uref
   U2   = max(Uref*Uref, 1.d-20)

   isdual = 0
   if(Time_Method .eq. Time_Dual_LU_SGS) isdual = 1

   call ac_have_pressure_bc(nMesh, mBlock, have_pbc)

   if(my_id .eq. 0) print*, ' AC solver (LS_Algorithm=3) block', mBlock, &
      ' nx,ny,nz=',nx,ny,nz,' beta=',beta,' CFL=',AC_CFL,' dual=',isdual, &
      ' flux=',AC_Flux,' recon=',AC_Recon,' lim=',AC_Limiter

   rq0 = 0.d0; rm0 = 0.d0
   do iter = 1, AC_Max_Iter
     call ac_fill_ghost(nMesh, mBlock, uin_x, uin_y, uin_z)
!    Same-class block interfaces (BC_Inner, -1): refresh the ghost
!    buffers from the neighbours on EVERY pseudo-step, not only once per
!    outer step.  The lagged interface is what destabilises the coupled
!    AC march.  Restricted to a single process: a per-pseudo-step MPI
!    exchange would require both sides to iterate the same number of
!    times, which they do not (each block exits on AC_Tol).
     if(Total_proc .eq. 1) call update_buffer_onemesh(nMesh)
     call ac_load_state(nMesh, mBlock)
     call ac_wall_pressure(nMesh, mBlock)   ! wall-face pressure (needs fresh XW)
     call ac_boundary_flux(nMesh, mBlock, beta, uin_x, uin_y, uin_z)
     call ac_internal_flux(nMesh, mBlock, beta)
     if(If_viscous .eq. 1 .and. LS_mu .gt. 0.d0) call ac_viscous_res(nMesh, mBlock)
     if(isdual .eq. 1) call ac_dual_res(nMesh, mBlock, Sfac)
     call ac_compute_dt(nMesh, mBlock, beta)
     call ac_res_norm(nMesh, mBlock, beta, U2, L, rq, rm)
     if(iter .eq. 1) then
       rq0 = max(rq, 1.d-30); rm0 = max(rm, 1.d-30)
       if(my_id .eq. 0) print*, '  AC iter', iter, ' res_q=', rq, ' res_m=', rm
     endif
     call ac_lusgs_sweep(nMesh, mBlock, beta, isdual, Sfac1)
!    write back primitive state
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
       if(my_id .eq. 0) print*, '  AC iter', iter, ' res_q=', rq, ' res_m=', rm
     endif
     if(iter .ge. 20 .and. rq .lt. AC_Tol .and. rm .lt. AC_Tol) exit
   enddo
   if(my_id .eq. 0) print*, ' AC solver block', mBlock, ': iterations =', iter, &
       ' final res_q=', rq, ' res_m=', rm
!   Temperature (passive scalar / conjugate): after the flow converges, relax
!   T (=U(5)) to steady state using the converged AC velocity field when the
!   fluid conducts heat (LS_k > 0).  Same discretisation as SIMPLE.
    if(LS_k .gt. 0.d0 .or. LS_T_wall .gt. 0.d0) call ac_lowspeed_energy(nMesh, mBlock)

  end subroutine lowspeed_ac_solver_one_block
!===============================================================================
! XW <- (q,u,v,w) from the primitive state (incl. ghost cells)
!===============================================================================
  subroutine ac_load_state(nMesh, mBlock)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   do k = 1-LAP, B%nz+LAP-1
   do j = 1-LAP, B%ny+LAP-1
   do i = 1-LAP, B%nx+LAP-1
     XW(1,i,j,k) = B%p(i,j,k)/LS_rho
     XW(2,i,j,k) = B%U(2,i,j,k)
     XW(3,i,j,k) = B%U(3,i,j,k)
     XW(4,i,j,k) = B%U(4,i,j,k)
   enddo; enddo; enddo
  end subroutine ac_load_state

!===============================================================================
! Has at least one physical boundary prescribing an absolute pressure level?
!===============================================================================
  subroutine ac_have_pressure_bc(nMesh, mBlock, flag)
   use Global_var
   use const_var
   implicit none
   integer :: nMesh, mBlock, ksub, mb
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   logical :: flag
   B => Mesh(nMesh)%Block(mBlock)
   flag = .false.
   do ksub = 1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle
     if(associated(B%bc_msg2)) then
       if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
     endif
     if(Bc%bc .eq. BC_Outflow .or. Bc%bc .eq. BC_LS_Outlet .or. &
        Bc%bc .eq. BC_Farfield) flag = .true.
     if((Bc%bc .eq. BC_Inflow .or. Bc%bc .eq. BC_LS_Inlet) .and. &
        LS_Inlet_Type .eq. 3) flag = .true.
   enddo
!  Coupled same-class blocks carry ONE pressure field: a pressure BC on any
!  low-speed block of the mesh anchors the level of all of them through the
!  interface ghost, so NO block may remove its own mean pressure.  Removing the
!  per-block mean of an unanchored inlet block injects a spurious pressure jump
!  across each interface and lets that pressure mode grow (block-1 style
!  "converges then diverges"); this is what killed the multi-block AC runs.
   if(.not. flag) then
     do mb = 1, Mesh(nMesh)%Num_Block
       if(mb .eq. mBlock) cycle
       if(Mesh(nMesh)%Block(mb)%Block_type .ne. BLOCK_LOWSPEED) cycle
       if(.not. associated(Mesh(nMesh)%Block(mb)%bc_msg)) cycle
       B => Mesh(nMesh)%Block(mb)
       do ksub = 1, B%subface
         Bc => B%bc_msg(ksub)
         if(is_interface_bc(Bc%bc)) cycle
         if(Bc%bc .eq. BC_Outflow .or. Bc%bc .eq. BC_LS_Outlet .or. &
            Bc%bc .eq. BC_Farfield) flag = .true.
         if((Bc%bc .eq. BC_Inflow .or. Bc%bc .eq. BC_LS_Inlet) .and. &
            LS_Inlet_Type .eq. 3) flag = .true.
       enddo
     enddo
     B => Mesh(nMesh)%Block(mBlock)
   endif
  end subroutine ac_have_pressure_bc

!===============================================================================
! Remove the spatial mean of p (reference-pressure anchoring, no pressure BC)
!===============================================================================
  subroutine ac_zero_mean_p(nMesh, mBlock)
   use Global_var
   implicit none
   integer :: nMesh, mBlock, i,j,k, ncell
   Type (Block_TYPE),pointer:: B
   real(PRE_EC) :: pmean
   B => Mesh(nMesh)%Block(mBlock)
   ncell = 0; pmean = 0.d0
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     pmean = pmean + B%p(i,j,k); ncell = ncell + 1
   enddo; enddo; enddo
   if(ncell .gt. 0) pmean = pmean / real(ncell)
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     B%p(i,j,k) = B%p(i,j,k) - pmean
   enddo; enddo; enddo
  end subroutine ac_zero_mean_p

!===============================================================================
! Rusanov flux for the AC system along unit normal (anx,any,anz)
! F/area = ( beta*un,  u*un+q*nx,  v*un+q*ny,  w*un+q*nz )
!
! DISSIPATION SPLIT BY FIELD (AC_MomDiss=1, default)
! The AC Jacobian has eigenvalues {un, un, un+c, un-c}, c = sqrt(un^2+beta): the
! pressure/mass pair travels at un+-c, the transverse velocity components at un.
! Using the single largest eigenvalue lam = max(|un|+c) for ALL four equations
! (the classic lumped Rusanov form) over-dissipates the momentum equations by
! O(c/|un|) >> 1 in a low-speed flow, i.e. it adds an artificial viscosity
! ~sqrt(beta)*dx to the momentum.  That is exactly the numerical "dead water"
! measured behind the inviscid cylinder (wall Cp plateaus at +0.39 instead of
! +1.0 for 0<theta<18 deg although the exact solution does not separate, and the
! velocity on the rear symmetry ray is only 0.17/0.46 at r=1.7/2.1 R instead of
! 0.67/0.78).  AC_MomDiss=1 keeps |un|+c on the mass equation (this is what
! couples and damps the collocated pressure field -- the AC analogue of the
! Rhie-Chow damping) and uses |un| on the three momentum equations, i.e. plain
! convective upwinding.  AC_MomDiss=0 restores the previous behaviour exactly.
!===============================================================================
  subroutine ac_flux_rusanov(qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta,F)
   use precision_EC
   use Global_var
   implicit none
   real(PRE_EC),intent(in) :: qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta
   real(PRE_EC),intent(out):: F(4)
   real(PRE_EC) :: unL,unR,cL,cR,lam,lamm
   unL = uL*anx + vL*any + wL*anz
   unR = uR*anx + vR*any + wR*anz
   cL  = sqrt(unL*unL + beta)
   cR  = sqrt(unR*unR + beta)
   lam = max(abs(unL)+cL, abs(unR)+cR, 1.d-12)
!  AC_MomDiss: momentum dissipation coefficient
!    0 = lumped   : lam = max(|un|+c)            (DEFAULT, previously validated)
!    1 = blended  : lam = max(|un|) + AC_MomFrac*max(c)   (AC_MomFrac default 0.2)
!    2 = pure     : lam = max(|un|)              (unstable at stagnation points)
!  The AC eigenvalues are {un, un, un+-c}: the mass equation keeps the full
!  |un|+c coefficient (this is what couples and damps the collocated pressure
!  field) while the momentum equations only need the convective scale |un|.
!  Option 0 applies the full |un|+c to the momentum as well, i.e. it adds an
!  artificial viscosity ~sqrt(beta)*dx there.
!  MEASURED TRADE-OFF (2026-09-17, docs/工作日志.md):
!    option 1 IMPROVES wall-bounded / viscous flows -- channel u_max error
!    0.303% -> 0.153% and profile RMS 0.184% -> 0.130% (same iteration count),
!    cylinder Re=40 Cd 1.5706 -> 1.5589 (Rogers 1.549) --
!    but it DEGRADES the fully inviscid cylinder (max|dCp| 0.5413 -> 0.8712 with
!    AC_MomFrac=0.2, 0.6209 with 0.5), where the numerical dissipation is the only
!    mechanism damping the collocated modes of the (physically non-dissipative)
!    Euler solution.  Option 0 therefore stays the default so that the documented
!    validation results are unchanged; use option 1 per case for viscous/internal
!    flows and re-validate.
!  Option 2 is NOT stable anywhere tested: at a stagnation point un -> 0 removes
!  the momentum dissipation and the collocated coupling diverges (inviscid
!  cylinder: res_q 0.15 at 6k, 0.24 at 14k, 0.40 at 16k pseudo steps).
   lamm = lam
   if(AC_MomDiss .eq. 1) then
     lamm = max(abs(unL), abs(unR), 1.d-12) + AC_MomFrac*max(cL, cR)
   else if(AC_MomDiss .eq. 2) then
     lamm = max(abs(unL), abs(unR), 1.d-12)
   endif
   F(1) = 0.5d0*beta*(unL+unR) - 0.5d0*lam*(qR-qL)
   F(2) = 0.5d0*((uL*unL+qL*anx) + (uR*unR+qR*anx)) - 0.5d0*lamm*(uR-uL)
   F(3) = 0.5d0*((vL*unL+qL*any) + (vR*unR+qR*any)) - 0.5d0*lamm*(vR-vL)
   F(4) = 0.5d0*((wL*unL+qL*anz) + (wR*unR+qR*anz)) - 0.5d0*lamm*(wR-wL)
  end subroutine ac_flux_rusanov

!===============================================================================
! Steger-Warming (characteristic) flux splitting for the AC system, selected by
! AC_Flux=2.  Along a unit normal n the Jacobian d(F.n)/dW (W=(q,u,v,w)) has
! eigenvalues {un, un, un+c, un-c}, c=sqrt(un^2+beta); the right eigenvectors
! are (0,0,1,0),(0,0,0,1),(beta,un+c,0,0),(beta,un-c,0,0) and the left ones
! (0,0,1,0),(0,0,0,1),(1,un+c,0,0)/(beta+(un+c)^2),(1,un-c,0,0)/(beta+(un-c)^2).
! The split flux is  F+(W)=0.5*(F(W)+|A|W),  F-(W)=0.5*(F(W)-|A|W)  with
! |A|W obtained by projecting W on the characteristic variables; the transverse
! part collapses to  un^+- * (u - un*n).  Face flux = F+(WL) + F-(WR).
!===============================================================================
   subroutine ac_flux_sw(qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta,F)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta
    real(PRE_EC),intent(out):: F(4)
    real(PRE_EC) :: Fp(4),Fm(4)
    call ac_sw_half(qL,uL,vL,wL,anx,any,anz,beta, 1,Fp)
    call ac_sw_half(qR,uR,vR,wR,anx,any,anz,beta,-1,Fm)
    F(1)=Fp(1)+Fm(1); F(2)=Fp(2)+Fm(2)
    F(3)=Fp(3)+Fm(3); F(4)=Fp(4)+Fm(4)
   end subroutine ac_flux_sw

!  one side of the Steger-Warming split: sgn=+1 -> F+, sgn=-1 -> F-
   subroutine ac_sw_half(q,u,v,w,anx,any,anz,beta,sgn,fh)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: q,u,v,w,anx,any,anz,beta
    integer,intent(in) :: sgn
    real(PRE_EC),intent(out):: fh(4)
    real(PRE_EC) :: un,c,cp,cm,s3,s4,awq,awu,fn,unp,ut1,ut2,ut3
    un = u*anx + v*any + w*anz
    c  = sqrt(un*un + beta)
    cp = un + c
    cm = un - c
    s3 = abs(cp)/(beta+cp*cp)
    s4 = abs(cm)/(beta+cm*cm)
    awq = beta*( s3*(q+cp*un) + s4*(q+cm*un) )
    awu = s3*cp*(q+cp*un) + s4*cm*(q+cm*un)
    ut1 = u - un*anx
    ut2 = v - un*any
    ut3 = w - un*anz
    if(sgn .gt. 0) then
      fh(1) = 0.5d0*(beta*un + awq)
      fn    = 0.5d0*(un*un + q + awu)
      unp   = 0.5d0*(un + abs(un))
    else
      fh(1) = 0.5d0*(beta*un - awq)
      fn    = 0.5d0*(un*un + q - awu)
      unp   = 0.5d0*(un - abs(un))
    endif
    fh(2) = fn*anx + unp*ut1
    fh(3) = fn*any + unp*ut2
    fh(4) = fn*anz + unp*ut3
   end subroutine ac_sw_half


!===============================================================================
! AUSM+ (AC_Flux=3) was REMOVED on 2026-09-17.
!  Reasons: (a) with the AC equations there is no physical sound speed, so the
!  Mach splitting is not defined consistently with the AC spectral radius
!  c=sqrt(un^2+beta) used by the LU-SGS implicit operator; (b) the plain AUSM+
!  splitting implemented here carries no low-Mach pressure dissipation (no
!  AUSM+-up style p_u term), which is exactly the mechanism that stabilises the
!  collocated AC pressure-velocity coupling; (c) empirically it diverged to NaN
!  within 2000 steps on the inviscid cylinder while Rusanov stayed stable
!  (cases/cylinder_re40_half/inv_ausm/run_ausm.log, docs/工作日志.md 2026-09-16).
!  Supported inviscid fluxes: AC_Flux=1 (Rusanov/LLF, default) and AC_Flux=2
!  (Steger-Warming characteristic splitting).  Other values are clamped to 1 by
!  ac_check_controls() with a warning.
!===============================================================================

! b and the right cell c, using the 5-point stencil a,b,c,d,e.
! Returns the face-left state wL and the face-right state wR.
!===============================================================================
   subroutine ac_weno5_lr(a,b,c,d,e,wL,wR)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: a,b,c,d,e
    real(PRE_EC),intent(out):: wL,wR
    real(PRE_EC),parameter :: g0=0.1d0, g1=0.6d0, g2=0.3d0, epw=1.d-6
    real(PRE_EC) :: u0,u1,u2, b0,b1,b2, a0,a1,a2, s
    u0 = ( 1.d0/3.d0)*a - (7.d0/6.d0)*b + (11.d0/6.d0)*c
    u1 = -(1.d0/6.d0)*b + (5.d0/6.d0)*c + ( 1.d0/3.d0)*d
    u2 = ( 1.d0/3.d0)*c + (5.d0/6.d0)*d - ( 1.d0/6.d0)*e
    b0 = (13.d0/12.d0)*(a-2.d0*b+c)**2 + 0.25d0*(a-4.d0*b+3.d0*c)**2
    b1 = (13.d0/12.d0)*(b-2.d0*c+d)**2 + 0.25d0*(b-d)**2
    b2 = (13.d0/12.d0)*(c-2.d0*d+e)**2 + 0.25d0*(3.d0*c-4.d0*d+e)**2
    a0 = g0/(epw+b0)**2; a1 = g1/(epw+b1)**2; a2 = g2/(epw+b2)**2
    s  = a0+a1+a2
    wL = (a0*u0 + a1*u1 + a2*u2)/s
    u0 = ( 1.d0/3.d0)*e - (7.d0/6.d0)*d + (11.d0/6.d0)*c
    u1 = -(1.d0/6.d0)*d + (5.d0/6.d0)*c + ( 1.d0/3.d0)*b
    u2 = ( 1.d0/3.d0)*c + (5.d0/6.d0)*b - ( 1.d0/6.d0)*a
    b0 = (13.d0/12.d0)*(e-2.d0*d+c)**2 + 0.25d0*(e-4.d0*d+3.d0*c)**2
    b1 = (13.d0/12.d0)*(d-2.d0*c+b)**2 + 0.25d0*(d-b)**2
    b2 = (13.d0/12.d0)*(c-2.d0*b+a)**2 + 0.25d0*(3.d0*c-4.d0*b+a)**2
    a0 = g0/(epw+b0)**2; a1 = g1/(epw+b1)**2; a2 = g2/(epw+b2)**2
    s  = a0+a1+a2
    wR = (a0*u0 + a1*u1 + a2*u2)/s
   end subroutine ac_weno5_lr

!===============================================================================
! 3rd-order WENO (Jiang-Shu) reconstruction (more dissipative/robust than WENO5)
! for the face between the left cell c and the right cell d; stencil a,b,c,d,e.
!===============================================================================
   subroutine ac_weno3_lr(a,b,c,d,e,wL,wR)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: a,b,c,d,e
    real(PRE_EC),intent(out):: wL,wR
    real(PRE_EC),parameter :: g0=1.d0/3.d0, g1=2.d0/3.d0, epw=1.d-6
    real(PRE_EC) :: u0,u1, bt0,bt1, a0,a1, s
    u0 = -0.5d0*b + 1.5d0*c
    u1 =  0.5d0*c + 0.5d0*d
    bt0 = (c-b)**2
    bt1 = (d-c)**2
    a0 = g0/(epw+bt0)**2; a1 = g1/(epw+bt1)**2
    s  = a0+a1
    wL = (a0*u0 + a1*u1)/s
    u0 = -0.5d0*e + 1.5d0*d
    u1 =  0.5d0*d + 0.5d0*c
    bt0 = (e-d)**2
    bt1 = (d-c)**2
    a0 = g0/(epw+bt0)**2; a1 = g1/(epw+bt1)**2
    s  = a0+a1
    wR = (a0*u0 + a1*u1)/s
   end subroutine ac_weno3_lr

!===============================================================================
! MUSCL slope limiter (AC_Limiter: 1=van Leer, 2=minmod).  minmod is the most
! dissipative/monotone choice and gives the most robust MUSCL reconstruction.
!===============================================================================
   real(PRE_EC) function ac_slope(a, b)
    use precision_EC
    use Global_var
    implicit none
    real(PRE_EC),intent(in) :: a, b
    if(a*b .le. 0.d0) then
      ac_slope = 0.d0
    else if(AC_Limiter .eq. 2) then
      ac_slope = sign(min(abs(a), abs(b)), a)
    else
      ac_slope = 2.d0*a*b/(a+b+1.d-30)
    endif
   end function ac_slope

!===============================================================================
! Face reconstruction of the left/right states from the 5-point stencil:
!   a,b = two cells left of the face ; c = face-left cell ;
!   d = face-right cell ; e = one cell right of d.
!   AC_Recon = 1 -> 2nd-order MUSCL (limiter AC_Limiter) ; 2 -> WENO5 ; 3 -> WENO3
!===============================================================================
   subroutine ac_recon_face(a,b,c,d,e,wL,wR)
    use precision_EC
    use Global_var
    implicit none
    real(PRE_EC),intent(in) :: a,b,c,d,e
    real(PRE_EC),intent(out):: wL,wR
    real(PRE_EC), external :: ac_slope
    if(AC_Recon .eq. 2) then
      call ac_weno5_lr(a,b,c,d,e,wL,wR)
    else if(AC_Recon .eq. 3) then
      call ac_weno3_lr(a,b,c,d,e,wL,wR)
    else
      wL = c + 0.5d0*ac_slope(c-b, d-c)
      wR = d - 0.5d0*ac_slope(d-c, e-d)
    endif
   end subroutine ac_recon_face

!===============================================================================
! Flux on one interior face between cells L and R=L+e_dir (normal points L->R).
! Reconstructs (q,u,v,w) with ac_recon_face then applies the selected inviscid
! flux (AC_Flux: 1=Rusanov, 2=Steger-Warming).  RAC(L) -= F*A ; RAC(R) += F*A.
!
! NEAR-WALL / CORNER DE-BIASING (AC_WallRecon=1, default)
! On a wall or symmetry plane the ghost state is a mirror of the interior one, so
! for a SLIP wall (If_viscous=0) or a symmetry plane the TANGENTIAL velocity of
! the ghost equals the cell value.  The two differences used by the MUSCL/WENO
! reconstruction of the tangential components then have one zero member and the
! limiter returns a zero slope: the wall-side face state degenerates to the cell
! value (1st order) while the interior side is extrapolated.  That 1st-order bias
! in the first layer is what produced the 8% first-cell velocity deficit and the
! total-pressure loss measured on the inviscid cylinder, and it de-biases every
! corner where a wall meets a symmetry plane or another wall (rear stagnation of
! the half cylinder, T-junctions).
! For such faces the tangential components are reconstructed from a ONE-SIDED
! interior stencil: the outer stencil value (the ghost) is replaced by the linear
! continuation of the interior field, a_eff = 2*c - d  (e_eff = 2*d - c on the
! high side).  For a no-slip wall with a linear near-wall profile this is
! identical to the mirror-ghost stencil, so viscous cases are not degraded.
! The wall-NORMAL component and the pressure keep the mirror-ghost stencil, so
! the no-penetration parity and therefore the mass flux are unchanged.
!===============================================================================
   subroutine ac_face_flux(B, iL,jL,kL, iR,jR,kR, dir, beta)
    use Global_var
    use lows_ac_work
    implicit none
    Type (Block_TYPE),pointer:: B
    integer,intent(in):: iL,jL,kL,iR,jR,kR,dir
    real(PRE_EC),intent(in):: beta
    real(PRE_EC) :: FLX(4), A, nx1,ny1,nz1
    real(PRE_EC) :: FLX1(4), sb
    real(PRE_EC) :: qL,uL,vL,wL, qR,uR,vR,wR
    real(PRE_EC) :: uLw,vLw,wLw, uRw,vRw,wRw, dummy
    logical :: wallL, wallR
    integer :: im2,jm2,km2, im1,jm1,km1, ip2,jp2,kp2
    im2=iL; jm2=jL; km2=kL; im1=iL; jm1=jL; km1=kL; ip2=iR; jp2=jR; kp2=kR
    if(dir .eq. 1) then
      im2=iL-2; im1=iL-1; ip2=iR+1
      A=B%Si(iR,jR,kR); nx1=B%ni1(iR,jR,kR); ny1=B%ni2(iR,jR,kR); nz1=B%ni3(iR,jR,kR)
    else if(dir .eq. 2) then
      jm2=jL-2; jm1=jL-1; jp2=jR+1
      A=B%Sj(iR,jR,kR); nx1=B%nj1(iR,jR,kR); ny1=B%nj2(iR,jR,kR); nz1=B%nj3(iR,jR,kR)
    else
      km2=kL-2; km1=kL-1; kp2=kR+1
      A=B%Sk(iR,jR,kR); nx1=B%nk1(iR,jR,kR); ny1=B%nk2(iR,jR,kR); nz1=B%nk3(iR,jR,kR)
    endif
!   ---- is one side of this face a wall/symmetry plane? ---------------------
    wallL = .false.; wallR = .false.
    if(AC_WallRecon .eq. 1 .and. allocated(ac_wallface)) then
      if(dir .eq. 1) then
        if(iL .eq. 1)      wallL = (ac_wallface(iL,jL,kL,1) .eq. 1)
        if(iR .eq. B%nx-1) wallR = (ac_wallface(iR,jR,kR,1) .eq. 2)
      else if(dir .eq. 2) then
        if(jL .eq. 1)      wallL = (ac_wallface(iL,jL,kL,2) .eq. 1)
        if(jR .eq. B%ny-1) wallR = (ac_wallface(iR,jR,kR,2) .eq. 2)
      else
        if(kL .eq. 1)      wallL = (ac_wallface(iL,jL,kL,3) .eq. 1)
        if(kR .eq. B%nz-1) wallR = (ac_wallface(iR,jR,kR,3) .eq. 2)
      endif
    endif
    call ac_recon_face(XW(1,im2,jm2,km2),XW(1,im1,jm1,km1),XW(1,iL,jL,kL), &
                       XW(1,iR,jR,kR),XW(1,ip2,jp2,kp2), qL,qR)
    call ac_recon_face(XW(2,im2,jm2,km2),XW(2,im1,jm1,km1),XW(2,iL,jL,kL), &
                       XW(2,iR,jR,kR),XW(2,ip2,jp2,kp2), uL,uR)
    call ac_recon_face(XW(3,im2,jm2,km2),XW(3,im1,jm1,km1),XW(3,iL,jL,kL), &
                       XW(3,iR,jR,kR),XW(3,ip2,jp2,kp2), vL,vR)
    call ac_recon_face(XW(4,im2,jm2,km2),XW(4,im1,jm1,km1),XW(4,iL,jL,kL), &
                       XW(4,iR,jR,kR),XW(4,ip2,jp2,kp2), wL,wR)
!   ---- wall side: one-sided tangential reconstruction -----------------------
    if(wallL) then
      call ac_recon_face(2.d0*XW(2,iL,jL,kL)-XW(2,iR,jR,kR),XW(2,im1,jm1,km1), &
           XW(2,iL,jL,kL),XW(2,iR,jR,kR),XW(2,ip2,jp2,kp2), uLw,dummy)
      call ac_recon_face(2.d0*XW(3,iL,jL,kL)-XW(3,iR,jR,kR),XW(3,im1,jm1,km1), &
           XW(3,iL,jL,kL),XW(3,iR,jR,kR),XW(3,ip2,jp2,kp2), vLw,dummy)
      call ac_recon_face(2.d0*XW(4,iL,jL,kL)-XW(4,iR,jR,kR),XW(4,im1,jm1,km1), &
           XW(4,iL,jL,kL),XW(4,iR,jR,kR),XW(4,ip2,jp2,kp2), wLw,dummy)
      call ac_tang_merge(uL,vL,wL, uLw,vLw,wLw, nx1,ny1,nz1)
    endif
    if(wallR) then
      call ac_recon_face(XW(2,im2,jm2,km2),XW(2,im1,jm1,km1),XW(2,iL,jL,kL), &
           XW(2,iR,jR,kR),2.d0*XW(2,iR,jR,kR)-XW(2,iL,jL,kL), dummy,uRw)
      call ac_recon_face(XW(3,im2,jm2,km2),XW(3,im1,jm1,km1),XW(3,iL,jL,kL), &
           XW(3,iR,jR,kR),2.d0*XW(3,iR,jR,kR)-XW(3,iL,jL,kL), dummy,vRw)
      call ac_recon_face(XW(4,im2,jm2,km2),XW(4,im1,jm1,km1),XW(4,iL,jL,kL), &
           XW(4,iR,jR,kR),2.d0*XW(4,iR,jR,kR)-XW(4,iL,jL,kL), dummy,wRw)
      call ac_tang_merge(uR,vR,wR, uRw,vRw,wRw, nx1,ny1,nz1)
    endif
    if(AC_Flux .eq. 2) then
      call ac_flux_sw(qL,uL,vL,wL, qR,uR,vR,wR, nx1,ny1,nz1, beta, FLX)
    else
      call ac_flux_rusanov(qL,uL,vL,wL, qR,uR,vR,wR, nx1,ny1,nz1, beta, FLX)
    endif
!   WENO5 reconstruction removes the low-order Rusanov dissipation that damps
!   the odd-even (checkerboard) pressure-velocity mode of the collocated AC
!   scheme, which makes it unstable.  AC_WenoBlend>0 restores a controlled
!   fraction (0..1) of the first-order (cell-value) Rusanov flux as a safety
!   dissipation while keeping the 5th-order reconstruction for the rest.
    if(AC_Recon .ge. 2 .and. AC_WenoBlend .gt. 0.d0) then
      call ac_flux_rusanov(XW(1,iL,jL,kL),XW(2,iL,jL,kL),XW(3,iL,jL,kL),XW(4,iL,jL,kL), &
                           XW(1,iR,jR,kR),XW(2,iR,jR,kR),XW(3,iR,jR,kR),XW(4,iR,jR,kR), &
                           nx1,ny1,nz1, beta, FLX1)
      sb = min(max(AC_WenoBlend, 0.d0), 1.d0)
      FLX(1:4) = (1.d0-sb)*FLX(1:4) + sb*FLX1(1:4)
    endif
    RAC(1:4,iL,jL,kL) = RAC(1:4,iL,jL,kL) - FLX(1:4)*A
    RAC(1:4,iR,jR,kR) = RAC(1:4,iR,jR,kR) + FLX(1:4)*A
   end subroutine ac_face_flux

!===============================================================================
! Merge a wall-side face state: the TANGENTIAL part comes from the one-sided
! interior reconstruction (uw,vw,ww), the WALL-NORMAL part from the mirror-ghost
! based state (u,v,w).  The normal part carries the no-penetration parity and the
! mass flux, so only the tangential momentum is de-biased (AC_WallRecon).
!===============================================================================
   subroutine ac_tang_merge(u,v,w, uw,vw,ww, nx1,ny1,nz1)
    use precision_EC
    implicit none
    real(PRE_EC),intent(inout) :: u,v,w
    real(PRE_EC),intent(in) :: uw,vw,ww, nx1,ny1,nz1
    real(PRE_EC) :: un, unw
    un  = u *nx1 + v *ny1 + w *nz1
    unw = uw*nx1 + vw*ny1 + ww*nz1
    u = (uw - unw*nx1) + un*nx1
    v = (vw - unw*ny1) + un*ny1
    w = (ww - unw*nz1) + un*nz1
   end subroutine ac_tang_merge

!===============================================================================
! Validate / clamp the AC controls once per run: a loud fallback instead of a
! silently wrong setting.  AC_Flux=3 (AUSM+) was removed on 2026-09-17 (see the
! note above ac_flux_rusanov); AC_Recon/AC_Limiter/AC_WenoBlend/AC_WallRecon are
! clamped to their documented ranges.
!===============================================================================
   subroutine ac_check_controls()
    use Global_var
    use lows_ac_work
    implicit none
    if(ac_ctl_checked) return
    ac_ctl_checked = .true.
    if(AC_Flux .ne. 1 .and. AC_Flux .ne. 2) then
      if(my_id .eq. 0) print*, ' WARNING: AC_Flux=',AC_Flux, &
        ' unsupported (1=Rusanov, 2=Steger-Warming; 3=AUSM+ removed 2026-09-17)' &
        ,' -> using AC_Flux=1'
      AC_Flux = 1
    endif
    if(AC_Recon .lt. 1 .or. AC_Recon .gt. 3) then
      if(my_id .eq. 0) print*, ' WARNING: AC_Recon=',AC_Recon, &
        ' out of range (1=MUSCL, 2=WENO5, 3=WENO3) -> using AC_Recon=1'
      AC_Recon = 1
    endif
    if(AC_Limiter .lt. 1 .or. AC_Limiter .gt. 2) AC_Limiter = 1
    AC_WenoBlend = min(max(AC_WenoBlend, 0.d0), 1.d0)
    if(AC_WallRecon .lt. 0 .or. AC_WallRecon .gt. 1) AC_WallRecon = 1
    if(AC_WallP .lt. 0 .or. AC_WallP .gt. 1) AC_WallP = 1
    if(AC_MomDiss .lt. 0 .or. AC_MomDiss .gt. 2) AC_MomDiss = 1
    AC_MomFrac = min(max(AC_MomFrac, 0.d0), 1.d0)
    if(my_id .eq. 0) print*, ' AC controls: flux=',AC_Flux,' recon=',AC_Recon, &
      ' limiter=',AC_Limiter,' weno_blend=',AC_WenoBlend,' wall_recon=',AC_WallRecon, &
      ' wall_p=',AC_WallP,' mom_diss=',AC_MomDiss,' mom_frac=',AC_MomFrac
   end subroutine ac_check_controls

!===============================================================================
! First-order Rusanov (LLF) flux on ONE face between cells L and R.
! Used on the same-class block-interface planes: the interface is closed with
! a robust dissipative flux built DIRECTLY from the two cell states, with no
! reconstruction.  A reconstructed (MUSCL/WENO) interface flux driven by
! ghost data that is only refreshed once per outer step is what excites the
! odd-even instability of the collocated AC scheme; the plain Rusanov flux
! adds the missing dissipation.  RAC(L) -= F*A ; RAC(R) += F*A.
!===============================================================================
   subroutine ac_face_flux_rus1(B, iL,jL,kL, iR,jR,kR, dir, beta)
    use Global_var
    use lows_ac_work
    implicit none
    Type(Block_TYPE),pointer:: B
    integer,intent(in):: iL,jL,kL,iR,jR,kR,dir
    real(PRE_EC),intent(in):: beta
    real(PRE_EC) :: FLX(4), A, nx1,ny1,nz1
    if(dir .eq. 1) then
      A=B%Si(iR,jR,kR); nx1=B%ni1(iR,jR,kR); ny1=B%ni2(iR,jR,kR); nz1=B%ni3(iR,jR,kR)
    else if(dir .eq. 2) then
      A=B%Sj(iR,jR,kR); nx1=B%nj1(iR,jR,kR); ny1=B%nj2(iR,jR,kR); nz1=B%nj3(iR,jR,kR)
    else
      A=B%Sk(iR,jR,kR); nx1=B%nk1(iR,jR,kR); ny1=B%nk2(iR,jR,kR); nz1=B%nk3(iR,jR,kR)
    endif
    call ac_flux_rusanov(XW(1,iL,jL,kL),XW(2,iL,jL,kL),XW(3,iL,jL,kL),XW(4,iL,jL,kL), &
                         XW(1,iR,jR,kR),XW(2,iR,jR,kR),XW(3,iR,jR,kR),XW(4,iR,jR,kR), &
                         nx1,ny1,nz1, beta, FLX)
    RAC(1:4,iL,jL,kL) = RAC(1:4,iL,jL,kL) - FLX(1:4)*A
    RAC(1:4,iR,jR,kR) = RAC(1:4,iR,jR,kR) + FLX(1:4)*A
   end subroutine ac_face_flux_rus1

!===============================================================================
! Same-class block-interface registry (owner = block with the smaller number).
! ===============================================================================
   subroutine ac_build_ifc_registry(nMesh)
    use Global_var
    use const_var
    use lows_ac_work
    implicit none
    integer :: nMesh, mBlock, ksub, nb, mx, my, mz
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc
    if(allocated(FSH)) deallocate(FSH)
    ac_ifc_n = 0
    mx = 1; my = 1; mz = 1
    do mBlock = 1, Mesh(nMesh)%Num_Block
      B => Mesh(nMesh)%Block(mBlock)
      if(B%Block_type .ne. BLOCK_LOWSPEED) cycle
      if(.not. associated(B%bc_msg)) cycle
      do ksub = 1, B%subface
        Bc => B%bc_msg(ksub)
        if(.not. is_interface_bc(Bc%bc)) cycle
        nb = Bc%nb1
        if(nb .le. 0 .or. nb .gt. Mesh(nMesh)%Num_Block) cycle
        if(Block_Type_List(nb) .ne. BLOCK_LOWSPEED) cycle
        if(mBlock .ge. nb) cycle
        ac_ifc_n = ac_ifc_n + 1
        if(ac_ifc_n .le. AC_IFC_MAX) then
          ac_ifc_own(ac_ifc_n) = mBlock
          ac_ifc_ks (ac_ifc_n) = ksub
          ac_ifc_nb (ac_ifc_n) = nb
        endif
        mx = max(mx, B%nx); my = max(my, B%ny); mz = max(mz, B%nz)
      enddo
    enddo
    if(ac_ifc_n .gt. AC_IFC_MAX) then
      if(my_id .eq. 0) write(*,*) ' AC: >', AC_IFC_MAX, ' interfaces: sharing cut'
      ac_ifc_n = AC_IFC_MAX
    endif
    if(ac_ifc_n .gt. 0) then
      allocate( FSH(4, mx, my, mz, ac_ifc_n) )
      FSH = 0.d0
    endif
    ac_ifc_mesh = nMesh
    if(my_id .eq. 0 .and. ac_ifc_n .gt. 0) &
      write(*,*) ' AC: shared same-class interfaces =', ac_ifc_n
   end subroutine ac_build_ifc_registry

!===============================================================================
! (mBlock, ksub) -> id, is_owner.  id = 0 -> fall back to the two-sided flux.
! ===============================================================================
   subroutine ac_ifc_find(nMesh, mBlock, ksub, id, is_owner)
    use Global_var
    use const_var
    use lows_ac_work
    implicit none
    integer, intent(in) :: nMesh, mBlock, ksub
    integer, intent(out) :: id
    logical, intent(out) :: is_owner
    integer :: k
    Type (Block_TYPE),pointer:: Bo
    TYPE (BC_MSG_TYPE),pointer:: Bc, Bco
    id = 0; is_owner = .false.
    do k = 1, ac_ifc_n
      if(ac_ifc_own(k) .eq. mBlock .and. ac_ifc_ks(k) .eq. ksub) then
        id = k; is_owner = .true.; return
      endif
    enddo
    if(mBlock .lt. 1 .or. mBlock .gt. Mesh(nMesh)%Num_Block) return
    if(.not. associated(Mesh(nMesh)%Block(mBlock)%bc_msg)) return
    Bc => Mesh(nMesh)%Block(mBlock)%bc_msg(ksub)
    do k = 1, ac_ifc_n
      if(ac_ifc_nb(k) .ne. mBlock) cycle
      Bo => Mesh(nMesh)%Block(ac_ifc_own(k))
      if(.not. associated(Bo%bc_msg)) cycle
      Bco => Bo%bc_msg(ac_ifc_ks(k))
      if(Bco%face1 .ne. Bc%face) cycle
      if(Bco%ib1 .ne. Bc%ib .or. Bco%ie1 .ne. Bc%ie) cycle
      if(Bco%jb1 .ne. Bc%jb .or. Bco%je1 .ne. Bc%je) cycle
      if(Bco%kb1 .ne. Bc%kb .or. Bco%ke1 .ne. Bc%ke) cycle
      id = k; is_owner = .false.; return
    enddo
   end subroutine ac_ifc_find

!===============================================================================
! Owner cell opposite the given ghost point, same index convention as the
! buffer exchange (Umessage_send_mpi) -> rotated interfaces stay consistent.
! ===============================================================================
   subroutine ac_ifc_src_cell(nMesh, id, ig,jg,kg, i1,j1,k1)
    use Global_var
    use const_var
    use lows_ac_work
    implicit none
    integer, intent(in) :: nMesh, id, ig,jg,kg
    integer, intent(out) :: i1,j1,k1
    integer :: k, d, kb(3), ke(3), kb1(3), ke1(3), ks(3), ka(3), L(3), P(3), g(3)
    TYPE (BC_MSG_TYPE),pointer:: Bc
    Bc => Mesh(nMesh)%Block(ac_ifc_own(id))%bc_msg(ac_ifc_ks(id))
    kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
    d = mod(Bc%face-1,3)+1
    if(Bc%face .gt. 3) kb(d) = kb(d) - LAP
    ke(d) = kb(d) + LAP - 1
    kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
    d = mod(Bc%face1-1,3)+1
    if(Bc%face1 .le. 3) kb1(d) = kb1(d) - LAP
    ke1(d) = kb1(d) + LAP - 1
    L(1)=abs(Bc%L1); P(1)=sign(1,Bc%L1)
    L(2)=abs(Bc%L2); P(2)=sign(1,Bc%L2)
    L(3)=abs(Bc%L3); P(3)=sign(1,Bc%L3)
    do k = 1, 3
      if(P(k) .gt. 0) then
        ks(k)=kb(k)
      else
        ks(k)=ke(k)
      endif
    enddo
    g(1)=ig; g(2)=jg; g(3)=kg
    do k = 1, 3
      ka(k) = g(k) - kb1(k)
    enddo
    i1 = ks(1) + ka(L(1))*P(1)
    j1 = ks(2) + ka(L(2))*P(2)
    k1 = ks(3) + ka(L(3))*P(3)
   end subroutine ac_ifc_src_cell

!===============================================================================
! Owner: single valued interface flux -> FSH(owner cell); RAC(owner cell) -= F*A
! ===============================================================================
   subroutine ac_ifc_flux_own(B, id, iO,jO,kO, iG,jG,kG, dir, beta)
    use Global_var
    use lows_ac_work
    implicit none
    Type (Block_TYPE),pointer:: B
    integer,intent(in):: id, iO,jO,kO, iG,jG,kG, dir
    real(PRE_EC),intent(in):: beta
    real(PRE_EC) :: FLX(4), A, nx1,ny1,nz1
    integer :: iF,jF,kF
    iF=max(iO,iG); jF=max(jO,jG); kF=max(kO,kG)          ! face plane index
    if(dir .eq. 1) then
      A=B%Si(iF,jF,kF); nx1=B%ni1(iF,jF,kF); ny1=B%ni2(iF,jF,kF); nz1=B%ni3(iF,jF,kF)
    else if(dir .eq. 2) then
      A=B%Sj(iF,jF,kF); nx1=B%nj1(iF,jF,kF); ny1=B%nj2(iF,jF,kF); nz1=B%nj3(iF,jF,kF)
    else
      A=B%Sk(iF,jF,kF); nx1=B%nk1(iF,jF,kF); ny1=B%nk2(iF,jF,kF); nz1=B%nk3(iF,jF,kF)
    endif
    call ac_flux_rusanov(XW(1,iO,jO,kO),XW(2,iO,jO,kO),XW(3,iO,jO,kO),XW(4,iO,jO,kO), &
                         XW(1,iG,jG,kG),XW(2,iG,jG,kG),XW(3,iG,jG,kG),XW(4,iG,jG,kG), &
                         nx1,ny1,nz1, beta, FLX)
    FSH(1:4,iO,jO,kO,id) = FLX(1:4)
    RAC(1:4,iO,jO,kO) = RAC(1:4,iO,jO,kO) - FLX(1:4)*A
!IFC DIAG
    if(jO .eq. 1 .and. kO .eq. 1) then
      ac_dbg_own = ac_dbg_own + 1
      if(ac_dbg .eq. 1 .and. (ac_dbg_own .le. 3 .or. mod(ac_dbg_own,500) .eq. 0)) then
        write(*,*) ' DBG own id=',id,' cell=',iO,jO,kO,' qvu=',XW(1:4,iO,jO,kO), &
          ' ghost=',XW(1:4,iG,jG,kG),' FLX=',FLX(1:4),' A=',A
      endif
    endif
!IFC DIAG END
   end subroutine ac_ifc_flux_own

!===============================================================================
! Neighbour: NO local flux -- read the single stored owner flux.
! F points from the owner cell into this cell -> RAC(own cell) += F*A
! ===============================================================================
   subroutine ac_ifc_flux_nb(nMesh, B, id, iO,jO,kO, iG,jG,kG, dir)
    use Global_var
    use lows_ac_work
    implicit none
    integer,intent(in):: nMesh, id, iO,jO,kO, iG,jG,kG, dir
    Type (Block_TYPE),pointer:: B
    real(PRE_EC) :: A
    integer :: i1,j1,k1, iF,jF,kF
    call ac_ifc_src_cell(nMesh, id, iG,jG,kG, i1,j1,k1)
    iF=max(iO,iG); jF=max(jO,jG); kF=max(kO,kG)
    if(dir .eq. 1) then
      A=B%Si(iF,jF,kF)
    else if(dir .eq. 2) then
      A=B%Sj(iF,jF,kF)
    else
      A=B%Sk(iF,jF,kF)
    endif
    RAC(1:4,iO,jO,kO) = RAC(1:4,iO,jO,kO) + FSH(1:4,i1,j1,k1,id)*A
!IFC DIAG
    if(jO .eq. 1 .and. kO .eq. 1) then
      ac_dbg_nb = ac_dbg_nb + 1
      if(ac_dbg .eq. 1 .and. (ac_dbg_nb .le. 3 .or. mod(ac_dbg_nb,500) .eq. 0)) then
        write(*,*) ' DBG nb  id=',id,' cell=',iO,jO,kO,' qvu=',XW(1:4,iO,jO,kO), &
          ' ghost=',XW(1:4,iG,jG,kG),' map=',i1,j1,k1, &
          ' FSHread=',FSH(1:4,i1,j1,k1,id),' A=',A
      endif
    endif
!IFC DIAG END
   end subroutine ac_ifc_flux_nb

!===============================================================================
! One interface face point: owner computes+stores, neighbour reads, and only if
! the pair could not be resolved does it fall back to the two-sided flux.
! ===============================================================================
   subroutine ac_ifc_face(nMesh, B, id, is_owner, iO,jO,kO, iG,jG,kG, dir, beta)
    use Global_var
    use precision_EC
    implicit none
    integer,intent(in):: nMesh, id, iO,jO,kO, iG,jG,kG, dir
    logical,intent(in):: is_owner
    real(PRE_EC),intent(in):: beta
    Type (Block_TYPE),pointer:: B
    interface
      subroutine ac_ifc_flux_own(B, id, iO,jO,kO, iG,jG,kG, dir, beta)
        use Global_var
        use precision_EC
        Type (Block_TYPE),pointer:: B
        integer,intent(in):: id, iO,jO,kO, iG,jG,kG, dir
        real(PRE_EC),intent(in):: beta
      end subroutine ac_ifc_flux_own
      subroutine ac_ifc_flux_nb(nMesh, B, id, iO,jO,kO, iG,jG,kG, dir)
        use Global_var
        use precision_EC
        integer,intent(in):: nMesh, id, iO,jO,kO, iG,jG,kG, dir
        Type (Block_TYPE),pointer:: B
      end subroutine ac_ifc_flux_nb
      subroutine ac_face_flux_rus1(B, iL,jL,kL, iR,jR,kR, dir, beta)
        use Global_var
        use precision_EC
        Type (Block_TYPE),pointer:: B
        integer,intent(in):: iL,jL,kL,iR,jR,kR,dir
        real(PRE_EC),intent(in):: beta
      end subroutine ac_face_flux_rus1
    end interface
    if(id .le. 0) then
      call ac_face_flux_rus1(B, iO,jO,kO, iG,jG,kG, dir, beta)
    else if(is_owner) then
      call ac_ifc_flux_own(B, id, iO,jO,kO, iG,jG,kG, dir, beta)
    else
      call ac_ifc_flux_nb(nMesh, B, id, iO,jO,kO, iG,jG,kG, dir)
    endif
   end subroutine ac_ifc_face

!===============================================================================
! All same-class interface faces of one block.  Blocks must be visited in
! ascending block order inside a pseudo-step (the AC drivers do), so the owner
! value that the other side reads is always from the current pseudo-step.
! ===============================================================================
   subroutine ac_interface_fluxes(nMesh, mBlock, beta)
    use Global_var
    use const_var
    use lows_ac_work
    implicit none
    integer,intent(in):: nMesh, mBlock
    real(PRE_EC),intent(in):: beta
    integer :: ksub, id, face_s, ib,ie,jb,je,kb,ke, i,j,k
    logical :: is_owner
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc
    interface
      subroutine ac_build_ifc_registry(nMesh)
        integer :: nMesh
      end subroutine ac_build_ifc_registry
      subroutine ac_ifc_find(nMesh, mBlock, ksub, id, is_owner)
        integer, intent(in) :: nMesh, mBlock, ksub
        integer, intent(out) :: id
        logical, intent(out) :: is_owner
      end subroutine ac_ifc_find
      subroutine ac_ifc_face(nMesh, B, id, is_owner, iO,jO,kO, iG,jG,kG, dir, beta)
        use Global_var
        use precision_EC
        integer,intent(in):: nMesh, id, iO,jO,kO, iG,jG,kG, dir
        logical,intent(in):: is_owner
        real(PRE_EC),intent(in):: beta
        Type (Block_TYPE),pointer:: B
      end subroutine ac_ifc_face
    end interface
    if(ac_ifc_mesh .ne. nMesh) call ac_build_ifc_registry(nMesh)
    B => Mesh(nMesh)%Block(mBlock)
    if(.not. associated(B%bc_msg)) return
    do ksub = 1, B%subface
      Bc => B%bc_msg(ksub)
      if(.not. is_interface_bc(Bc%bc)) cycle
      if(Bc%nb1 .le. 0 .or. Bc%nb1 .gt. Mesh(nMesh)%Num_Block) cycle
      if(Block_Type_List(Bc%nb1) .ne. BLOCK_LOWSPEED) cycle
      face_s = Bc%face
      ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke
      call ac_ifc_find(nMesh, mBlock, ksub, id, is_owner)
!IFC DIAG
      if(ac_dbg .eq. 1 .and. ac_dbg_find .lt. 14) then
        ac_dbg_find = ac_dbg_find + 1
        write(*,*) ' DBG ifc blk=',mBlock,' ksub=',ksub,' bc=',Bc%bc, &
          ' nb1=',Bc%nb1,' face=',Bc%face,' face_s=',face_s,' rng=',ib,ie,jb,je,kb,ke, &
          ' own=',is_owner,' id=',id
      endif
!IFC DIAG END
      select case(face_s)
      case(1)                                    ! i- : cell ib   , ghost ib-1
        do k=kb,ke-1; do j=jb,je-1
          call ac_ifc_face(nMesh, B, id, is_owner, ib,j,k, ib-1,j,k, 1, beta)
        enddo; enddo
      case(4)                                    ! i+ : cell ie-1 , ghost ie
        do k=kb,ke-1; do j=jb,je-1
          call ac_ifc_face(nMesh, B, id, is_owner, ie-1,j,k, ie,j,k, 1, beta)
        enddo; enddo
      case(2)                                    ! j- : cell jb   , ghost jb-1
        do k=kb,ke-1; do i=ib,ie-1
          call ac_ifc_face(nMesh, B, id, is_owner, i,jb,k, i,jb-1,k, 2, beta)
        enddo; enddo
      case(5)                                    ! j+ : cell je-1 , ghost je
        do k=kb,ke-1; do i=ib,ie-1
          call ac_ifc_face(nMesh, B, id, is_owner, i,je-1,k, i,je,k, 2, beta)
        enddo; enddo
      case(3)                                    ! k- : cell kb   , ghost kb-1
        do j=jb,je-1; do i=ib,ie-1
          call ac_ifc_face(nMesh, B, id, is_owner, i,j,kb, i,j,kb-1, 3, beta)
        enddo; enddo
      case(6)                                    ! k+ : cell ke-1 , ghost ke
        do j=jb,je-1; do i=ib,ie-1
          call ac_ifc_face(nMesh, B, id, is_owner, i,j,ke-1, i,j,ke, 3, beta)
        enddo; enddo
      end select
    enddo
   end subroutine ac_interface_fluxes

!===============================================================================
! Internal-face inviscid fluxes -> RAC.  RAC must be zeroed before use.
! Plane between cells L=(f-1) and R=f :  L gets -F*A, R gets +F*A
!===============================================================================
  subroutine ac_internal_flux(nMesh, mBlock, beta)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, ii, i, j, k
   real(PRE_EC),intent(in) :: beta
   Type (Block_TYPE),pointer:: B
   interface
     subroutine ac_face_flux(B, iL,jL,kL, iR,jR,kR, dir, beta)
       use Global_var
       use precision_EC
       Type (Block_TYPE),pointer:: B
       integer,intent(in):: iL,jL,kL,iR,jR,kR,dir
       real(PRE_EC),intent(in):: beta
     end subroutine ac_face_flux
     subroutine ac_interface_fluxes(nMesh, mBlock, beta)
       use precision_EC
       integer,intent(in):: nMesh, mBlock
       real(PRE_EC),intent(in):: beta
     end subroutine ac_interface_fluxes
   end interface
   B => Mesh(nMesh)%Block(mBlock)

!  Same-class (LOWSPEED<->LOWSPEED) interface planes: owner computes / both
!  sides share the single flux value (see ac_interface_fluxes below).
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do ii = 2, B%nx-1                  ! i-faces between interior cells
     call ac_face_flux(B, ii-1,j,k, ii,j,k, 1, beta)
   enddo; enddo; enddo

   do k = 1, B%nz-1
   do i = 1, B%nx-1
   do ii = 2, B%ny-1                  ! j-faces between interior cells
     call ac_face_flux(B, i,ii-1,k, i,ii,k, 2, beta)
   enddo; enddo; enddo

   do j = 1, B%ny-1
   do i = 1, B%nx-1
   do ii = 2, B%nz-1                  ! k-faces between interior cells
     call ac_face_flux(B, i,j,ii-1, i,j,ii, 3, beta)
   enddo; enddo; enddo
!  --- same-class block-interface planes: single-valued (owner/shared) flux ---
   call ac_interface_fluxes(nMesh, mBlock, beta)
  end subroutine ac_internal_flux

!===============================================================================
! van-Leer-type slope limiter
!===============================================================================
   real(PRE_EC) function ac_lim(a, b)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: a, b
    real(PRE_EC) :: eps
    eps = 1.d-30
    if(a*b .le. 0.d0) then
      ac_lim = 0.d0
    else
      ac_lim = 2.d0*a*b/(a+b+eps)
    endif
   end function ac_lim

!===============================================================================
! Physical-boundary face fluxes -> RAC (explicit BC flux, no ghost reliance)
! Also zeroes RAC at entry (internal-face routine only accumulates).
!   wall / symmetry : F = (0, q*nx, q*ny, q*nz)            (u_n = 0)
!   velocity/mass inlet : F = F(state=uin, q=cell q)
!   pressure inlet      : F = F(state=cell, q=LS_P_in/rho)
!   outlet / farfield   : F = F(state=cell, q=LS_P_out/rho)
!===============================================================================
!===============================================================================
! Wall/symmetry face pressure for AC_WallP=1.
! The ghost cell centre is the mirror of the first interior cell, so a linear
! reconstruction of the interior field through the face gives the face value
!    q_wall = 1.5 q_1 - 0.5 q_2
! with q_2 the next cell in the direction of the boundary face.  The classic
! "q_wall = q_1" (zero normal pressure gradient) drops the term (dr/2)*dq/dn; on
! a curved wall the normal momentum balance gives dq/dn = u_t^2/R, which is the
! leading part of the inviscid-cylinder suction-peak Cp error.  On a flat wall /
! in a boundary layer dq/dn = 0 and the two forms are identical, so the change is
! neutral exactly where the previous form is already right.
! Cells without a wall/symmetry face keep their own pressure (q2 = q1).
!===============================================================================
   subroutine ac_wall_pressure(nMesh, mBlock)
    use Global_var
    use lows_ac_work
    implicit none
    integer :: nMesh, mBlock, i,j,k, nx,ny,nz
    real(PRE_EC) :: q1, q2
    Type (Block_TYPE),pointer:: B
    if(.not. allocated(QWF)) return       ! AC work arrays not set up (yet)
    B => Mesh(nMesh)%Block(mBlock)
    nx = B%nx; ny = B%ny; nz = B%nz
    do k = 1, nz-1
    do j = 1, ny-1
    do i = 1, nx-1
      q1 = XW(1,i,j,k)
      q2 = q1
      if(AC_WallP .eq. 1) then
        if(ac_wallface(i,j,k,1) .eq. 1) then
          q2 = XW(1,min(i+1,nx-1),j,k)
        else if(ac_wallface(i,j,k,1) .eq. 2) then
          q2 = XW(1,max(i-1,1),j,k)
        else if(ac_wallface(i,j,k,2) .eq. 1) then
          q2 = XW(1,i,min(j+1,ny-1),k)
        else if(ac_wallface(i,j,k,2) .eq. 2) then
          q2 = XW(1,i,max(j-1,1),k)
        else if(ac_wallface(i,j,k,3) .eq. 1) then
          q2 = XW(1,i,j,min(k+1,nz-1))
        else if(ac_wallface(i,j,k,3) .eq. 2) then
          q2 = XW(1,i,j,max(k-1,1))
        endif
      endif
      QWF(i,j,k) = 1.5d0*q1 - 0.5d0*q2
    enddo; enddo; enddo
   end subroutine ac_wall_pressure

  subroutine ac_boundary_flux(nMesh, mBlock, beta, uin_x, uin_y, uin_z)
   use Global_var
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, ksub, i,j,k
   real(PRE_EC) :: beta, uin_x, uin_y, uin_z
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer :: face_s, ib,ie,jb,je,kb,ke, ic,jc,kc
   real(PRE_EC) :: F(4), A, nx1,ny1,nz1
   integer :: qb

   B => Mesh(nMesh)%Block(mBlock)
   RAC = 0.d0

   do ksub = 1, B%subface
     Bc => B%bc_msg(ksub)
     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke

!    --- classify BC and build face state ---
!    No-penetration interfaces (LOWSPEED<->SOLID 13, LOWSPEED<->POROUS 16):
!    these subfaces carry no physical BC but DO need a flux on the adjacent
!    interior cell.  Previously they were `cycle`d here and are also outside
!    the ac_internal_flux interior-face loops, so NO inviscid (pressure /
!    convective / beta*un) flux was applied while the VISCOUS flux was
!    (ac_viscous_res loops all i/j/k planes).  The interface cell's normal
!    momentum was then balanced only by a spurious viscous drag, producing a
!    spurious interface pressure spike and normal velocity (BJ case: v ~ -5e-3
!    on the interface row, p discontinuity ~ 4e-3) and a wrong dp/dx.
!    On a flat, aligned interface the correct kinematic condition is
!    no-penetration with the pressure force transmitted -> use the same
!    pressure-only flux as a wall; the tangential shear still comes from the
!    viscous term using the neighbour ghost.  Interfaces that carry mass flux
!    (e.g. 19 blowing) keep the old (no-flux) path for now.
     qb = -1.d0
     if(Bc%bc .eq. BC_Interface_LowSolid .or. &
        Bc%bc .eq. BC_Interface_LowPorous) qb = 0.d0
     if(associated(B%bc_msg2)) then
       if(B%bc_msg2(ksub)%bc .eq. BC_Interface_LowSolid .or. &
          B%bc_msg2(ksub)%bc .eq. BC_Interface_LowPorous) qb = 0.d0
     endif
     if(qb .lt. 0.d0) then
!      other interfaces (12/19..., negative inner codes) keep old no-flux path
       if(is_interface_bc(Bc%bc)) cycle
       if(associated(B%bc_msg2)) then
         if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
       endif
       if(Bc%bc .eq. BC_Wall .or. Bc%bc .eq. BC_Symmetry) then
         qb = 0.d0                      ! pressure-only flux marker
       else if(Bc%bc .eq. BC_Outflow .or. Bc%bc .eq. BC_LS_Outlet .or. &
               Bc%bc .eq. BC_Farfield) then
         qb = -2.d0                     ! prescribed outlet pressure
       else if(Bc%bc .eq. BC_Inflow .or. Bc%bc .eq. BC_LS_Inlet) then
         if(LS_Inlet_Type .eq. 3) then
           qb = -3.d0                   ! pressure inlet
         else
           qb = -4.d0                   ! velocity inlet
         endif
       endif
     endif

     select case(face_s)
     case(1)                          ! i-  (right cell ib)
       do k=kb,ke-1; do j=jb,je-1
         ic=ib; jc=j; kc=k
         A=B%Si(ic,jc,kc); nx1=B%ni1(ic,jc,kc); ny1=B%ni2(ic,jc,kc); nz1=B%ni3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(4)                          ! i+  (left cell ie-1)
       do k=kb,ke-1; do j=jb,je-1
         ic=ie-1; jc=j; kc=k
         A=B%Si(ie,jc,kc); nx1=B%ni1(ie,jc,kc); ny1=B%ni2(ie,jc,kc); nz1=B%ni3(ie,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) - F(1:4)*A
       enddo; enddo
     case(2)                          ! j-  (right cell jb)
       do k=kb,ke-1; do i=ib,ie-1
         ic=i; jc=jb; kc=k
         A=B%Sj(ic,jc,kc); nx1=B%nj1(ic,jc,kc); ny1=B%nj2(ic,jc,kc); nz1=B%nj3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(5)                          ! j+  (left cell je-1)
       do k=kb,ke-1; do i=ib,ie-1
         ic=i; jc=je-1; kc=k
         A=B%Sj(ic,je,kc); nx1=B%nj1(ic,je,kc); ny1=B%nj2(ic,je,kc); nz1=B%nj3(ic,je,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) - F(1:4)*A
       enddo; enddo
     case(3)                          ! k-  (right cell kb)
       do j=jb,je-1; do i=ib,ie-1
         ic=i; jc=j; kc=kb
         A=B%Sk(ic,jc,kc); nx1=B%nk1(ic,jc,kc); ny1=B%nk2(ic,jc,kc); nz1=B%nk3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(6)                          ! k+  (left cell ke-1)
       do j=jb,je-1; do i=ib,ie-1
         ic=i; jc=j; kc=ke-1
         A=B%Sk(ic,jc,ke); nx1=B%nk1(ic,jc,ke); ny1=B%nk2(ic,jc,ke); nz1=B%nk3(ic,jc,ke)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, QWF(ic,jc,kc), F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) - F(1:4)*A
       enddo; enddo
     end select
   enddo
  end subroutine ac_boundary_flux

!===============================================================================
! Face flux for a physical boundary, mode:
!   0 = wall/symmetry  (pressure-only: F = (0, q n));  u_n set to 0
!  -2 = outlet pressure, -3 = pressure inlet, -4 = velocity inlet (default)
!===============================================================================
  subroutine ac_bc_flux(mode, qc,uc,vc,wc, uin_x,uin_y,uin_z, nx1,ny1,nz1, &
                        beta, pin, pout, qw, F)
   use precision_EC
   use Global_var
   implicit none
   integer,intent(in) :: mode
   real(PRE_EC),intent(in) :: qc,uc,vc,wc, uin_x,uin_y,uin_z
   real(PRE_EC),intent(in) :: nx1,ny1,nz1, beta, pin, pout, qw
   real(PRE_EC),intent(out):: F(4)
   real(PRE_EC) :: qb, ub,vb,wb, un

   if(mode .eq. 0) then                ! wall / symmetry
!    qw is the wall-face pressure (AC_WallP): 1.5q_1-0.5q_2, i.e. the linear
!    reconstruction of the interior field, which retains the wall-normal
!    variation dq/dn that "q_wall = q_1" drops (see ac_wall_pressure).
     qb = qw
     F(1) = 0.d0
     F(2) = qb*nx1; F(3) = qb*ny1; F(4) = qb*nz1
     return
   endif
   if(mode .eq. -2) then               ! outlet: extrapolated velocity + p_out
     qb = pout/max(LS_rho,1.d-20); ub = uc; vb = vc; wb = wc
   else if(mode .eq. -3) then          ! pressure inlet
     qb = pin/max(LS_rho,1.d-20); ub = uc; vb = vc; wb = wc
   else                                ! velocity / mass-flow inlet
     qb = qc
     ub = uin_x; vb = uin_y; wb = uin_z
   endif
   un = ub*nx1 + vb*ny1 + wb*nz1
   F(1) = beta*un
   F(2) = ub*un + qb*nx1
   F(3) = vb*un + qb*ny1
   F(4) = wb*un + qb*nz1
  end subroutine ac_bc_flux

!===============================================================================
! Pseudo time step per cell:
!   dt_con  = AC_CFL*Vol/SIG, SIG = sum over 6 faces of (|un|+c)*A
!   dt_vis  = AC_CFLv*Vol/sum(nu*A/d)
!===============================================================================
  subroutine ac_compute_dt(nMesh, mBlock, beta)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k
   real(PRE_EC) :: beta, un, c, sv, nu
   Type (Block_TYPE),pointer:: B
   interface
     real(PRE_EC) function dist_ac(B,i,j,k,dir)
       use precision_EC
       use Global_var
       Type (Block_TYPE),pointer:: B
       integer,intent(in):: i,j,k,dir
     end function dist_ac
   end interface
   B => Mesh(nMesh)%Block(mBlock)
   nu = LS_mu/max(LS_rho,1.d-20)

   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     un = XW(2,i,j,k)*B%ni1(i,j,k) + XW(3,i,j,k)*B%ni2(i,j,k) + XW(4,i,j,k)*B%ni3(i,j,k)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = (abs(un)+c)*B%Si(i,j,k)
     un = XW(2,i,j,k)*B%ni1(i+1,j,k) + XW(3,i,j,k)*B%ni2(i+1,j,k) + XW(4,i,j,k)*B%ni3(i+1,j,k)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = SIG(i,j,k) + (abs(un)+c)*B%Si(i+1,j,k)
     un = XW(2,i,j,k)*B%nj1(i,j,k) + XW(3,i,j,k)*B%nj2(i,j,k) + XW(4,i,j,k)*B%nj3(i,j,k)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = SIG(i,j,k) + (abs(un)+c)*B%Sj(i,j,k)
     un = XW(2,i,j,k)*B%nj1(i,j+1,k) + XW(3,i,j,k)*B%nj2(i,j+1,k) + XW(4,i,j,k)*B%nj3(i,j+1,k)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = SIG(i,j,k) + (abs(un)+c)*B%Sj(i,j+1,k)
     un = XW(2,i,j,k)*B%nk1(i,j,k) + XW(3,i,j,k)*B%nk2(i,j,k) + XW(4,i,j,k)*B%nk3(i,j,k)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = SIG(i,j,k) + (abs(un)+c)*B%Sk(i,j,k)
     un = XW(2,i,j,k)*B%nk1(i,j,k+1) + XW(3,i,j,k)*B%nk2(i,j,k+1) + XW(4,i,j,k)*B%nk3(i,j,k+1)
     c  = sqrt(un*un + beta)
     SIG(i,j,k) = SIG(i,j,k) + (abs(un)+c)*B%Sk(i,j,k+1)

     DTAC(i,j,k) = AC_CFL * B%Vol(i,j,k) / max(SIG(i,j,k), 1.d-30)

     if(If_viscous .eq. 1 .and. nu .gt. 0.d0) then
       sv = 0.d0
       sv = sv + B%Si(i,j,k)/max(dist_ac(B,i,j,k,1),1.d-30)
       sv = sv + B%Si(i+1,j,k)/max(dist_ac(B,i+1,j,k,1),1.d-30)
       sv = sv + B%Sj(i,j,k)/max(dist_ac(B,i,j,k,2),1.d-30)
       sv = sv + B%Sj(i,j+1,k)/max(dist_ac(B,i,j+1,k,2),1.d-30)
       sv = sv + B%Sk(i,j,k)/max(dist_ac(B,i,j,k,3),1.d-30)
       sv = sv + B%Sk(i,j,k+1)/max(dist_ac(B,i,j,k+1,3),1.d-30)
       DTAC(i,j,k) = min(DTAC(i,j,k), AC_CFLv*B%Vol(i,j,k)/(nu*sv+1.d-30))
     endif
   enddo; enddo; enddo
  end subroutine ac_compute_dt

!  face distance between the two cells across face (d along the coordinate line,
!  projected on the face normal when the grid is skewed; see ac_viscous_res)
   real(PRE_EC) function dist_ac(B, i, j, k, dir)
    use precision_EC
    use Global_var
    Type (Block_TYPE),pointer:: B
    integer,intent(in):: i,j,k,dir
    real(PRE_EC):: dx,dy,dz, dm, dn
    if(dir .eq. 1) then
      dx = B%x(i+1,j,k)-B%x(i-1,j,k); dy = B%y(i+1,j,k)-B%y(i-1,j,k); dz = B%z(i+1,j,k)-B%z(i-1,j,k)
    else if(dir .eq. 2) then
      dx = B%x(i,j+1,k)-B%x(i,j-1,k); dy = B%y(i,j+1,k)-B%y(i,j-1,k); dz = B%z(i,j+1,k)-B%z(i,j-1,k)
    else
      dx = B%x(i,j,k+1)-B%x(i,j,k-1); dy = B%y(i,j,k+1)-B%y(i,j,k-1); dz = B%z(i,j,k+1)-B%z(i,j,k-1)
    endif
    dm = 0.5d0*sqrt(dx*dx+dy*dy+dz*dz)
    if(dir .eq. 1) then
      dn = 0.5d0*abs(dx*B%ni1(i,j,k)+dy*B%ni2(i,j,k)+dz*B%ni3(i,j,k))
    else if(dir .eq. 2) then
      dn = 0.5d0*abs(dx*B%nj1(i,j,k)+dy*B%nj2(i,j,k)+dz*B%nj3(i,j,k))
    else
      dn = 0.5d0*abs(dx*B%nk1(i,j,k)+dy*B%nk2(i,j,k)+dz*B%nk3(i,j,k))
    endif
    if(dn .gt. 0.2d0*dm) then
      dist_ac = dn
    else
      dist_ac = dm
    endif
   end function dist_ac

!===============================================================================
! Dimensionless residual norms (monitoring / convergence):
!   res_q = max_cell |RAC1|/Vol/beta * L/Uref   (~ |div u| normalized)
!   res_m = max_cell |RAC(2:4)|/Vol * L/Uref^2  (~ |momentum residual| norm)
!===============================================================================
  subroutine ac_res_norm(nMesh, mBlock, beta, U2, L, rq, rm)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k
   real(PRE_EC) :: beta, U2, L, rq, rm, vq, vm, rr
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   rq = 0.d0; rm = 0.d0
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     vq = abs(RAC(1,i,j,k))/(max(B%Vol(i,j,k),1.d-30)*beta) * L / sqrt(U2)
     vm = sqrt(RAC(2,i,j,k)**2 + RAC(3,i,j,k)**2 + RAC(4,i,j,k)**2) &
          /max(B%Vol(i,j,k),1.d-30) * L / U2
     rr = abs(RAC(1,i,j,k))/max(B%Vol(i,j,k),1.d-30) * L / U2
     rq = max(rq, vq, 1.d-3*rr)
     rm = max(rm, vm)
   enddo; enddo; enddo
  end subroutine ac_res_norm

!===============================================================================
! Dual-time physical source for momentum (BDF2, coefficients Sfac from main):
!   RAC(2:4) -= (3u - 4*u^n + u^n-1) * Vol * Sfac      (u^n from B%Un)
!===============================================================================
  subroutine ac_dual_res(nMesh, mBlock, Sfac)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k, m
   real(PRE_EC) :: Sfac
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   if(Sfac .le. 0.d0) return
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     RAC(2,i,j,k) = RAC(2,i,j,k) - (3.d0*XW(2,i,j,k) - 4.d0*B%Un(2,i,j,k) + B%Un1(2,i,j,k))*B%Vol(i,j,k)*Sfac
     RAC(3,i,j,k) = RAC(3,i,j,k) - (3.d0*XW(3,i,j,k) - 4.d0*B%Un(3,i,j,k) + B%Un1(3,i,j,k))*B%Vol(i,j,k)*Sfac
     RAC(4,i,j,k) = RAC(4,i,j,k) - (3.d0*XW(4,i,j,k) - 4.d0*B%Un(4,i,j,k) + B%Un1(4,i,j,k))*B%Vol(i,j,k)*Sfac
   enddo; enddo; enddo
  end subroutine ac_dual_res

!===============================================================================
! AC implicit LU-SGS sweeps (forward + backward).  DU4 stores the update;
! XW is refreshed with DU4 at the end.
!===============================================================================
  subroutine ac_lusgs_sweep(nMesh, mBlock, beta, isdual, Sfac1)
   use Global_var
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, isdual, i,j,k, plane, m
   real(PRE_EC) :: beta, Sfac1, w
   real(PRE_EC) :: alfa(4), dui(4), duj(4), duk(4), DF(4)
   real(PRE_EC) :: A, lamf, un, c
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   w = max(AC_w, 0.5d0)

!  ---------------- forward sweep (plane = i+j+k increasing) ----------------
   do plane = 3, (B%nx-1)+(B%ny-1)+(B%nz-1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,i,m,alfa,dui,duj,duk,DF,A,lamf,un,c)
   do k = 1, B%nz-1
   do j = 1, B%ny-1
     i = plane - k - j
     if(i .lt. 1 .or. i .gt. B%nx-1) cycle
     alfa(1) = B%Vol(i,j,k)/max(DTAC(i,j,k),1.d-30) + w*SIG(i,j,k)
     alfa(2) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     alfa(3) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     alfa(4) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     if(isdual .eq. 1) then
       alfa(2) = alfa(2) + Sfac1*B%Vol(i,j,k)
       alfa(3) = alfa(3) + Sfac1*B%Vol(i,j,k)
       alfa(4) = alfa(4) + Sfac1*B%Vol(i,j,k)
     endif

     dui(1:4)=0.d0; duj(1:4)=0.d0; duk(1:4)=0.d0

     if(i .gt. 1) then
       A = B%Si(i,j,k)
       call ac_df_ac(XW(1,i-1,j,k), DU4(1,i-1,j,k), &
                     B%ni1(i,j,k),B%ni2(i,j,k),B%ni3(i,j,k), beta, DF)
       un = XW(2,i-1,j,k)*B%ni1(i,j,k) + XW(3,i-1,j,k)*B%ni2(i,j,k) + XW(4,i-1,j,k)*B%ni3(i,j,k)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       dui(1:4) = 0.5d0*(DF(1:4)*A + w*lamf*DU4(1:4,i-1,j,k))
     endif
     if(j .gt. 1) then
       A = B%Sj(i,j,k)
       call ac_df_ac(XW(1,i,j-1,k), DU4(1,i,j-1,k), &
                     B%nj1(i,j,k),B%nj2(i,j,k),B%nj3(i,j,k), beta, DF)
       un = XW(2,i,j-1,k)*B%nj1(i,j,k) + XW(3,i,j-1,k)*B%nj2(i,j,k) + XW(4,i,j-1,k)*B%nj3(i,j,k)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       duj(1:4) = 0.5d0*(DF(1:4)*A + w*lamf*DU4(1:4,i,j-1,k))
     endif
     if(k .gt. 1) then
       A = B%Sk(i,j,k)
       call ac_df_ac(XW(1,i,j,k-1), DU4(1,i,j,k-1), &
                     B%nk1(i,j,k),B%nk2(i,j,k),B%nk3(i,j,k), beta, DF)
       un = XW(2,i,j,k-1)*B%nk1(i,j,k) + XW(3,i,j,k-1)*B%nk2(i,j,k) + XW(4,i,j,k-1)*B%nk3(i,j,k)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       duk(1:4) = 0.5d0*(DF(1:4)*A + w*lamf*DU4(1:4,i,j,k-1))
     endif
     do m = 1, 4
       DU4(m,i,j,k) = (RAC(m,i,j,k) + dui(m) + duj(m) + duk(m))/alfa(m)
     enddo
   enddo; enddo
!$OMP END PARALLEL DO
   enddo

!  ---------------- backward sweep (plane = i+j+k decreasing) ----------------
   do plane = (B%nx-1)+(B%ny-1)+(B%nz-1), 3, -1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,i,m,alfa,dui,duj,duk,DF,A,lamf,un,c)
   do k = B%nz-1, 1, -1
   do j = B%ny-1, 1, -1
     i = plane - k - j
     if(i .lt. 1 .or. i .gt. B%nx-1) cycle
     alfa(1) = B%Vol(i,j,k)/max(DTAC(i,j,k),1.d-30) + w*SIG(i,j,k)
     alfa(2) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     alfa(3) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     alfa(4) = alfa(1) + dragc(i,j,k)*B%Vol(i,j,k)
     if(isdual .eq. 1) then
       alfa(2) = alfa(2) + Sfac1*B%Vol(i,j,k)
       alfa(3) = alfa(3) + Sfac1*B%Vol(i,j,k)
       alfa(4) = alfa(4) + Sfac1*B%Vol(i,j,k)
     endif
     dui(1:4)=0.d0; duj(1:4)=0.d0; duk(1:4)=0.d0

     if(i .lt. B%nx-1) then
       A = B%Si(i+1,j,k)
       call ac_df_ac(XW(1,i+1,j,k), DU4(1,i+1,j,k), &
                     B%ni1(i+1,j,k),B%ni2(i+1,j,k),B%ni3(i+1,j,k), beta, DF)
       un = XW(2,i+1,j,k)*B%ni1(i+1,j,k) + XW(3,i+1,j,k)*B%ni2(i+1,j,k) + XW(4,i+1,j,k)*B%ni3(i+1,j,k)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       dui(1:4) = -0.5d0*(DF(1:4)*A - w*lamf*DU4(1:4,i+1,j,k))
     endif
     if(j .lt. B%ny-1) then
       A = B%Sj(i,j+1,k)
       call ac_df_ac(XW(1,i,j+1,k), DU4(1,i,j+1,k), &
                     B%nj1(i,j+1,k),B%nj2(i,j+1,k),B%nj3(i,j+1,k), beta, DF)
       un = XW(2,i,j+1,k)*B%nj1(i,j+1,k) + XW(3,i,j+1,k)*B%nj2(i,j+1,k) + XW(4,i,j+1,k)*B%nj3(i,j+1,k)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       duj(1:4) = -0.5d0*(DF(1:4)*A - w*lamf*DU4(1:4,i,j+1,k))
     endif
     if(k .lt. B%nz-1) then
       A = B%Sk(i,j,k+1)
       call ac_df_ac(XW(1,i,j,k+1), DU4(1,i,j,k+1), &
                     B%nk1(i,j,k+1),B%nk2(i,j,k+1),B%nk3(i,j,k+1), beta, DF)
       un = XW(2,i,j,k+1)*B%nk1(i,j,k+1) + XW(3,i,j,k+1)*B%nk2(i,j,k+1) + XW(4,i,j,k+1)*B%nk3(i,j,k+1)
       c = sqrt(un*un+beta); lamf = (abs(un)+c)*A
       duk(1:4) = -0.5d0*(DF(1:4)*A - w*lamf*DU4(1:4,i,j,k+1))
     endif
     do m = 1, 4
       DU4(m,i,j,k) = DU4(m,i,j,k) + (dui(m) + duj(m) + duk(m))/alfa(m)
     enddo
   enddo; enddo
!$OMP END PARALLEL DO
   enddo

!  update XW (interior + also refresh ghost copies from nearest interior so the
!  next boundary fill starts from a consistent state)
   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do i = 1, B%nx-1
     XW(1:4,i,j,k) = XW(1:4,i,j,k) + DU4(1:4,i,j,k)
   enddo; enddo; enddo
  end subroutine ac_lusgs_sweep

!===============================================================================
! Flux increment of the AC system for a state jump DU (implicit operator):
!   DF = F(W+DU).n - F(W).n     (per unit area, normal n)
!===============================================================================
  subroutine ac_df_ac(WV, DWV, anx, any, anz, beta, DF)
   use precision_EC
   implicit none
   real(PRE_EC),intent(in) :: WV(4), DWV(4), anx,any,anz, beta
   real(PRE_EC),intent(out):: DF(4)
   real(PRE_EC) :: q0,u0,v0,w0, dq,du,dv,dw, q1,u1,v1,w1, un0,un1
   q0 = WV(1); u0 = WV(2); v0 = WV(3); w0 = WV(4)
   dq = DWV(1); du = DWV(2); dv = DWV(3); dw = DWV(4)
   q1 = q0+dq; u1 = u0+du; v1 = v0+dv; w1 = w0+dw
   un0 = u0*anx + v0*any + w0*anz
   un1 = u1*anx + v1*any + w1*anz
   DF(1) = beta*(un1-un0)
   DF(2) = (u1*un1 + q1*anx) - (u0*un0 + q0*anx)
   DF(3) = (v1*un1 + q1*any) - (v0*un0 + q0*any)
   DF(4) = (w1*un1 + q1*anz) - (w0*un0 + q0*anz)
  end subroutine ac_df_ac

!===============================================================================
! Viscous momentum terms -> RAC(2:4). Orthogonal-grid Laplacian through face
! gradients H = nu*(phi_R-phi_L)/d.  Accumulation opposite to the inviscid one:
! left cell gets -H*A, right cell gets +H*A  (i.e. + Laplacian).  nu=LS_mu/rho.
!
! d is the CELL-CENTRE distance across the face, d = |0.5*(x_{f+1}-x_{f-1})|.
! At a wall face this equals exactly 2*(distance wall -> first cell centre)
! because the ghost centres are stored as mirrored nodes, which is what makes the
! no-slip shear (u_ghost-u_1)/d = 2u_1/d_wall second-order.  On a skewed
! (non-orthogonal) grid the gradient is taken along the index line, which
! overestimates the face-normal derivative by 1/cos(skew), so d is projected on
! the face normal (B%ni/nj/nk are unit normals) with the index distance kept as
! a fallback for degenerate or strongly skewed faces.  On an orthogonal mesh the
! projection is identical, so those cases are unchanged to round-off.
!===============================================================================
  subroutine ac_viscous_res(nMesh, mBlock)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k, f, m
   real(PRE_EC) :: nu, H, A, d, dxv,dyv,dzv, dproj
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   nu = LS_mu/max(LS_rho,1.d-20)

   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do f = 1, B%nx                    ! i-planes (boundary ones use ghost states)
     A = B%Si(f,j,k)
     dxv = 0.5d0*(B%x(f+1,j,k)-B%x(f-1,j,k))
     dyv = 0.5d0*(B%y(f+1,j,k)-B%y(f-1,j,k))
     dzv = 0.5d0*(B%z(f+1,j,k)-B%z(f-1,j,k))
     d = sqrt(dxv*dxv+dyv*dyv+dzv*dzv)
     dproj = abs(dxv*B%ni1(f,j,k)+dyv*B%ni2(f,j,k)+dzv*B%ni3(f,j,k))
     if(dproj .gt. 0.2d0*d) d = dproj
     d = max(d,1.d-30)
     do m = 2, 4
       H = nu*(XW(m,f,j,k)-XW(m,f-1,j,k))/d
       if(f-1 .ge. 1) RAC(m,f-1,j,k) = RAC(m,f-1,j,k) + H*A
       if(f   .le. B%nx-1) RAC(m,f  ,j,k) = RAC(m,f  ,j,k) - H*A
     enddo
   enddo; enddo; enddo

   do k = 1, B%nz-1
   do i = 1, B%nx-1
   do f = 1, B%ny                    ! j-planes
     A = B%Sj(i,f,k)
     dxv = 0.5d0*(B%x(i,f+1,k)-B%x(i,f-1,k))
     dyv = 0.5d0*(B%y(i,f+1,k)-B%y(i,f-1,k))
     dzv = 0.5d0*(B%z(i,f+1,k)-B%z(i,f-1,k))
     d = sqrt(dxv*dxv+dyv*dyv+dzv*dzv)
     dproj = abs(dxv*B%nj1(i,f,k)+dyv*B%nj2(i,f,k)+dzv*B%nj3(i,f,k))
     if(dproj .gt. 0.2d0*d) d = dproj
     d = max(d,1.d-30)
     do m = 2, 4
       H = nu*(XW(m,i,f,k)-XW(m,i,f-1,k))/d
       if(f-1 .ge. 1) RAC(m,i,f-1,k) = RAC(m,i,f-1,k) + H*A
       if(f   .le. B%ny-1) RAC(m,i,f  ,k) = RAC(m,i,f  ,k) - H*A
     enddo
   enddo; enddo; enddo

   do j = 1, B%ny-1
   do i = 1, B%nx-1
   do f = 1, B%nz                    ! k-planes
     A = B%Sk(i,j,f)
     dxv = 0.5d0*(B%x(i,j,f+1)-B%x(i,j,f-1))
     dyv = 0.5d0*(B%y(i,j,f+1)-B%y(i,j,f-1))
     dzv = 0.5d0*(B%z(i,j,f+1)-B%z(i,j,f-1))
     d = sqrt(dxv*dxv+dyv*dyv+dzv*dzv)
     dproj = abs(dxv*B%nk1(i,j,f)+dyv*B%nk2(i,j,f)+dzv*B%nk3(i,j,f))
     if(dproj .gt. 0.2d0*d) d = dproj
     d = max(d,1.d-30)
     do m = 2, 4
       H = nu*(XW(m,i,j,f)-XW(m,i,j,f-1))/d
       if(f-1 .ge. 1) RAC(m,i,j,f-1) = RAC(m,i,j,f-1) + H*A
       if(f   .le. B%nz-1) RAC(m,i,j,f  ) = RAC(m,i,j,f  ) - H*A
     enddo
   enddo; enddo; enddo
  end subroutine ac_viscous_res

!===============================================================================
! AC ghost-cell filling.  no-slip walls reverse ALL velocity components (point
! mirror about the resting / moving wall); symmetry reverses only the normal
! component; inlet/outlet use Dirichlet/zero-gradient velocity as appropriate.
! Pressure ghost is a copy (not used by the explicit boundary flux).
!===============================================================================
  subroutine ac_fill_ghost(nMesh, mBlock, uin_x, uin_y, uin_z)
   use Global_var
   use const_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, ksub, n, i,j,k
   real(PRE_EC) :: uin_x, uin_y, uin_z
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer :: face_s, ib,ie,jb,je,kb,ke, i1,i2,j1,j2,k1,k2, wsym
   real(PRE_EC) :: nx1,ny1,nz1, ud(3), uw(3), un
   interface
     subroutine ac_set_ghost_cell(B, bc, face_s, i1,j1,k1, i2,j2,k2, &
                                  nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
       use Global_var
       use const_var
       Type (Block_TYPE),pointer:: B
       integer,intent(in) :: bc, face_s, i1,j1,k1, i2,j2,k2
       real(PRE_EC),intent(in) :: nx1,ny1,nz1, uw(3), uin_x,uin_y,uin_z
     end subroutine ac_set_ghost_cell
   end interface

   B => Mesh(nMesh)%Block(mBlock)
!  wall/symmetry face map used by the near-wall reconstruction (AC_WallRecon):
!  rebuilt from the subfaces on every pseudo step, so partial boundary faces are
!  resolved cell by cell.
   if(allocated(ac_wallface)) ac_wallface = 0
   do ksub = 1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle
     if(associated(B%bc_msg2)) then
       if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
     endif
     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke
     wsym = 0
     if(Bc%bc .eq. BC_Wall .or. Bc%bc .eq. BC_Symmetry) wsym = 1

     uw(1)=0.d0; uw(2)=0.d0; uw(3)=0.d0
     if(Bc%bc .eq. BC_Wall .and. face_s .eq. 5) uw(1) = LS_U_lid

     select case(face_s)
     case(1)                                ! i- (plane ib)
       do n=1,LAP
         i1=ib-n; i2=ib+n-1
         do k=kb,ke-1; do j=jb,je-1
           nx1=B%ni1(ib,j,k); ny1=B%ni2(ib,j,k); nz1=B%ni3(ib,j,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i1,j,k, i2,j,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(ib,j,k,1) = 1
         enddo; enddo
       enddo
     case(4)                                ! i+ (plane ie)
       do n=1,LAP
         i1=ie+n-1; i2=ie-n
         do k=kb,ke-1; do j=jb,je-1
           nx1=B%ni1(ie,j,k); ny1=B%ni2(ie,j,k); nz1=B%ni3(ie,j,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i1,j,k, i2,j,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(ie-1,j,k,1) = 2
         enddo; enddo
       enddo
     case(2)                                ! j- (plane jb)
       do n=1,LAP
         j1=jb-n; j2=jb+n-1
         do k=kb,ke-1; do i=ib,ie-1
           nx1=B%nj1(i,jb,k); ny1=B%nj2(i,jb,k); nz1=B%nj3(i,jb,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j1,k, i,j2,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(i,jb,k,2) = 1
         enddo; enddo
       enddo
     case(5)                                ! j+ (plane je)
       do n=1,LAP
         j1=je+n-1; j2=je-n
         do k=kb,ke-1; do i=ib,ie-1
           nx1=B%nj1(i,je,k); ny1=B%nj2(i,je,k); nz1=B%nj3(i,je,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j1,k, i,j2,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(i,je-1,k,2) = 2
         enddo; enddo
       enddo
     case(3)                                ! k- (plane kb)
       do n=1,LAP
         k1=kb-n; k2=kb+n-1
         do j=jb,je-1; do i=ib,ie-1
           nx1=B%nk1(i,j,kb); ny1=B%nk2(i,j,kb); nz1=B%nk3(i,j,kb)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j,k1, i,j,k2, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(i,j,kb,3) = 1
         enddo; enddo
       enddo
     case(6)                                ! k+ (plane ke)
       do n=1,LAP
         k1=ke+n-1; k2=ke-n
         do j=jb,je-1; do i=ib,ie-1
           nx1=B%nk1(i,j,ke); ny1=B%nk2(i,j,ke); nz1=B%nk3(i,j,ke)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j,k1, i,j,k2, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
           if(wsym .eq. 1) ac_wallface(i,j,ke-1,3) = 2
         enddo; enddo
       enddo
     end select
   enddo
  end subroutine ac_fill_ghost

!===============================================================================
! Fill a single ghost cell (velocity/T/p) for an AC boundary
!===============================================================================
  subroutine ac_set_ghost_cell(B, bc, face_s, i1,j1,k1, i2,j2,k2, &
                               nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
   use Global_var
   use const_var
   implicit none
   Type (Block_TYPE),pointer:: B
   integer,intent(in) :: bc, face_s, i1,j1,k1, i2,j2,k2
   real(PRE_EC),intent(in) :: nx1,ny1,nz1, uw(3), uin_x,uin_y,uin_z
   real(PRE_EC) :: ui,vi,wi, un

   ui = B%U(2,i2,j2,k2); vi = B%U(3,i2,j2,k2); wi = B%U(4,i2,j2,k2)

   select case(bc)
   case(BC_Wall)                              ! no-slip (moving wall supported)
!    Inviscid (Euler) wall: If_viscous=0 means there is no viscous momentum
!    flux, so a no-slip ghost would only corrupt the MUSCL reconstruction of
!    the first interior face.  Use the slip condition (mirror the normal
!    component only), exactly as the compressible solver does for
!    BC_Wall + If_viscous=0 (boundary_Symmetry_or_SlideWall).
     if(If_viscous .eq. 0) then
       un = ui*nx1 + vi*ny1 + wi*nz1
       B%U(2,i1,j1,k1) = ui - 2.d0*un*nx1
       B%U(3,i1,j1,k1) = vi - 2.d0*un*ny1
       B%U(4,i1,j1,k1) = wi - 2.d0*un*nz1
     else
       B%U(2,i1,j1,k1) = 2.d0*uw(1) - ui
       B%U(3,i1,j1,k1) = 2.d0*uw(2) - vi
       B%U(4,i1,j1,k1) = 2.d0*uw(3) - wi
     endif
     if(LS_T_wall > 0.d0) then
       B%U(5,i1,j1,k1) = 2.d0*LS_T_wall - B%U(5,i2,j2,k2)
     else
       B%U(5,i1,j1,k1) = B%U(5,i2,j2,k2)
     endif
   case(BC_Symmetry)                          ! normal-component mirror
     un = ui*nx1 + vi*ny1 + wi*nz1
     B%U(2,i1,j1,k1) = ui - 2.d0*un*nx1
     B%U(3,i1,j1,k1) = vi - 2.d0*un*ny1
     B%U(4,i1,j1,k1) = wi - 2.d0*un*nz1
     B%U(5,i1,j1,k1) = B%U(5,i2,j2,k2)
   case(BC_Inflow, BC_LS_Inlet)
     if(LS_Inlet_Type .eq. 3) then
       B%U(2,i1,j1,k1) = ui; B%U(3,i1,j1,k1) = vi; B%U(4,i1,j1,k1) = wi
     else
       B%U(2,i1,j1,k1) = uin_x; B%U(3,i1,j1,k1) = uin_y; B%U(4,i1,j1,k1) = uin_z
     endif
     if(B%Block_type .eq. BLOCK_POROUS .and. Porous_T_in .gt. 0.d0) then
       B%U(5,i1,j1,k1) = Porous_T_in
     else
       B%U(5,i1,j1,k1) = LS_T_ref
     endif
   case default                                ! outlet / farfield: zero gradient
     B%U(2,i1,j1,k1) = ui; B%U(3,i1,j1,k1) = vi; B%U(4,i1,j1,k1) = wi
     B%U(5,i1,j1,k1) = B%U(5,i2,j2,k2)
   end select
   B%p(i1,j1,k1)  = B%p(i2,j2,k2)
   B%U(1,i1,j1,k1) = LS_rho
  end subroutine ac_set_ghost_cell

!===============================================================================
! AC low-speed temperature (passive scalar) solver -- runs after the AC flow
! has converged so the fluid temperature can couple with solids / low-speed
! interfaces (conjugate cases, LS_k > 0).
!
! The face mass fluxes are reconstructed from the converged AC velocity field
! and the steady convection-diffusion equation for T (=B%U(5)) is relaxed with
! the SAME finite-volume discretisation as the SIMPLE path
! (src/sub_lowspeed.f90::lowspeed_energy):
!   rho*cp u.grad(T) = div(k_f grad T),   k_f = LS_k
! Ghost cells: physical faces refreshed by ac_fill_ghost (isothermal walls use
! LS_T_wall, inlets LS_T_ref, outlets/extrapolation); interface (cross-block)
! faces are skipped by ac_fill_ghost so any manual interface ghost set by the
! caller (e.g. staggered conjugate coupling) is preserved.
!===============================================================================
  subroutine ac_lowspeed_energy(nMesh, mBlock)
   use Global_var
   use const_var
   use lowspeed_work
   implicit none
   integer :: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer :: nx,ny,nz, i,j,k, iter
   real(PRE_EC) :: rho, uin_x, uin_y, uin_z, dT, Told_v
   real(PRE_EC), allocatable :: Told(:,:,:)
   interface
     subroutine lowspeed_energy(nMesh, mBlock)
       use precision_EC
       implicit none
       integer, intent(in) :: nMesh, mBlock
     end subroutine lowspeed_energy
   end interface

   B => Mesh(nMesh)%Block(mBlock)
   nx = B%nx; ny = B%ny; nz = B%nz
   rho = max(LS_rho, 1.d-20)

!  face-flux / deferred-correction work arrays (same shape as the SIMPLE path)
   if(.not. allocated(Fi) .or. nxw /= nx .or. nyw /= ny .or. nzw /= nz) then
     if(allocated(Fi)) deallocate(Fi,Fj,Fk,apu,apv,apw,su_nb,sv_nb,sw_nb,du,dv,dw,pp,conv_src)
     allocate(Fi(nx,ny,nz), Fj(nx,ny,nz), Fk(nx,ny,nz))
     allocate(apu(nx,ny,nz), apv(nx,ny,nz), apw(nx,ny,nz))
     allocate(su_nb(nx,ny,nz), sv_nb(nx,ny,nz), sw_nb(nx,ny,nz))
     allocate(du(nx,ny,nz), dv(nx,ny,nz), dw(nx,ny,nz))
     allocate(pp(0:nx,0:ny,0:nz))
     allocate(conv_src(nx,ny,nz))
     nxw=nx; nyw=ny; nzw=nz
   endif
   Fi=0.d0; Fj=0.d0; Fk=0.d0; conv_src=0.d0

!  face mass fluxes (kg/s) from the converged AC velocity field; Fi/Fj/Fk are
!  the flux on face plane i/j/k between cell (plane-1) and (plane), so include
!  the boundary planes 1 and nx (ny, nz) with the ghost-cell velocity -- for
!  walls the mirror ghost gives zero net flux, matching the SIMPLE convention.
   call lowspeed_inlet_velocity(B, uin_x, uin_y, uin_z)
   do k = 1, nz-1
   do j = 1, ny-1
   do i = 1, nx
     Fi(i,j,k) = rho*B%Si(i,j,k)*0.5d0*(B%U(2,i-1,j,k)+B%U(2,i,j,k))
   enddo; enddo; enddo
   do k = 1, nz-1
   do i = 1, nx-1
   do j = 1, ny
     Fj(i,j,k) = rho*B%Sj(i,j,k)*0.5d0*(B%U(3,i,j-1,k)+B%U(3,i,j,k))
   enddo; enddo; enddo
   do j = 1, ny-1
   do i = 1, nx-1
   do k = 1, nz
     Fk(i,j,k) = rho*B%Sk(i,j,k)*0.5d0*(B%U(4,i,j,k-1)+B%U(4,i,j,k))
   enddo; enddo; enddo

!  relax the steady temperature equation (lowspeed_energy, same discretisation
!  as SIMPLE) until the change over a window is small
   allocate(Told(nx,ny,nz)); Told = 0.d0
   dT = 0.d0
   do iter = 1, max(LS_Max_Iter, 2000)
     if(mod(iter,50) .eq. 1) then
       do k = 1, nz-1; do j = 1, ny-1; do i = 1, nx-1
         Told(i,j,k) = B%U(5,i,j,k)
       enddo; enddo; enddo
     endif
     call ac_fill_ghost(nMesh, mBlock, uin_x, uin_y, uin_z)  ! refresh physical T ghosts
     call lowspeed_energy(nMesh, mBlock)
     if(mod(iter,50) .eq. 0) then
       dT = 0.d0
       do k = 1, nz-1; do j = 1, ny-1; do i = 1, nx-1
         dT = max(dT, abs(B%U(5,i,j,k)-Told(i,j,k)))
       enddo; enddo; enddo
       if(my_id .eq. 0 .and. mod(iter,1000) .eq. 0) then
         print*, '  AC-lowspeed T iter', iter, ' max|dT|/50=', dT, &
                 ' T[', minval(B%U(5,1:nx-1,1:ny-1,1:nz-1)), ',', &
                       maxval(B%U(5,1:nx-1,1:ny-1,1:nz-1)), '] K'
       endif
       if(dT .lt. 1.d-3) exit
     endif
   enddo
   deallocate(Told)
   if(my_id .eq. 0) print*, ' AC-lowspeed T block', mBlock, ': sweeps =', iter, &
       ' final max|dT|(50-window)=', dT
  end subroutine ac_lowspeed_energy


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
     dragc = 0.d0
   endif
   dragc = 0.d0   ! clear residual drag coefficient (porous AC sets it per iteration)

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
      ' nx,ny,nz=',nx,ny,nz,' beta=',beta,' CFL=',AC_CFL,' dual=',isdual

   rq0 = 0.d0; rm0 = 0.d0
   do iter = 1, AC_Max_Iter
     call ac_fill_ghost(nMesh, mBlock, uin_x, uin_y, uin_z)
     call ac_load_state(nMesh, mBlock)
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
   integer :: nMesh, mBlock, ksub
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   logical :: flag
   B => Mesh(nMesh)%Block(mBlock)
   flag = .false.
   do ksub = 1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle
     if(associated(B%bc_msg2)) then
       if(B%bc_msg2(ksub)%bc < 0) cycle
     endif
     if(Bc%bc .eq. BC_Outflow .or. Bc%bc .eq. BC_LS_Outlet .or. &
        Bc%bc .eq. BC_Farfield) flag = .true.
     if((Bc%bc .eq. BC_Inflow .or. Bc%bc .eq. BC_LS_Inlet) .and. &
        LS_Inlet_Type .eq. 3) flag = .true.
   enddo
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
!===============================================================================
  subroutine ac_flux_rusanov(qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta,F)
   use precision_EC
   implicit none
   real(PRE_EC),intent(in) :: qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta
   real(PRE_EC),intent(out):: F(4)
   real(PRE_EC) :: unL,unR,cL,cR,lam
   unL = uL*anx + vL*any + wL*anz
   unR = uR*anx + vR*any + wR*anz
   cL  = sqrt(unL*unL + beta)
   cR  = sqrt(unR*unR + beta)
   lam = max(abs(unL)+cL, abs(unR)+cR, 1.d-12)
   F(1) = 0.5d0*beta*(unL+unR) - 0.5d0*lam*(qR-qL)
   F(2) = 0.5d0*((uL*unL+qL*anx) + (uR*unR+qR*anx)) - 0.5d0*lam*(uR-uL)
   F(3) = 0.5d0*((vL*unL+qL*any) + (vR*unR+qR*any)) - 0.5d0*lam*(vR-vL)
   F(4) = 0.5d0*((wL*unL+qL*anz) + (wR*unR+qR*anz)) - 0.5d0*lam*(wR-wL)
  end subroutine ac_flux_rusanov

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
   real(PRE_EC) :: FLX(4), A
   real(PRE_EC) :: qL,uL,vL,wL, qR,uR,vR,wR
   real(PRE_EC), external :: ac_lim
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)

   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do ii = 2, B%nx-1                  ! i-faces between interior cells
     A = B%Si(ii,j,k)
     qL = XW(1,ii-1,j,k) + 0.5d0*ac_lim(XW(1,ii-1,j,k)-XW(1,ii-2,j,k), XW(1,ii,j,k)-XW(1,ii-1,j,k))
     qR = XW(1,ii  ,j,k) - 0.5d0*ac_lim(XW(1,ii,j,k)-XW(1,ii-1,j,k), XW(1,ii+1,j,k)-XW(1,ii,j,k))
     uL = XW(2,ii-1,j,k) + 0.5d0*ac_lim(XW(2,ii-1,j,k)-XW(2,ii-2,j,k), XW(2,ii,j,k)-XW(2,ii-1,j,k))
     uR = XW(2,ii  ,j,k) - 0.5d0*ac_lim(XW(2,ii,j,k)-XW(2,ii-1,j,k), XW(2,ii+1,j,k)-XW(2,ii,j,k))
     vL = XW(3,ii-1,j,k) + 0.5d0*ac_lim(XW(3,ii-1,j,k)-XW(3,ii-2,j,k), XW(3,ii,j,k)-XW(3,ii-1,j,k))
     vR = XW(3,ii  ,j,k) - 0.5d0*ac_lim(XW(3,ii,j,k)-XW(3,ii-1,j,k), XW(3,ii+1,j,k)-XW(3,ii,j,k))
     wL = XW(4,ii-1,j,k) + 0.5d0*ac_lim(XW(4,ii-1,j,k)-XW(4,ii-2,j,k), XW(4,ii,j,k)-XW(4,ii-1,j,k))
     wR = XW(4,ii  ,j,k) - 0.5d0*ac_lim(XW(4,ii,j,k)-XW(4,ii-1,j,k), XW(4,ii+1,j,k)-XW(4,ii,j,k))
     call ac_flux_rusanov(qL,uL,vL,wL, qR,uR,vR,wR, &
                          B%ni1(ii,j,k), B%ni2(ii,j,k), B%ni3(ii,j,k), beta, FLX)
     RAC(1:4,ii-1,j,k) = RAC(1:4,ii-1,j,k) - FLX(1:4)*A
     RAC(1:4,ii  ,j,k) = RAC(1:4,ii  ,j,k) + FLX(1:4)*A
   enddo; enddo; enddo

   do k = 1, B%nz-1
   do i = 1, B%nx-1
   do ii = 2, B%ny-1                  ! j-faces between interior cells
     A = B%Sj(i,ii,k)
     qL = XW(1,i,ii-1,k) + 0.5d0*ac_lim(XW(1,i,ii-1,k)-XW(1,i,ii-2,k), XW(1,i,ii,k)-XW(1,i,ii-1,k))
     qR = XW(1,i,ii  ,k) - 0.5d0*ac_lim(XW(1,i,ii,k)-XW(1,i,ii-1,k), XW(1,i,ii+1,k)-XW(1,i,ii,k))
     uL = XW(2,i,ii-1,k) + 0.5d0*ac_lim(XW(2,i,ii-1,k)-XW(2,i,ii-2,k), XW(2,i,ii,k)-XW(2,i,ii-1,k))
     uR = XW(2,i,ii  ,k) - 0.5d0*ac_lim(XW(2,i,ii,k)-XW(2,i,ii-1,k), XW(2,i,ii+1,k)-XW(2,i,ii,k))
     vL = XW(3,i,ii-1,k) + 0.5d0*ac_lim(XW(3,i,ii-1,k)-XW(3,i,ii-2,k), XW(3,i,ii,k)-XW(3,i,ii-1,k))
     vR = XW(3,i,ii  ,k) - 0.5d0*ac_lim(XW(3,i,ii,k)-XW(3,i,ii-1,k), XW(3,i,ii+1,k)-XW(3,i,ii,k))
     wL = XW(4,i,ii-1,k) + 0.5d0*ac_lim(XW(4,i,ii-1,k)-XW(4,i,ii-2,k), XW(4,i,ii,k)-XW(4,i,ii-1,k))
     wR = XW(4,i,ii  ,k) - 0.5d0*ac_lim(XW(4,i,ii,k)-XW(4,i,ii-1,k), XW(4,i,ii+1,k)-XW(4,i,ii,k))
     call ac_flux_rusanov(qL,uL,vL,wL, qR,uR,vR,wR, &
                          B%nj1(i,ii,k), B%nj2(i,ii,k), B%nj3(i,ii,k), beta, FLX)
     RAC(1:4,i,ii-1,k) = RAC(1:4,i,ii-1,k) - FLX(1:4)*A
     RAC(1:4,i,ii  ,k) = RAC(1:4,i,ii  ,k) + FLX(1:4)*A
   enddo; enddo; enddo

   do j = 1, B%ny-1
   do i = 1, B%nx-1
   do ii = 2, B%nz-1                  ! k-faces between interior cells
     A = B%Sk(i,j,ii)
     qL = XW(1,i,j,ii-1) + 0.5d0*ac_lim(XW(1,i,j,ii-1)-XW(1,i,j,ii-2), XW(1,i,j,ii)-XW(1,i,j,ii-1))
     qR = XW(1,i,j,ii  ) - 0.5d0*ac_lim(XW(1,i,j,ii)-XW(1,i,j,ii-1), XW(1,i,j,ii+1)-XW(1,i,j,ii))
     uL = XW(2,i,j,ii-1) + 0.5d0*ac_lim(XW(2,i,j,ii-1)-XW(2,i,j,ii-2), XW(2,i,j,ii)-XW(2,i,j,ii-1))
     uR = XW(2,i,j,ii  ) - 0.5d0*ac_lim(XW(2,i,j,ii)-XW(2,i,j,ii-1), XW(2,i,j,ii+1)-XW(2,i,j,ii))
     vL = XW(3,i,j,ii-1) + 0.5d0*ac_lim(XW(3,i,j,ii-1)-XW(3,i,j,ii-2), XW(3,i,j,ii)-XW(3,i,j,ii-1))
     vR = XW(3,i,j,ii  ) - 0.5d0*ac_lim(XW(3,i,j,ii)-XW(3,i,j,ii-1), XW(3,i,j,ii+1)-XW(3,i,j,ii))
     wL = XW(4,i,j,ii-1) + 0.5d0*ac_lim(XW(4,i,j,ii-1)-XW(4,i,j,ii-2), XW(4,i,j,ii)-XW(4,i,j,ii-1))
     wR = XW(4,i,j,ii  ) - 0.5d0*ac_lim(XW(4,i,j,ii)-XW(4,i,j,ii-1), XW(4,i,j,ii+1)-XW(4,i,j,ii))
     call ac_flux_rusanov(qL,uL,vL,wL, qR,uR,vR,wR, &
                          B%nk1(i,j,ii), B%nk2(i,j,ii), B%nk3(i,j,ii), beta, FLX)
     RAC(1:4,i,j,ii-1) = RAC(1:4,i,j,ii-1) - FLX(1:4)*A
     RAC(1:4,i,j,ii  ) = RAC(1:4,i,j,ii  ) + FLX(1:4)*A
   enddo; enddo; enddo
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
     if(is_interface_bc(Bc%bc)) cycle
     if(associated(B%bc_msg2)) then
       if(B%bc_msg2(ksub)%bc < 0) cycle
     endif
     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke

!    --- classify BC and build face state ---
     qb = -1.d0
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

     select case(face_s)
     case(1)                          ! i-  (right cell ib)
       do k=kb,ke-1; do j=jb,je-1
         ic=ib; jc=j; kc=k
         A=B%Si(ic,jc,kc); nx1=B%ni1(ic,jc,kc); ny1=B%ni2(ic,jc,kc); nz1=B%ni3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(4)                          ! i+  (left cell ie-1)
       do k=kb,ke-1; do j=jb,je-1
         ic=ie-1; jc=j; kc=k
         A=B%Si(ie,jc,kc); nx1=B%ni1(ie,jc,kc); ny1=B%ni2(ie,jc,kc); nz1=B%ni3(ie,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) - F(1:4)*A
       enddo; enddo
     case(2)                          ! j-  (right cell jb)
       do k=kb,ke-1; do i=ib,ie-1
         ic=i; jc=jb; kc=k
         A=B%Sj(ic,jc,kc); nx1=B%nj1(ic,jc,kc); ny1=B%nj2(ic,jc,kc); nz1=B%nj3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(5)                          ! j+  (left cell je-1)
       do k=kb,ke-1; do i=ib,ie-1
         ic=i; jc=je-1; kc=k
         A=B%Sj(ic,je,kc); nx1=B%nj1(ic,je,kc); ny1=B%nj2(ic,je,kc); nz1=B%nj3(ic,je,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) - F(1:4)*A
       enddo; enddo
     case(3)                          ! k-  (right cell kb)
       do j=jb,je-1; do i=ib,ie-1
         ic=i; jc=j; kc=kb
         A=B%Sk(ic,jc,kc); nx1=B%nk1(ic,jc,kc); ny1=B%nk2(ic,jc,kc); nz1=B%nk3(ic,jc,kc)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
         RAC(1:4,ic,jc,kc) = RAC(1:4,ic,jc,kc) + F(1:4)*A
       enddo; enddo
     case(6)                          ! k+  (left cell ke-1)
       do j=jb,je-1; do i=ib,ie-1
         ic=i; jc=j; kc=ke-1
         A=B%Sk(ic,jc,ke); nx1=B%nk1(ic,jc,ke); ny1=B%nk2(ic,jc,ke); nz1=B%nk3(ic,jc,ke)
         call ac_bc_flux(qb, XW(1,ic,jc,kc), XW(2,ic,jc,kc), XW(3,ic,jc,kc), XW(4,ic,jc,kc), &
              uin_x,uin_y,uin_z, nx1,ny1,nz1, beta, LS_P_in, LS_P_out, F)
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
                        beta, pin, pout, F)
   use precision_EC
   use Global_var
   implicit none
   integer,intent(in) :: mode
   real(PRE_EC),intent(in) :: qc,uc,vc,wc, uin_x,uin_y,uin_z
   real(PRE_EC),intent(in) :: nx1,ny1,nz1, beta, pin, pout
   real(PRE_EC),intent(out):: F(4)
   real(PRE_EC) :: qb, ub,vb,wb, un

   if(mode .eq. 0) then                ! wall / symmetry
     qb = qc
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

!  face distance between the two cells across face (d along the coordinate line)
   real(PRE_EC) function dist_ac(B, i, j, k, dir)
    use precision_EC
    use Global_var
    Type (Block_TYPE),pointer:: B
    integer,intent(in):: i,j,k,dir
    real(PRE_EC):: dx,dy,dz
    if(dir .eq. 1) then
      dx = B%x(i+1,j,k)-B%x(i-1,j,k); dy = B%y(i+1,j,k)-B%y(i-1,j,k); dz = B%z(i+1,j,k)-B%z(i-1,j,k)
    else if(dir .eq. 2) then
      dx = B%x(i,j+1,k)-B%x(i,j-1,k); dy = B%y(i,j+1,k)-B%y(i,j-1,k); dz = B%z(i,j+1,k)-B%z(i,j-1,k)
    else
      dx = B%x(i,j,k+1)-B%x(i,j,k-1); dy = B%y(i,j,k+1)-B%y(i,j,k-1); dz = B%z(i,j,k+1)-B%z(i,j,k-1)
    endif
    dist_ac = 0.5d0*sqrt(dx*dx+dy*dy+dz*dz)
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
!===============================================================================
  subroutine ac_viscous_res(nMesh, mBlock)
   use Global_var
   use lows_ac_work
   implicit none
   integer :: nMesh, mBlock, i,j,k, f, m
   real(PRE_EC) :: nu, H, A, d
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   nu = LS_mu/max(LS_rho,1.d-20)

   do k = 1, B%nz-1
   do j = 1, B%ny-1
   do f = 1, B%nx                    ! i-planes (boundary ones use ghost states)
     A = B%Si(f,j,k)
     d = 0.5d0*sqrt( (B%x(f+1,j,k)-B%x(f-1,j,k))**2 &
                   + (B%y(f+1,j,k)-B%y(f-1,j,k))**2 &
                   + (B%z(f+1,j,k)-B%z(f-1,j,k))**2 )
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
     d = 0.5d0*sqrt( (B%x(i,f+1,k)-B%x(i,f-1,k))**2 &
                   + (B%y(i,f+1,k)-B%y(i,f-1,k))**2 &
                   + (B%z(i,f+1,k)-B%z(i,f-1,k))**2 )
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
     d = 0.5d0*sqrt( (B%x(i,j,f+1)-B%x(i,j,f-1))**2 &
                   + (B%y(i,j,f+1)-B%y(i,j,f-1))**2 &
                   + (B%z(i,j,f+1)-B%z(i,j,f-1))**2 )
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
   implicit none
   integer :: nMesh, mBlock, ksub, n, i,j,k
   real(PRE_EC) :: uin_x, uin_y, uin_z
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer :: face_s, ib,ie,jb,je,kb,ke, i1,i2,j1,j2,k1,k2
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
   do ksub = 1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle
     if(associated(B%bc_msg2)) then
       if(B%bc_msg2(ksub)%bc < 0) cycle
     endif
     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke

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
         enddo; enddo
       enddo
     case(4)                                ! i+ (plane ie)
       do n=1,LAP
         i1=ie+n-1; i2=ie-n
         do k=kb,ke-1; do j=jb,je-1
           nx1=B%ni1(ie,j,k); ny1=B%ni2(ie,j,k); nz1=B%ni3(ie,j,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i1,j,k, i2,j,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
     case(2)                                ! j- (plane jb)
       do n=1,LAP
         j1=jb-n; j2=jb+n-1
         do k=kb,ke-1; do i=ib,ie-1
           nx1=B%nj1(i,jb,k); ny1=B%nj2(i,jb,k); nz1=B%nj3(i,jb,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j1,k, i,j2,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
     case(5)                                ! j+ (plane je)
       do n=1,LAP
         j1=je+n-1; j2=je-n
         do k=kb,ke-1; do i=ib,ie-1
           nx1=B%nj1(i,je,k); ny1=B%nj2(i,je,k); nz1=B%nj3(i,je,k)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j1,k, i,j2,k, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
     case(3)                                ! k- (plane kb)
       do n=1,LAP
         k1=kb-n; k2=kb+n-1
         do j=jb,je-1; do i=ib,ie-1
           nx1=B%nk1(i,j,kb); ny1=B%nk2(i,j,kb); nz1=B%nk3(i,j,kb)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j,k1, i,j,k2, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
     case(6)                                ! k+ (plane ke)
       do n=1,LAP
         k1=ke+n-1; k2=ke-n
         do j=jb,je-1; do i=ib,ie-1
           nx1=B%nk1(i,j,ke); ny1=B%nk2(i,j,ke); nz1=B%nk3(i,j,ke)
           call ac_set_ghost_cell(B, Bc%bc, face_s, i,j,k1, i,j,k2, &
                nx1,ny1,nz1, uw, uin_x,uin_y,uin_z)
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
     B%U(2,i1,j1,k1) = 2.d0*uw(1) - ui
     B%U(3,i1,j1,k1) = 2.d0*uw(2) - vi
     B%U(4,i1,j1,k1) = 2.d0*uw(3) - wi
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
     B%U(5,i1,j1,k1) = LS_T_ref
   case default                                ! outlet / farfield: zero gradient
     B%U(2,i1,j1,k1) = ui; B%U(3,i1,j1,k1) = vi; B%U(4,i1,j1,k1) = wi
     B%U(5,i1,j1,k1) = B%U(5,i2,j2,k2)
   end select
   B%p(i1,j1,k1)  = B%p(i2,j2,k2)
   B%U(1,i1,j1,k1) = LS_rho
  end subroutine ac_set_ghost_cell


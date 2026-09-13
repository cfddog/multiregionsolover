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
      ' nx,ny,nz=',nx,ny,nz,' beta=',beta,' CFL=',AC_CFL,' dual=',isdual, &
      ' flux=',AC_Flux,' recon=',AC_Recon,' lim=',AC_Limiter

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
       if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
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
! AUSM+ (Liou) flux for the AC system, selected by AC_Flux=3.
! The AC flux separates exactly into a convected part and a pressure part:
!   F.n = un*(beta,u,v,w) + q*(0,nx,ny,nz)
! AUSM+ builds the interface convective speed from the Mach splitting
!   M~ = M+(M_L) + M-(M_R),  u_face = c0*M~
! (Phi taken from the upwind side) and the interface pressure from
!   q_face = P+(M_L)*q_L + P-(M_R)*q_R
! (Liou, J.Comput.Phys. 129 (1996) 364).  Reference speed c0=sqrt(beta+un^2/2)
! follows the AC spectral radius used by the LU-SGS.  For equal states
! M~=M and P++P-=1, so the flux reduces exactly to F(W).
!===============================================================================
   subroutine ac_flux_ausm(qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta,F)
    use precision_EC
    implicit none
    real(PRE_EC),intent(in) :: qL,uL,vL,wL,qR,uR,vR,wR,anx,any,anz,beta
    real(PRE_EC),intent(out):: F(4)
    real(PRE_EC) :: unL,unR,c0,ML,MR,MpL,MmR,Mt,uf,ppL,pmR,pf,b1,b2,b3,b4
    unL = uL*anx + vL*any + wL*anz
    unR = uR*anx + vR*any + wR*anz
    c0  = sqrt(max(beta,1.d-30) + 0.5d0*(unL*unL + unR*unR))
    ML  = unL/c0
    MR  = unR/c0
!   AUSM+ convected-speed splitting
    if(abs(ML) .ge. 1.d0) then
      MpL = 0.5d0*(ML + abs(ML))
    else
      MpL = 0.25d0*(ML + 1.d0)**2
    endif
    if(abs(MR) .ge. 1.d0) then
      MmR = 0.5d0*(MR - abs(MR))
    else
      MmR = -0.25d0*(MR - 1.d0)**2
    endif
    Mt  = MpL + MmR
    uf  = c0*Mt
!   upwind convected vector Phi=(beta,u,v,w)
    if(uf .ge. 0.d0) then
      b1 = beta; b2 = uL; b3 = vL; b4 = wL
    else
      b1 = beta; b2 = uR; b3 = vR; b4 = wR
    endif
!   AUSM+ pressure splitting
    if(abs(ML) .ge. 1.d0) then
      ppL = 0.5d0*(1.d0 + sign(1.d0,ML))
    else
      ppL = 0.25d0*(ML + 1.d0)**2*(2.d0 - ML)
    endif
    if(abs(MR) .ge. 1.d0) then
      pmR = 0.5d0*(1.d0 - sign(1.d0,MR))
    else
      pmR = 0.25d0*(MR - 1.d0)**2*(2.d0 + MR)
    endif
    pf = ppL*qL + pmR*qR
    F(1) = uf*b1
    F(2) = uf*b2 + pf*anx
    F(3) = uf*b3 + pf*any
    F(4) = uf*b4 + pf*anz
   end subroutine ac_flux_ausm

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
! flux (AC_Flux: 1=Rusanov, 2=Steger-Warming, 3=AUSM+).  RAC(L) -= F*A ; RAC(R) += F*A.
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
    call ac_recon_face(XW(1,im2,jm2,km2),XW(1,im1,jm1,km1),XW(1,iL,jL,kL), &
                       XW(1,iR,jR,kR),XW(1,ip2,jp2,kp2), qL,qR)
    call ac_recon_face(XW(2,im2,jm2,km2),XW(2,im1,jm1,km1),XW(2,iL,jL,kL), &
                       XW(2,iR,jR,kR),XW(2,ip2,jp2,kp2), uL,uR)
    call ac_recon_face(XW(3,im2,jm2,km2),XW(3,im1,jm1,km1),XW(3,iL,jL,kL), &
                       XW(3,iR,jR,kR),XW(3,ip2,jp2,kp2), vL,vR)
    call ac_recon_face(XW(4,im2,jm2,km2),XW(4,im1,jm1,km1),XW(4,iL,jL,kL), &
                       XW(4,iR,jR,kR),XW(4,ip2,jp2,kp2), wL,wR)
    if(AC_Flux .eq. 2) then
      call ac_flux_sw(qL,uL,vL,wL, qR,uR,vR,wR, nx1,ny1,nz1, beta, FLX)
    else if(AC_Flux .eq. 3) then
      call ac_flux_ausm(qL,uL,vL,wL, qR,uR,vR,wR, nx1,ny1,nz1, beta, FLX)
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
   end interface
   B => Mesh(nMesh)%Block(mBlock)

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
       if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
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


!==============================================================================
! Porous-media (incompressible) flow solver
! SIMPLE algorithm on collocated rectangular grids (non-uniform allowed).
! Local-thermal-non-equilibrium (LTNE) two-temperature model:
!   - fluid phase:  U(1..5) = [rho, u, v, w, Tf],  pressure in B%p
!   - solid frame:  B%Ts
! Difference vs. the low-speed solver: Darcy + Forchheimer drag in the
! momentum equations, volumetric inter-phase heat exchange h_v*(Ts-Tf) in the
! fluid energy equation, plus the solid-frame energy equation with (1-eps)*k_s.
! Steady state, no turbulence model.  Fluid properties = low-speed LS_* set.
!==============================================================================
  module porous_work
   use precision_EC
   implicit none
   real(PRE_EC), allocatable, dimension(:,:,:) :: Fi, Fj, Fk   ! face mass fluxes (rho*u_n*A)
   real(PRE_EC), allocatable, dimension(:,:,:) :: apu, apv, apw ! momentum diagonal coeffs (no relaxation)
   real(PRE_EC), allocatable, dimension(:,:,:) :: su_nb, sv_nb, sw_nb ! sum of neighbor coeffs (for SIMPLEC)
   real(PRE_EC), allocatable, dimension(:,:,:) :: du, dv, dw    ! Vol/ap (for Rhie-Chow & velocity correction)
   real(PRE_EC), allocatable, dimension(:,:,:) :: pp            ! pressure correction p'
   real(PRE_EC), allocatable, dimension(:,:,:) :: conv_src      ! high-order convection correction (deferred)
   integer :: nxw, nyw, nzw
  end module porous_work

!==============================================================================
! Initialize a porous block with primitive variables (fluid + solid frame)
!==============================================================================
  subroutine init_porous_block(B)
    use Global_Var
    implicit none
    Type (Block_TYPE),pointer:: B
    integer:: i,j,k
    do k=1-LAP,B%nz+LAP-1
    do j=1-LAP,B%ny+LAP-1
    do i=1-LAP,B%nx+LAP-1
      B%U(1,i,j,k)=LS_rho
      B%U(2,i,j,k)=0.d0
      B%U(3,i,j,k)=0.d0
      B%U(4,i,j,k)=0.d0
      B%U(5,i,j,k)=LS_T_ref
      B%p(i,j,k)=LS_P_out
      B%Ts(i,j,k)=Porous_T_ref
    enddo; enddo; enddo
  end subroutine init_porous_block

!==============================================================================
! Low-speed boundary conditions: fill ghost cells (primitive variables + p)
! and set boundary face mass fluxes (wall/symmetry = 0, inlet/outlet/farfield).
!==============================================================================
  subroutine porous_boundary_block(nMesh, mBlock)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer:: i,j,k,n,ksub,face_s,ib,ie,jb,je,kb,ke,i1,j1,k1,i2,j2,k2
   real(PRE_EC):: uin_x, uin_y, uin_z, Ain_tot, rho

   interface
     subroutine porous_ghost(bc, face_s, ig,jg,kg, i2,j2,k2, B, uin_x,uin_y,uin_z)
       use precision_EC
       use Global_Var
       integer:: bc, face_s, ig,jg,kg, i2,j2,k2
       Type (Block_TYPE),pointer:: B
       real(PRE_EC):: uin_x,uin_y,uin_z
     end subroutine porous_ghost
     subroutine porous_face_flux_bc(bc, face_s, i1,j1,k1, B, uin_x,uin_y,uin_z)
       use precision_EC
       use Global_Var
       integer:: bc, face_s, i1,j1,k1
       Type (Block_TYPE),pointer:: B
       real(PRE_EC):: uin_x,uin_y,uin_z
     end subroutine porous_face_flux_bc
      subroutine porous_inlet_velocity(B, uin_x, uin_y, uin_z)
        use precision_EC
        use Global_Var
        Type (Block_TYPE),pointer:: B
        real(PRE_EC):: uin_x,uin_y,uin_z
      end subroutine porous_inlet_velocity
   end interface

   B=>Mesh(nMesh)%Block(mBlock)
   rho = LS_rho

!  ensure face-flux work arrays are allocated (boundary may be called before solver)
   if(.not. allocated(Fi) .or. nxw /= B%nx .or. nyw /= B%ny .or. nzw /= B%nz) then
     if(allocated(Fi)) deallocate(Fi,Fj,Fk,apu,apv,apw,su_nb,sv_nb,sw_nb,du,dv,dw,pp,conv_src)
     allocate(Fi(B%nx,B%ny,B%nz), Fj(B%nx,B%ny,B%nz), Fk(B%nx,B%ny,B%nz))
     allocate(apu(B%nx,B%ny,B%nz), apv(B%nx,B%ny,B%nz), apw(B%nx,B%ny,B%nz))
     allocate(su_nb(B%nx,B%ny,B%nz), sv_nb(B%nx,B%ny,B%nz), sw_nb(B%nx,B%ny,B%nz))
     allocate(du(B%nx,B%ny,B%nz), dv(B%nx,B%ny,B%nz), dw(B%nx,B%ny,B%nz))
     allocate(pp(0:B%nx,0:B%ny,0:B%nz))
     allocate(conv_src(B%nx,B%ny,B%nz))
     nxw=B%nx; nyw=B%ny; nzw=B%nz
     Fi=0.d0; Fj=0.d0; Fk=0.d0
   endif

!  Compute inlet velocity components (velocity inlet / mass-flow inlet / pressure inlet)
   call porous_inlet_velocity(B, uin_x, uin_y, uin_z)

   do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle                       ! internal interface -> buffer exchange
     if(associated(B%bc_msg2)) then
       if(B%bc_msg2(ksub)%bc < 0) cycle        ! interface connection handled elsewhere
     endif

     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke

     select case(face_s)
     case(1)  ! i- face
       do n=1, LAP
         i1 = ib - n          ! ghost
         i2 = ib + n - 1      ! interior
         do k=kb,ke-1; do j=jb,je-1
           call porous_ghost(Bc%bc, face_s, i1,j,k, i2,j,k, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do k=kb,ke-1; do j=jb,je-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=ib, j1=j, k1=k, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     case(4)  ! i+ face
       do n=1, LAP
         i1 = ie + n - 1
         i2 = ie - n
         do k=kb,ke-1; do j=jb,je-1
           call porous_ghost(Bc%bc, face_s, i1,j,k, i2,j,k, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do k=kb,ke-1; do j=jb,je-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=ie, j1=j, k1=k, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     case(2)  ! j- face
       do n=1, LAP
         j1 = jb - n
         j2 = jb + n - 1
         do k=kb,ke-1; do i=ib,ie-1
           call porous_ghost(Bc%bc, face_s, i,j1,k, i,j2,k, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do k=kb,ke-1; do i=ib,ie-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=i, j1=jb, k1=k, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     case(5)  ! j+ face
       do n=1, LAP
         j1 = je + n - 1
         j2 = je - n
         do k=kb,ke-1; do i=ib,ie-1
           call porous_ghost(Bc%bc, face_s, i,j1,k, i,j2,k, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do k=kb,ke-1; do i=ib,ie-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=i, j1=je, k1=k, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     case(3)  ! k- face
       do n=1, LAP
         k1 = kb - n
         k2 = kb + n - 1
         do j=jb,je-1; do i=ib,ie-1
           call porous_ghost(Bc%bc, face_s, i,j,k1, i,j,k2, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do j=jb,je-1; do i=ib,ie-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=i, j1=j, k1=kb, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     case(6)  ! k+ face
       do n=1, LAP
         k1 = ke + n - 1
         k2 = ke - n
         do j=jb,je-1; do i=ib,ie-1
           call porous_ghost(Bc%bc, face_s, i,j,k1, i,j,k2, B, uin_x,uin_y,uin_z)
         enddo; enddo
       enddo
       do j=jb,je-1; do i=ib,ie-1
         call porous_face_flux_bc(Bc%bc, face_s, i1=i, j1=j, k1=ke, B=B, &
                                    uin_x=uin_x, uin_y=uin_y, uin_z=uin_z)
       enddo; enddo
     end select
   enddo
  end subroutine porous_boundary_block

!==============================================================================
! Compute inlet velocity components according to LS_Inlet_Type
!   1 = velocity inlet (LS_U_in, LS_V_in, LS_W_in)
!   2 = mass-flow inlet (normal velocity = mdot/(rho*A), tangential = 0)
!   3 = pressure inlet (velocity is extrapolated, not used here)
!==============================================================================
  subroutine porous_inlet_velocity(B, uin_x, uin_y, uin_z)
   use Global_Var
   implicit none
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: uin_x, uin_y, uin_z
    integer:: ksub, face_s, ib,ie,jb,je,kb,ke, i,j,k
    real(PRE_EC):: A_tot
   TYPE (BC_MSG_TYPE),pointer:: Bc

   uin_x = LS_U_in; uin_y = LS_V_in; uin_z = LS_W_in

   if(LS_Inlet_Type == 2) then
!    mass-flow inlet: total inlet area then normal velocity
     A_tot = 0.d0
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc == BC_Inflow .or. Bc%bc == BC_LS_Inlet) then
         face_s = Bc%face
         ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke
         select case(face_s)
         case(1,4)
           do k=kb,ke-1; do j=jb,je-1
             A_tot = A_tot + B%Si(ib,j,k)
           enddo; enddo
         case(2,5)
           do k=kb,ke-1; do i=ib,ie-1
             A_tot = A_tot + B%Sj(i,jb,k)
           enddo; enddo
         case(3,6)
           do j=jb,je-1; do i=ib,ie-1
             A_tot = A_tot + B%Sk(i,j,kb)
           enddo; enddo
         end select
       endif
     enddo
     if(A_tot > 1.d-30) then
       uin_x = LS_Mdot_in / (LS_rho * A_tot)
     else
       uin_x = 0.d0
     endif
     uin_y = 0.d0; uin_z = 0.d0
   endif
  end subroutine porous_inlet_velocity

!==============================================================================
! Fill one ghost cell (ig,jg,kg) from interior cell (i2,j2,k2) for a given
! boundary type. face_s: 1=i-, 2=j-, 3=k-, 4=i+, 5=j+, 6=k+
!==============================================================================
  subroutine porous_ghost(bc, face_s, ig,jg,kg, i2,j2,k2, B, uin_x,uin_y,uin_z)
   use Global_Var
   use const_var
   implicit none
   integer:: bc, face_s, ig,jg,kg, i2,j2,k2
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: uin_x, uin_y, uin_z
   real(PRE_EC):: uw, vw, ww    ! wall velocity components (lid-driven)
    real(PRE_EC):: Tw_val, Qw_val, dx_g, htc_val, Tinf_val, ks_eff
    integer:: ii
    logical:: found_bc


   interface
     subroutine porous_set_face_speed(face_s, ig,jg,kg, i2,j2,k2, B, wall_mode, uw,vw,ww)
       use precision_EC
       use Global_Var
       integer:: face_s, ig,jg,kg, i2,j2,k2
       Type (Block_TYPE),pointer:: B
       logical:: wall_mode
       real(PRE_EC):: uw, vw, ww
     end subroutine porous_set_face_speed

   end interface

!  Look up a per-face thermal BC (from solid_bc.inp, reused for porous blocks).
!  Tw>0    -> isothermal at Tw (K);
!  Tw<0    -> heat flux Qw (W/m2);
!  Tw==0 & htc>0 -> convective (Robin) q = htc*(T - Tinf).
!  The BC is applied to the SOLID-FRAME temperature Ts on ANY physical face
!  (wall / inlet / outlet / symmetry) that carries a solid_bc.inp entry, so the
!  frame and the fluid can have independent boundary conditions at one face.
   found_bc = .false.; Tw_val = 0.d0; Qw_val = 0.d0
   htc_val = 0.d0; Tinf_val = 0.d0
   if(associated(B%solid_bc_face_no)) then
     do ii=1, B%solid_bc_nface
       if(B%solid_bc_face_no(ii) == face_s) then
         found_bc = .true.
         Tw_val   = B%solid_bc_Tw(ii)
         Qw_val   = B%solid_bc_Qw(ii)
         htc_val  = B%solid_bc_htc(ii)
         Tinf_val = B%solid_bc_Tinf(ii)
         exit
       endif
     enddo
   endif

!  Wall tangential velocity (lid-driven: j+ face moves at LS_U_lid in x)
   uw = 0.d0; vw = 0.d0; ww = 0.d0
   if(bc == BC_Wall .and. face_s == 5) uw = LS_U_lid

   select case(bc)
   case(BC_Wall)   ! no-slip wall (possibly lid-driven)
     call porous_set_face_speed(face_s, ig,jg,kg, i2,j2,k2, B, &
          .true., uw, vw, ww)
!    temperature: adiabatic (T_g = T_i) or isothermal (T_g = 2*Tw - T_i)
     if(LS_T_wall > 0.d0) then
       B%U(5,ig,jg,kg) = 2.d0*LS_T_wall - B%U(5,i2,j2,k2)
     else
       B%U(5,ig,jg,kg) = B%U(5,i2,j2,k2)
     endif
     B%p(ig,jg,kg) = B%p(i2,j2,k2)
   case(BC_Symmetry)  ! symmetry: normal velocity mirror, tangential copy
     call porous_set_face_speed(face_s, ig,jg,kg, i2,j2,k2, B, &
          .false., uw, vw, ww)
     B%U(5,ig,jg,kg) = B%U(5,i2,j2,k2)
     B%p(ig,jg,kg) = B%p(i2,j2,k2)
   case(BC_Inflow, BC_LS_Inlet)   ! inlet
     if(LS_Inlet_Type == 3) then
!      pressure inlet: velocity extrapolated, pressure specified
       B%U(2,ig,jg,kg) = B%U(2,i2,j2,k2)
       B%U(3,ig,jg,kg) = B%U(3,i2,j2,k2)
       B%U(4,ig,jg,kg) = B%U(4,i2,j2,k2)
       B%p(ig,jg,kg) = 2.d0*LS_P_in - B%p(i2,j2,k2)
     else
!      velocity / mass-flow inlet: Dirichlet velocity, zero pressure gradient
       B%U(2,ig,jg,kg) = 2.d0*uin_x - B%U(2,i2,j2,k2)
       B%U(3,ig,jg,kg) = 2.d0*uin_y - B%U(3,i2,j2,k2)
       B%U(4,ig,jg,kg) = 2.d0*uin_z - B%U(4,i2,j2,k2)
       B%p(ig,jg,kg) = B%p(i2,j2,k2)
     endif
     B%U(5,ig,jg,kg) = 2.d0*LS_T_ref - B%U(5,i2,j2,k2)
   case(BC_Outflow, BC_LS_Outlet)  ! pressure outlet
     B%U(2,ig,jg,kg) = B%U(2,i2,j2,k2)
     B%U(3,ig,jg,kg) = B%U(3,i2,j2,k2)
     B%U(4,ig,jg,kg) = B%U(4,i2,j2,k2)
     B%U(5,ig,jg,kg) = B%U(5,i2,j2,k2)
     B%p(ig,jg,kg) = 2.d0*LS_P_out - B%p(i2,j2,k2)
   case default     ! farfield / others: zero gradient
     B%U(2,ig,jg,kg) = B%U(2,i2,j2,k2)
     B%U(3,ig,jg,kg) = B%U(3,i2,j2,k2)
     B%U(4,ig,jg,kg) = B%U(4,i2,j2,k2)
     B%U(5,ig,jg,kg) = B%U(5,i2,j2,k2)
     B%p(ig,jg,kg) = B%p(i2,j2,k2)
   end select


!  --- porous solid-frame temperature ghost ---
!  Adiabatic (zero-gradient) by default.  A face with a solid_bc.inp entry
!  honours the per-phase thermal BC on the frame:
!    isothermal : Ts_g = 2*Tw - Ts_int
!    heat flux  : Ts_g = Ts_int + Qw*dx/kse  (linear ghost extrapolation)
!    convective : kse*(Ts_int - Ts_g)/dx = htc*((Ts_int+Ts_g)/2 - Tinf)
!  with kse = (1-eps)*k_s (effective skeleton conductivity).
!  (dx = distance between the interior and ghost cell centres).
   B%Ts(ig,jg,kg) = B%Ts(i2,j2,k2)
   if(found_bc) then
     if(Tw_val > 0.d0) then
       B%Ts(ig,jg,kg) = 2.d0*Tw_val - B%Ts(i2,j2,k2)
     else if(Tw_val == 0.d0 .and. htc_val > 0.d0 .and. Qw_val == 0.d0) then
       dx_g = sqrt( (B%xc(ig,jg,kg)-B%xc(i2,j2,k2))**2 &
                  + (B%yc(ig,jg,kg)-B%yc(i2,j2,k2))**2 &
                  + (B%zc(ig,jg,kg)-B%zc(i2,j2,k2))**2 ) * Lscale
       ! effective skeleton conductivity (Ts eq. is in total-volume units):
       ks_eff = max((1.d0 - B%porous_eps)*B%solid_k, 1.d-30)
       B%Ts(ig,jg,kg) = ( (ks_eff - 0.5d0*htc_val*dx_g)*B%Ts(i2,j2,k2) &
                        + htc_val*dx_g*Tinf_val ) &
                      / max(ks_eff + 0.5d0*htc_val*dx_g, 1.d-30)
     else if(Qw_val /= 0.d0) then
       dx_g = sqrt( (B%xc(ig,jg,kg)-B%xc(i2,j2,k2))**2 &
                  + (B%yc(ig,jg,kg)-B%yc(i2,j2,k2))**2 &
                  + (B%zc(ig,jg,kg)-B%zc(i2,j2,k2))**2 ) * Lscale
       ks_eff = max((1.d0 - B%porous_eps)*B%solid_k, 1.d-30)
       B%Ts(ig,jg,kg) = B%Ts(i2,j2,k2) + Qw_val*dx_g/ks_eff
     endif
   endif
!  density is constant
   B%U(1,ig,jg,kg) = LS_rho
  end subroutine porous_ghost

!==============================================================================
! Clip primitive-variable fields to physically plausible ranges.
! Purely-diagnostic safety net: SIMPLE transients with deferred correction
! can occasionally produce 10^100 spikes in isolated cells at high Re /
! fine grids; clipping at �10U_scale keeps subsequent steps sane rather
! than cascading to Inf/NaN.  T is clipped around the expected reference
! so pure-convection (k=0) cannot accumulate a single outlier Fi=10^100
! and poison the upwind stabilisation.
!==============================================================================
  subroutine porous_clip_fields(nMesh, mBlock)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz
   real(PRE_EC):: U_scale, T_lo, T_hi, p_scale
   ! Very wide hard limits: only truncate genuine SIMPLE overshoot that
   ! would cascade to Inf/NaN (10^100).  They never touch a well-tuned
   ! solution but keep isolated overflow from killing the run.
   U_scale = 1.d6
   T_lo = 1.d0;     T_hi = 1.d5
   p_scale = 1.d8
   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   do k=0,nz
   do j=0,ny
   do i=0,nx
     if(B%U(2,i,j,k) >  U_scale) B%U(2,i,j,k)= U_scale
     if(B%U(2,i,j,k) < -U_scale) B%U(2,i,j,k)=-U_scale
     if(B%U(3,i,j,k) >  U_scale) B%U(3,i,j,k)= U_scale
     if(B%U(3,i,j,k) < -U_scale) B%U(3,i,j,k)=-U_scale
     if(B%U(4,i,j,k) >  U_scale) B%U(4,i,j,k)= U_scale
     if(B%U(4,i,j,k) < -U_scale) B%U(4,i,j,k)=-U_scale
     if(B%U(5,i,j,k) > T_hi) B%U(5,i,j,k)=T_hi
     if(B%U(5,i,j,k) < T_lo) B%U(5,i,j,k)=T_lo
     if(B%p(i,j,k) >  p_scale) B%p(i,j,k)= p_scale
     if(B%p(i,j,k) >  p_scale) B%p(i,j,k)= p_scale
    if(B%p(i,j,k) < -p_scale) B%p(i,j,k)=-p_scale
    if(B%Ts(i,j,k) > T_hi) B%Ts(i,j,k)=T_hi
    if(B%Ts(i,j,k) < T_lo) B%Ts(i,j,k)=T_lo
   enddo; enddo; enddo
  end subroutine porous_clip_fields

!==============================================================================
! Set velocity ghost cell for wall (no-slip, mirror) or symmetry (normal mirror).
! wall_mode = .true. : no-slip wall (normal & tangential mirror, plus wall speed)
! wall_mode = .false.: symmetry (normal mirror, tangential copy)
!==============================================================================
  subroutine porous_set_face_speed(face_s, ig,jg,kg, i2,j2,k2, B, wall_mode, uw,vw,ww)
   use precision_EC
   use Global_Var
   implicit none
   integer:: face_s, ig,jg,kg, i2,j2,k2
   Type (Block_TYPE),pointer:: B
   logical:: wall_mode
   real(PRE_EC):: uw, vw, ww
   real(PRE_EC):: ui, vi, wi

   ui = B%U(2,i2,j2,k2); vi = B%U(3,i2,j2,k2); wi = B%U(4,i2,j2,k2)

   select case(face_s)
   case(1,4)  ! i-face: normal = x
     if(wall_mode) then
       B%U(2,ig,jg,kg) = 2.d0*uw - ui     ! normal
       B%U(3,ig,jg,kg) = 2.d0*vw - vi
       B%U(4,ig,jg,kg) = 2.d0*ww - wi
     else
       B%U(2,ig,jg,kg) = -ui              ! normal mirror
       B%U(3,ig,jg,kg) =  vi
       B%U(4,ig,jg,kg) =  wi
     endif
   case(2,5)  ! j-face: normal = y
     if(wall_mode) then
       B%U(2,ig,jg,kg) = 2.d0*uw - ui
       B%U(3,ig,jg,kg) = 2.d0*vw - vi
       B%U(4,ig,jg,kg) = 2.d0*ww - wi
     else
       B%U(2,ig,jg,kg) =  ui
       B%U(3,ig,jg,kg) = -vi
       B%U(4,ig,jg,kg) =  wi
     endif
   case(3,6)  ! k-face: normal = z
     if(wall_mode) then
       B%U(2,ig,jg,kg) = 2.d0*uw - ui
       B%U(3,ig,jg,kg) = 2.d0*vw - vi
       B%U(4,ig,jg,kg) = 2.d0*ww - wi
     else
       B%U(2,ig,jg,kg) =  ui
       B%U(3,ig,jg,kg) =  vi
       B%U(4,ig,jg,kg) = -wi
     endif
   end select
  end subroutine porous_set_face_speed

!==============================================================================
! Set boundary face mass flux (wall/symmetry = 0, inlet/outlet/farfield).
! (i1,j1,k1) is the face location (node index); face_s gives orientation.
!==============================================================================
  subroutine porous_face_flux_bc(bc, face_s, i1,j1,k1, B, uin_x,uin_y,uin_z)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: bc, face_s, i1,j1,k1
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: uin_x, uin_y, uin_z
   real(PRE_EC):: A, un

   select case(face_s)
   case(1,4); A = B%Si(i1,j1,k1)
   case(2,5); A = B%Sj(i1,j1,k1)
   case(3,6); A = B%Sk(i1,j1,k1)
   case default; A = 0.d0
   end select

   select case(bc)
   case(BC_Wall, BC_Symmetry)
!    no mass flux through wall / symmetry
     call porous_store_face_flux(face_s, i1,j1,k1, 0.d0)
   case(BC_Inflow, BC_LS_Inlet)
!    inlet normal velocity
     if(LS_Inlet_Type == 3) then
!      pressure inlet: velocity extrapolated (zero gradient)
       select case(face_s)
       case(1);   un =  B%U(2,i1,j1,k1)
       case(4);   un =  B%U(2,i1-1,j1,k1)
       case(2);   un =  B%U(3,i1,j1,k1)
       case(5);   un =  B%U(3,i1,j1-1,k1)
       case(3);   un =  B%U(4,i1,j1,k1)
       case(6);   un =  B%U(4,i1,j1,k1-1)
       case default; un = 0.d0
       end select
       call porous_store_face_flux(face_s, i1,j1,k1, LS_rho*A*un)
     else
!      Fi/Fj/Fk carry the flux in the +coordinate direction (interior faces
!      use ue = +x face velocity), so a velocity inlet must add +uin*A on
!      EVERY face regardless of orientation; the sign of uin decides inflow
!      vs outflow.  The former -uin_x on the i-/j-/k- faces turned an inlet
!      into a suction boundary and drove a spurious reverse Poiseuille flow.
       select case(face_s)
       case(1);   un =  uin_x
       case(4);   un =  uin_x
       case(2);   un =  uin_y
       case(5);   un =  uin_y
       case(3);   un =  uin_z
       case(6);   un =  uin_z
       case default; un = 0.d0
       end select
       call porous_store_face_flux(face_s, i1,j1,k1, LS_rho*A*un)
     endif
   case(BC_Outflow, BC_LS_Outlet, BC_Farfield)
!    outlet / farfield: velocity extrapolated (zero gradient)
     select case(face_s)
     case(1);   un = -B%U(2,i1,j1,k1)
     case(4);   un =  B%U(2,i1-1,j1,k1)
     case(2);   un = -B%U(3,i1,j1,k1)
     case(5);   un =  B%U(3,i1,j1-1,k1)
     case(3);   un = -B%U(4,i1,j1,k1)
     case(6);   un =  B%U(4,i1,j1,k1-1)
     case default; un = 0.d0
     end select
     call porous_store_face_flux(face_s, i1,j1,k1, LS_rho*A*un)
   case default
     call porous_store_face_flux(face_s, i1,j1,k1, 0.d0)
   end select
  end subroutine porous_face_flux_bc

!==============================================================================
! Store a face mass flux into Fi / Fj / Fk
!==============================================================================
  subroutine porous_store_face_flux(face_s, i1,j1,k1, flux)
   use porous_work
   implicit none
   integer:: face_s, i1,j1,k1
   real(PRE_EC):: flux
   select case(face_s)
   case(1,4); Fi(i1,j1,k1) = flux
   case(2,5); Fj(i1,j1,k1) = flux
   case(3,6); Fk(i1,j1,k1) = flux
   end select
  end subroutine porous_store_face_flux

!==============================================================================
! Compute face mass fluxes.
! Interior faces: Rhie-Chow momentum interpolation.
! Boundary faces: already set by porous_boundary_block.
!==============================================================================
  subroutine porous_face_flux(nMesh, mBlock)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz
   real(PRE_EC):: uL,uR,pL,pR,dL,dR,dxe,ue, rho, gpf
   real(PRE_EC):: x_LL, x_L, x_R, x_RR, y_LL, y_L, y_R, y_RR, z_LL, z_L, z_R, z_RR

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   rho = LS_rho

!  i-direction internal faces (i = 2 .. nx-1)
   do k=1,nz-1
   do j=1,ny-1
   do i=2,nx-1
     uL=B%U(2,i-1,j,k); uR=B%U(2,i,j,k)
     pL=B%p(i-1,j,k);   pR=B%p(i,j,k)
     dL=du(i-1,j,k);    dR=du(i,j,k)
     x_LL = B%xc(i-2,j,k); x_L = B%xc(i-1,j,k)
     x_R  = B%xc(i,j,k);   x_RR= B%xc(i+1,j,k)
     dxe = x_R - x_L
     if(abs(dxe) < 1.d-30) dxe = 1.d-30
!    4-point (i-2:i+1) quadratic pressure gradient at face i (midpoint of x_L,x_R).
!    Exact for quadratic p(x) on arbitrary non-uniform Cartesian grids.
     call porous_grad_face_4pt(pL, B%p(i-2,j,k), pR, B%p(i+1,j,k), x_L, x_LL, x_R, x_RR, gpf)
     ue = 0.5d0*(uL+uR) + 0.5d0*(dL+dR)*( (pL-pR)/dxe - gpf )
     Fi(i,j,k) = rho*B%Si(i,j,k)*ue
   enddo; enddo; enddo

!  j-direction internal faces (j = 2 .. ny-1)
   do k=1,nz-1
   do j=2,ny-1
   do i=1,nx-1
     uL=B%U(3,i,j-1,k); uR=B%U(3,i,j,k)
     pL=B%p(i,j-1,k);   pR=B%p(i,j,k)
     dL=dv(i,j-1,k);    dR=dv(i,j,k)
     y_LL = B%yc(i,j-2,k); y_L = B%yc(i,j-1,k)
     y_R  = B%yc(i,j,k);   y_RR= B%yc(i,j+1,k)
     dxe = y_R - y_L
     if(abs(dxe) < 1.d-30) dxe = 1.d-30
     call porous_grad_face_4pt(pL, B%p(i,j-2,k), pR, B%p(i,j+1,k), y_L, y_LL, y_R, y_RR, gpf)
     ue = 0.5d0*(uL+uR) + 0.5d0*(dL+dR)*( (pL-pR)/dxe - gpf )
     Fj(i,j,k) = rho*B%Sj(i,j,k)*ue
   enddo; enddo; enddo

!  k-direction internal faces (k = 2 .. nz-1)
   do k=2,nz-1
   do j=1,ny-1
   do i=1,nx-1
     uL=B%U(4,i,j,k-1); uR=B%U(4,i,j,k)
     pL=B%p(i,j,k-1);   pR=B%p(i,j,k)
     dL=dw(i,j,k-1);    dR=dw(i,j,k)
     z_LL = B%zc(i,j,k-2); z_L = B%zc(i,j,k-1)
     z_R  = B%zc(i,j,k);   z_RR= B%zc(i,j,k+1)
     dxe = z_R - z_L
     if(abs(dxe) < 1.d-30) dxe = 1.d-30
     call porous_grad_face_4pt(pL, B%p(i,j,k-2), pR, B%p(i,j,k+1), z_L, z_LL, z_R, z_RR, gpf)
     ue = 0.5d0*(uL+uR) + 0.5d0*(dL+dR)*( (pL-pR)/dxe - gpf )
     Fk(i,j,k) = rho*B%Sk(i,j,k)*ue
   enddo; enddo; enddo
  end subroutine porous_face_flux

!------------------------------------------------------------------------------
! 3-point Lagrange derivative at the center cell x_C, given cell values
!   phi_L at x_L (left neighbour), phi_C at x_C, phi_R at x_R (right nbr).
! Second-order accurate on arbitrary non-uniform grids; collapses to the
! standard central difference (phi_R-phi_L)/(x_R-x_L) when x_C = 0.5(x_L+x_R).
!------------------------------------------------------------------------------
  subroutine porous_grad_cell_3pt(phi_L, phi_C, phi_R, x_L, x_C, x_R, grad)
   use const_var, only: PRE_EC
   implicit none
   real(PRE_EC), intent(in) :: phi_L, phi_C, phi_R, x_L, x_C, x_R
   real(PRE_EC), intent(out):: grad
   real(PRE_EC):: dxL, dxR, dxT, aL, aR
   dxL = x_C - x_L;   dxR = x_R - x_C
   dxT = x_R - x_L
   if (abs(dxT) < 1.d-30) then
     grad = 0.d0
     return
   end if
   aL = dxR / dxT
   aR = dxL / dxT
   grad = aL*(phi_C - phi_L)/max(dxL,1.d-30) + aR*(phi_R - phi_C)/max(dxR,1.d-30)
  end subroutine porous_grad_cell_3pt

  subroutine porous_grad_face_4pt(p_L, p_LL, p_R, p_RR, x_L, x_LL, x_R, x_RR, grad)
   use const_var, only: PRE_EC
   implicit none
   real(PRE_EC), intent(in) :: p_L, p_LL, p_R, p_RR, x_L, x_LL, x_R, x_RR
   real(PRE_EC), intent(out):: grad
   real(PRE_EC):: gL, gR, dxL, dxR, invs
   call porous_grad_cell_3pt(p_LL, p_L, p_R,  x_LL, x_L, x_R, gL)
   call porous_grad_cell_3pt(p_L,  p_R, p_RR, x_L,  x_R, x_RR, gR)
!  Distance from cell L/R centers to the face midpoint. For a uniform grid
!  dxL=dxR so the blend is 0.5*(gL+gR) = same as the previous arithmetic
!  average; for non-uniform grids it properly weights the nearer cell more.
   dxL  = 0.5d0*(x_L + x_R) - x_L
   dxR  = x_R - 0.5d0*(x_L + x_R)
   invs = 1.d0 / max(dxL + dxR, 1.d-30)
   grad = (dxR*gL + dxL*gR) * invs
  end subroutine porous_grad_face_4pt

!==============================================================================
! High-order face correction delta = phi_hi - phi_upwind for ONE face
! (deferred-correction term). ph_up = upwind value, ph_up2 = second upwind,
! ph_dn = downwind value.
!   LS_Scheme = 2 : 2nd-order upwind (SOU): delta = 0.5*(ph_up - ph_up2)
!   LS_Scheme = 3 : MUSCL with Van Leer limiter
!==============================================================================
  subroutine porous_conv_face_delta(ph_up, ph_up2, ph_dn, delta)
   use Global_Var
   implicit none
   real(PRE_EC):: ph_up, ph_up2, ph_dn, delta
   real(PRE_EC):: r, psi
   if(LS_Scheme == 2) then
     delta = 0.5d0*(ph_up - ph_up2)
   else if(LS_Scheme == 3) then
     if(abs(ph_dn-ph_up) < 1.d-30) then
       delta = 0.d0
     else
       r = (ph_up - ph_up2)/(ph_dn - ph_up)
       psi = (r + abs(r))/(1.d0 + abs(r))
       delta = 0.5d0*psi*(ph_dn - ph_up)
     endif
   else
     delta = 0.d0
   endif
  end subroutine porous_conv_face_delta

!==============================================================================
! High-order convection correction source (deferred correction) for variable
! ivar (2=u, 3=v, 4=w, 5=T). Result stored in module array conv_src.
! Implicit part keeps 1st-order upwind (diagonal dominance); the difference
! between high-order and 1st-order face values is added explicitly here.
!==============================================================================
  subroutine porous_conv_high(nMesh, mBlock, ivar)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock, ivar
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz
   real(PRE_EC):: cf, Fe,Fw,Fn,Fs,Ft,Fb, conv, delta
   real(PRE_EC):: phP,phE,phW,phN,phS,phT,phB
   real(PRE_EC):: phEE,phWW,phNN,phSS,phTT,phBB

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   cf = 1.d0
   if(ivar == 5) cf = LS_Cp
   conv_src = 0.d0

   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
     phP = B%U(ivar,i,j,k)
     phE = B%U(ivar,i+1,j,k); phW = B%U(ivar,i-1,j,k)
     phN = B%U(ivar,i,j+1,k); phS = B%U(ivar,i,j-1,k)
     phT = B%U(ivar,i,j,k+1); phB = B%U(ivar,i,j,k-1)
     phEE= B%U(ivar,i+2,j,k); phWW= B%U(ivar,i-2,j,k)
     phNN= B%U(ivar,i,j+2,k); phSS= B%U(ivar,i,j-2,k)
     phTT= B%U(ivar,i,j,k+2); phBB= B%U(ivar,i,j,k-2)

     conv = 0.d0
!    e face
     Fe = cf*Fi(i+1,j,k)
     if(Fe > 0.d0) then
       call porous_conv_face_delta(phP, phW, phE, delta)
       conv = conv + Fe*delta
     else
       call porous_conv_face_delta(phE, phEE, phP, delta)
       conv = conv + Fe*delta
     endif
!    w face (outward normal is -x: accumulate with outward flux -Fw)
     Fw = cf*Fi(i,j,k)
     if(Fw > 0.d0) then
       call porous_conv_face_delta(phW, phWW, phP, delta)
       conv = conv - Fw*delta
     else
       call porous_conv_face_delta(phP, phE, phW, delta)
       conv = conv - Fw*delta
     endif
!    n face
     Fn = cf*Fj(i,j+1,k)
     if(Fn > 0.d0) then
       call porous_conv_face_delta(phP, phS, phN, delta)
       conv = conv + Fn*delta
     else
       call porous_conv_face_delta(phN, phNN, phP, delta)
       conv = conv + Fn*delta
     endif
!    s face (outward normal is -y: accumulate with outward flux -Fs)
     Fs = cf*Fj(i,j,k)
     if(Fs > 0.d0) then
       call porous_conv_face_delta(phS, phSS, phP, delta)
       conv = conv - Fs*delta
     else
       call porous_conv_face_delta(phP, phN, phS, delta)
       conv = conv - Fs*delta
     endif
!    t face
     Ft = cf*Fk(i,j,k+1)
     if(Ft > 0.d0) then
       call porous_conv_face_delta(phP, phB, phT, delta)
       conv = conv + Ft*delta
     else
       call porous_conv_face_delta(phT, phTT, phP, delta)
       conv = conv + Ft*delta
     endif
!    b face (outward normal is -z: accumulate with outward flux -Fb)
     Fb = cf*Fk(i,j,k)
     if(Fb > 0.d0) then
       call porous_conv_face_delta(phB, phBB, phP, delta)
       conv = conv - Fb*delta
     else
       call porous_conv_face_delta(phP, phT, phB, delta)
       conv = conv - Fb*delta
     endif
     conv_src(i,j,k) = conv
   enddo; enddo; enddo
  end subroutine porous_conv_high

!==============================================================================
! Solve one momentum component by Gauss-Seidel with first-order upwind
! convection and central-difference diffusion. dir = 1(x), 2(y), 3(z).
! Stores diagonal (apu/apv/apw) and du/dv/dw = Vol/ap for Rhie-Chow.
!==============================================================================
  subroutine porous_momentum(nMesh, mBlock, dir)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock, dir
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz, iter, it
    real(PRE_EC):: rho, mu, ap, aE,aW,aN,aS,aT,aB, src, De,Dwe,Dn,Ds,Dt,Db
   real(PRE_EC):: Fe,Fw,Fn,Fs,Ft,Fb, gx, vol, unew, uold, alpha, denom
       real(PRE_EC):: eps_p, porK, Fcoef_p, vm, sp
integer,parameter:: INNER=3

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   rho=LS_rho; mu=LS_mu; alpha=LS_alpha_u

!  high-order convection correction (deferred), computed once per sweep
   if(LS_Scheme > 1) call porous_conv_high(nMesh, mBlock, dir+1)
   do it=1, INNER
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       vol = B%Vol(i,j,k)

!      diffusion coefficients (face area / center distance)
       De = mu*B%Si(i+1,j,k)/max(B%xc(i+1,j,k)-B%xc(i,j,k), 1.d-30)
       Dwe = mu*B%Si(i,j,k)/max(B%xc(i,j,k)-B%xc(i-1,j,k), 1.d-30)
       Dn = mu*B%Sj(i,j+1,k)/max(B%yc(i,j+1,k)-B%yc(i,j,k), 1.d-30)
       Ds = mu*B%Sj(i,j,k)/max(B%yc(i,j,k)-B%yc(i,j-1,k), 1.d-30)
       Dt = mu*B%Sk(i,j,k+1)/max(B%zc(i,j,k+1)-B%zc(i,j,k), 1.d-30)
       Db = mu*B%Sk(i,j,k)/max(B%zc(i,j,k)-B%zc(i,j,k-1), 1.d-30)

!      convective fluxes
       Fe = Fi(i+1,j,k); Fw = Fi(i,j,k)
       Fn = Fj(i,j+1,k); Fs = Fj(i,j,k)
       Ft = Fk(i,j,k+1); Fb = Fk(i,j,k)

!      upwind coefficients
       aE = De + max(-Fe, 0.d0)
       aW = Dwe + max( Fw, 0.d0)
       aN = Dn + max(-Fn, 0.d0)
       aS = Ds + max( Fs, 0.d0)
       aT = Dt + max(-Ft, 0.d0)
       aB = Db + max( Fb, 0.d0)
       ap = aE+aW+aN+aS+aT+aB + (Fe-Fw+Fn-Fs+Ft-Fb)

!      Porous drag (Darcy + Forchheimer), volume-averaged momentum:
!        S_u = -(eps^2*mu/K)*u - (eps^3*F*rho*|V|/sqrt(K))*u
!      Linearised with a negative slope added to a_p; du=Vol/ap and
!      the Rhie-Chow / pressure-correction steps inherit the reduced d.
        eps_p   = B%porous_eps
        porK    = B%porous_dp**2 * eps_p**3 / max(150.d0*(1.d0-eps_p)**2, 1.d-30)
        Fcoef_p = 1.75d0 / max(sqrt(150.d0)*eps_p**1.5d0, 1.d-30)
        vm      = sqrt(B%U(2,i,j,k)**2 + B%U(3,i,j,k)**2 + B%U(4,i,j,k)**2)
!       per-unit-volume drag coefficient, integrated over the cell (like the
!       pressure-gradient source -gx*vol) before entering the a_p diagonal
        sp      = (eps_p*eps_p*mu/max(porK,1.d-30) &
                + eps_p**3 * Fcoef_p * rho * vm / max(sqrt(porK),1.d-30)) * vol
        ap = ap + sp

!      pressure gradient source. 3-point Lagrange derivative on arbitrary
!      non-uniform cell centers (p_{i-1},p_i,p_{i+1}) evaluated at cell i
!      center, formally second-order; collapses to central difference on
!      uniform grids (factor 1 exact).
       if(dir == 1) then
         call porous_grad_cell_3pt(B%p(i-1,j,k), B%p(i,j,k), B%p(i+1,j,k), &
                            B%xc(i-1,j,k),B%xc(i,j,k),B%xc(i+1,j,k), gx)
         src = -gx*vol
         if(LS_Scheme > 1) src = src - conv_src(i,j,k)
         uold = B%U(2,i,j,k)
         unew = (aE*B%U(2,i+1,j,k)+aW*B%U(2,i-1,j,k)+aN*B%U(2,i,j+1,k)+aS*B%U(2,i,j-1,k) &
                +aT*B%U(2,i,j,k+1)+aB*B%U(2,i,j,k-1)+src)/max(ap,1.d-30)
         B%U(2,i,j,k) = uold + alpha*(unew-uold)
         apu(i,j,k) = ap
         su_nb(i,j,k) = aE+aW+aN+aS+aT+aB
         if(LS_Algorithm == 2) then
!          SIMPLEC d factor: d_C = Vol*alpha_u/(ap - alpha_u*sum_a_nb).
!          The pressure update below is FULL (alpha_p=1), so no pressure
!          under-relaxation is needed (Vandoormaal & Raithby 1984).
           denom = ap - alpha*su_nb(i,j,k)
           if(denom < 0.1d0*ap) denom = 0.1d0*ap   ! floor: d cannot blow up
!          cap d at 3x the SIMPLE value: in cells with net inflow (ap close to
!          sum_a_nb) the raw SIMPLEC d can be 5-10x larger, which amplifies the
!          Rhie-Chow pressure-gradient truncation error and blows up Fi.
           du(i,j,k) = min(vol*alpha/max(denom, 1.d-30), 3.d0*vol/max(ap,1.d-30))
         else
!          SIMPLE: d = Vol/ap
           du(i,j,k) = vol/max(ap, 1.d-30)
         endif
       else if(dir == 2) then
         call porous_grad_cell_3pt(B%p(i,j-1,k), B%p(i,j,k), B%p(i,j+1,k), &
                            B%yc(i,j-1,k),B%yc(i,j,k),B%yc(i,j+1,k), gx)
         src = -gx*vol
!      Deferred correction: converged solution must satisfy the high-order
!      equation C_HO = Dform. Since C_HO = C_up + conv_src and the implicit
!      operator is the upwind one (ap*phi - sum a_nb*phi_nb = C_up - Dform),
!      the RHS source is -conv_src (see Ferziger & Peric, deferred correction).
         if(LS_Scheme > 1) src = src - conv_src(i,j,k)
         uold = B%U(3,i,j,k)
         unew = (aE*B%U(3,i+1,j,k)+aW*B%U(3,i-1,j,k)+aN*B%U(3,i,j+1,k)+aS*B%U(3,i,j-1,k) &
                +aT*B%U(3,i,j,k+1)+aB*B%U(3,i,j,k-1)+src)/max(ap,1.d-30)
         B%U(3,i,j,k) = uold + alpha*(unew-uold)
         apv(i,j,k) = ap
         sv_nb(i,j,k) = aE+aW+aN+aS+aT+aB
         if(LS_Algorithm == 2) then
           denom = ap - alpha*sv_nb(i,j,k)
           if(denom < 0.1d0*ap) denom = 0.1d0*ap
           dv(i,j,k) = min(vol*alpha/max(denom, 1.d-30), 3.d0*vol/max(ap,1.d-30))
         else
           dv(i,j,k) = vol/max(ap, 1.d-30)
         endif
       else
         call porous_grad_cell_3pt(B%p(i,j,k-1), B%p(i,j,k), B%p(i,j,k+1), &
                            B%zc(i,j,k-1),B%zc(i,j,k),B%zc(i,j,k+1), gx)
         src = -gx*vol
!      Deferred correction: converged solution must satisfy the high-order
!      equation C_HO = Dform. Since C_HO = C_up + conv_src and the implicit
!      operator is the upwind one (ap*phi - sum a_nb*phi_nb = C_up - Dform),
!      the RHS source is -conv_src (see Ferziger & Peric, deferred correction).
         if(LS_Scheme > 1) src = src - conv_src(i,j,k)
         uold = B%U(4,i,j,k)
         unew = (aE*B%U(4,i+1,j,k)+aW*B%U(4,i-1,j,k)+aN*B%U(4,i,j+1,k)+aS*B%U(4,i,j-1,k) &
                +aT*B%U(4,i,j,k+1)+aB*B%U(4,i,j,k-1)+src)/max(ap,1.d-30)
         B%U(4,i,j,k) = uold + alpha*(unew-uold)
         apw(i,j,k) = ap
         sw_nb(i,j,k) = aE+aW+aN+aS+aT+aB
         if(LS_Algorithm == 2) then
           denom = ap - alpha*sw_nb(i,j,k)
           if(denom < 0.1d0*ap) denom = 0.1d0*ap
           dw(i,j,k) = min(vol*alpha/max(denom, 1.d-30), 3.d0*vol/max(ap,1.d-30))
         else
           dw(i,j,k) = vol/max(ap, 1.d-30)
         endif
       endif
     enddo; enddo; enddo
   enddo
  end subroutine porous_momentum

!==============================================================================
! Set one ghost layer of the pressure-correction field pp.
!   - default: zero-gradient (Neumann), correct for walls, symmetry and
!     velocity/mass-flow inlets (pp does not couple through a Dirichlet p face)
!   - on faces where p itself is Dirichlet (pressure outlet BC_Outflow, or
!     pressure inlet LS_Inlet_Type=3) the correction must vanish: pp(ghost)=0
! Without this, pp(i-1)/pp(i,j-1)/pp(i,j,k-1) at the i-/j-/k- boundaries read
! out of bounds (pp was allocated 1:nx) and the p' solve carries an undefined
! (grid/layout dependent) boundary value -> p drifts and, with a pressure
! outlet, drives a spurious reverse flow.  Allocation is now (0:nx,0:ny,0:nz).
!==============================================================================
  subroutine porous_pp_ghost(nMesh, mBlock)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer:: i,j,k, nx,ny,nz, ksub, face_s, ib,ie,jb,je,kb,ke

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

!  default: Neumann (dp'/dn = 0) on all six ghost planes
   do k=0,nz; do j=0,ny
     pp(0,j,k)  = pp(1,j,k)
     pp(nx,j,k) = pp(nx-1,j,k)
   enddo; enddo
   do k=0,nz; do i=0,nx
     pp(i,0,k)  = pp(i,1,k)
     pp(i,ny,k) = pp(i,ny-1,k)
   enddo; enddo
   do j=0,ny; do i=0,nx
     pp(i,j,0)  = pp(i,j,1)
     pp(i,j,nz) = pp(i,j,nz-1)
   enddo; enddo

!  pressure-Dirichlet faces: p fixed on the face -> correction pp = 0 there
   do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     if(is_interface_bc(Bc%bc)) cycle
     if((Bc%bc == BC_Outflow .or. Bc%bc == BC_LS_Outlet) .or. &
       ((Bc%bc == BC_Inflow .or. Bc%bc == BC_LS_Inlet) .and. LS_Inlet_Type == 3)) then
       face_s = Bc%face
       ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke
       select case(face_s)
       case(1)   ! i- face: ghost cell i = ib-1
         do k=kb,ke-1; do j=jb,je-1
           pp(ib-1,j,k) = 0.d0
         enddo; enddo
       case(4)   ! i+ face: ghost cell i = ie
         do k=kb,ke-1; do j=jb,je-1
           pp(ie,j,k) = 0.d0
         enddo; enddo
       case(2)   ! j- face: ghost cell j = jb-1
         do k=kb,ke-1; do i=ib,ie-1
           pp(i,jb-1,k) = 0.d0
         enddo; enddo
       case(5)   ! j+ face: ghost cell j = je
         do k=kb,ke-1; do i=ib,ie-1
           pp(i,je,k) = 0.d0
         enddo; enddo
       case(3)   ! k- face: ghost cell k = kb-1
         do j=jb,je-1; do i=ib,ie-1
           pp(i,j,kb-1) = 0.d0
         enddo; enddo
       case(6)   ! k+ face: ghost cell k = ke
         do j=jb,je-1; do i=ib,ie-1
           pp(i,j,ke) = 0.d0
         enddo; enddo
       end select
     endif
   enddo
  end subroutine porous_pp_ghost

!==============================================================================
! Pressure-correction equation (Poisson) + velocity correction + pressure update.
! SOR (Successive Over-Relaxation) accelerates convergence on finer grids.
! Plain Gauss-Seidel spectral radius ~1-pi^2/(2N^2) is too slow for N>=80.
!==============================================================================
  subroutine porous_pressure_correction(nMesh, mBlock, res_p)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   real(PRE_EC):: res_p
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz, iter
    real(PRE_EC):: rho, ap, aE,aW,aN,aS,aT,aB, src, pnew, gx
   real(PRE_EC):: Fe,Fw,Fn,Fs,Ft,Fb, dxe,dxw,dyn,dys,dzt,dzb
   real(PRE_EC):: dbe,dbw,dbn,dbs,dbt,dbb
   real(PRE_EC):: omega, resid_max, pchange
   integer:: inner_max
   real(PRE_EC),parameter:: INNER_TOL=1.d-12

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   rho=LS_rho

!  Use plain Gauss-Seidel (omega=1.0) with enough iterations for 80x80 grid.
!  SOR (omega>1) causes divergence on this non-standard Poisson system.
!  GS spectral radius ~cos^2(pi/N) ~0.998 for N=80, needs ~400+ iters
!  to match 40x40 convergence level (200 iters at rho~0.997).
   omega = 1.0d0

!  Mass-conservation residual of the predicted velocity field = SIMPLE outer
!  convergence monitor.  Do NOT use max|pp| for this: |pp| stays O(1) even at
!  convergence (it is a pressure increment, not an imbalance).
   res_p = 0.d0
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
     res_p = max(res_p, abs(Fi(i+1,j,k)-Fi(i,j,k) &
                           +Fj(i,j+1,k)-Fj(i,j,k) &
                           +Fk(i,j,k+1)-Fk(i,j,k)))
   enddo; enddo; enddo

!  p' is solved fresh each outer sweep: pp holds only the current correction,
!  and after p/u are updated below it must NOT carry over to the next sweep
!  (otherwise p accumulates the same correction repeatedly).
   pp = 0.d0
   inner_max = 100
   if(LS_Algorithm == 2) inner_max = 500   ! SIMPLEC: p' benefits from a more
                                           ! converged Gauss-Seidel sweep
   do iter=1, inner_max
!    refresh the pp ghost layer from the latest interior values (Neumann
!    default, zero on pressure-Dirichlet faces)
     call porous_pp_ghost(nMesh, mBlock)
     resid_max = 0.d0
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
!      continuity imbalance (mass source)
       Fe = Fi(i+1,j,k); Fw = Fi(i,j,k)
       Fn = Fj(i,j+1,k); Fs = Fj(i,j,k)
       Ft = Fk(i,j,k+1); Fb = Fk(i,j,k)
       src = -(Fe-Fw+Fn-Fs+Ft-Fb)

!      pressure-correction coefficients (from Rhie-Chow face velocity)
       dxe = max(B%xc(i+1,j,k)-B%xc(i,j,k), 1.d-30)
       dxw = max(B%xc(i,j,k)-B%xc(i-1,j,k), 1.d-30)
       dyn = max(B%yc(i,j+1,k)-B%yc(i,j,k), 1.d-30)
       dys = max(B%yc(i,j,k)-B%yc(i,j-1,k), 1.d-30)
       dzt = max(B%zc(i,j,k+1)-B%zc(i,j,k), 1.d-30)
       dzb = max(B%zc(i,j,k)-B%zc(i,j,k-1), 1.d-30)
       dbe = 0.5d0*(du(i,j,k)+du(i+1,j,k))
       dbw = 0.5d0*(du(i,j,k)+du(i-1,j,k))
       dbn = 0.5d0*(dv(i,j,k)+dv(i,j+1,k))
       dbs = 0.5d0*(dv(i,j,k)+dv(i,j-1,k))
       dbt = 0.5d0*(dw(i,j,k)+dw(i,j,k+1))
       dbb = 0.5d0*(dw(i,j,k)+dw(i,j,k-1))
       aE = rho*B%Si(i+1,j,k)*dbe/dxe
       aW = rho*B%Si(i,j,k)*dbw/dxw
       aN = rho*B%Sj(i,j+1,k)*dbn/dyn
       aS = rho*B%Sj(i,j,k)*dbs/dys
       aT = rho*B%Sk(i,j,k+1)*dbt/dzt
       aB = rho*B%Sk(i,j,k)*dbb/dzb
       ap = aE+aW+aN+aS+aT+aB

       pnew = (aE*pp(i+1,j,k)+aW*pp(i-1,j,k)+aN*pp(i,j+1,k)+aS*pp(i,j-1,k) &
              +aT*pp(i,j,k+1)+aB*pp(i,j,k-1)+src)/max(ap,1.d-30)
!      SOR: over-relax the Gauss-Seidel update
       pchange = omega*(pnew - pp(i,j,k))
       pp(i,j,k) = pp(i,j,k) + pchange
       resid_max = max(resid_max, abs(pchange))
     enddo; enddo; enddo
!    early termination: solver converged when change is negligible
     if(resid_max < INNER_TOL .and. iter > 5) exit
   enddo

!  velocity correction + pressure update. Use the same 3-point cell gradient
!  for the pressure-correction difference as in the momentum source, to keep
!  the two parts of the SIMPLE split consistent on non-uniform meshes.
   call porous_pp_ghost(nMesh, mBlock)
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
     call porous_grad_cell_3pt(pp(i-1,j,k), pp(i,j,k), pp(i+1,j,k), &
                        B%xc(i-1,j,k),B%xc(i,j,k),B%xc(i+1,j,k), gx)
     B%U(2,i,j,k) = B%U(2,i,j,k) - du(i,j,k)*gx
     call porous_grad_cell_3pt(pp(i,j-1,k), pp(i,j,k), pp(i,j+1,k), &
                        B%yc(i,j-1,k),B%yc(i,j,k),B%yc(i,j+1,k), gx)
     B%U(3,i,j,k) = B%U(3,i,j,k) - dv(i,j,k)*gx
     call porous_grad_cell_3pt(pp(i,j,k-1), pp(i,j,k), pp(i,j,k+1), &
                        B%zc(i,j,k-1),B%zc(i,j,k),B%zc(i,j,k+1), gx)
     B%U(4,i,j,k) = B%U(4,i,j,k) - dw(i,j,k)*gx
!    Pressure update.  Pure SIMPLEC (d-factor above already accounts for the
!    neglected neighbor corrections) uses alpha_p = 1.0; a smaller alpha_p
!    behaves as SIMPLEC-d with under-relaxed pressure and is more robust when
!    the p' Gauss-Seidel sweep count is limited.
     B%p(i,j,k) = B%p(i,j,k) + LS_alpha_p*pp(i,j,k)
   enddo; enddo; enddo
  end subroutine porous_pressure_correction

!==============================================================================
! Energy (temperature) equation, first-order upwind + central diffusion.
!==============================================================================
  subroutine porous_energy(nMesh, mBlock)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz, iter
    real(PRE_EC):: rho, kcond, cph, alpha, ap, aE,aW,aN,aS,aT,aB, src, apstab, eps_p
    real(PRE_EC):: De,Dwe,Dn,Ds,Dt,Db, Fe,Fw,Fn,Fs,Ft,Fb, unew, uold, vol
   integer,parameter:: INNER=3

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   rho=LS_rho; kcond=LS_k; cph=LS_Cp; alpha=LS_alpha_T
    eps_p = B%porous_eps
    kcond = LS_k*eps_p   ! effective fluid-phase conductivity

!  high-order convection correction (deferred)
   if(LS_Scheme > 1) call porous_conv_high(nMesh, mBlock, 5)
   do iter=1, INNER
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       vol = B%Vol(i,j,k)
       De = kcond*B%Si(i+1,j,k)/max(B%xc(i+1,j,k)-B%xc(i,j,k), 1.d-30)
       Dwe = kcond*B%Si(i,j,k)/max(B%xc(i,j,k)-B%xc(i-1,j,k), 1.d-30)
       Dn = kcond*B%Sj(i,j+1,k)/max(B%yc(i,j+1,k)-B%yc(i,j,k), 1.d-30)
       Ds = kcond*B%Sj(i,j,k)/max(B%yc(i,j,k)-B%yc(i,j-1,k), 1.d-30)
       Dt = kcond*B%Sk(i,j,k+1)/max(B%zc(i,j,k+1)-B%zc(i,j,k), 1.d-30)
       Db = kcond*B%Sk(i,j,k)/max(B%zc(i,j,k)-B%zc(i,j,k-1), 1.d-30)

       Fe = cph*Fi(i+1,j,k); Fw = cph*Fi(i,j,k)
       Fn = cph*Fj(i,j+1,k); Fs = cph*Fj(i,j,k)
       Ft = cph*Fk(i,j,k+1); Fb = cph*Fk(i,j,k)

       aE = De + max(-Fe, 0.d0)
      aW = Dwe + max( Fw, 0.d0)
      aN = Dn + max(-Fn, 0.d0)
      aS = Ds + max( Fs, 0.d0)
      aT = Dt + max(-Ft, 0.d0)
      aB = Db + max( Fb, 0.d0)
!      Advection (skew-symmetric) form: only inflow enters the matrix, the
!      outflow term cancels with the T_P part. The net mass-flux imbalance
!      (Fe-Fw+...) must NOT appear here: during SIMPLE transients the fluxes
!      are not divergence-free and a divergence-form net term creates spurious
!      sources/sinks that drive T away from bounded values (Inf/NaN with k=0).
!      With this form the update is a convex combination of neighbor T and
!      T_old: bounded by construction, and T-uniform fields are preserved
!      exactly. At convergence (divergence-free fluxes) it reduces to the
!      standard steady convection-diffusion discretization.
      ap = aE+aW+aN+aS+aT+aB
      src = 0.d0
!      Deferred correction: outward-flux conv_src, applied as -conv_src
!      (same convention as momentum; see notes at lines 688/702/716).
      if(LS_Scheme > 1) src = src - conv_src(i,j,k)
      uold = B%U(5,i,j,k)

!      pseudo-transient stabilization: pure convection (k=0) can give ap=0.
!      Pseudo time step dt = rho*Cp*vol/sum|F| (CFL~1) guarantees ap>0.
      apstab = abs(Fe)+abs(Fw)+abs(Fn)+abs(Fs)+abs(Ft)+abs(Fb)
      ap = ap + apstab
      src = src + apstab*uold

!      Inter-phase heat exchange: h_v*Vol*(Ts - Tf), implicit slope.
      if(B%porous_hv > 0.d0) then
        ap = ap + B%porous_hv*vol
        src = src + B%porous_hv*vol*B%Ts(i,j,k)
      endif

       if(ap > 1.d-20) then
         unew = (aE*B%U(5,i+1,j,k)+aW*B%U(5,i-1,j,k)+aN*B%U(5,i,j+1,k)+aS*B%U(5,i,j-1,k) &
                +aT*B%U(5,i,j,k+1)+aB*B%U(5,i,j,k-1)+src)/ap
       else
         unew = uold
       endif
       B%U(5,i,j,k) = uold + alpha*(unew-uold)
     enddo; enddo; enddo
   enddo
  end subroutine porous_energy
!==============================================================================
! Solid-frame (skeleton) energy equation of the porous block (LTNE model).
!   (1-eps)*k_s * Laplacian(Ts) + h_v*(Tf - Ts) = 0
! First-order implicit Gauss-Seidel sweeps, no convection; under-relaxed with
! Porous_alpha_Ts.  Ghost cells of Ts are refreshed by porous_boundary_block.
!==============================================================================
  subroutine porous_energy_Ts(nMesh, mBlock)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k, nx,ny,nz, iter
   real(PRE_EC):: alpha, ap, aE,aW,aN,aS,aT,aB, src, vol
   real(PRE_EC):: De,Dwe,Dn,Ds,Dt,Db, kcond_s, Tsold, Tsnew, hv
   integer,parameter:: INNER=5

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   alpha=Porous_alpha_Ts
   kcond_s=(1.d0-B%porous_eps)*B%solid_k
   hv=B%porous_hv
   do iter=1, INNER
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       vol = B%Vol(i,j,k)
       De  = kcond_s*B%Si(i+1,j,k)/max(B%xc(i+1,j,k)-B%xc(i,j,k), 1.d-30)
       Dwe = kcond_s*B%Si(i,j,k)/max(B%xc(i,j,k)-B%xc(i-1,j,k), 1.d-30)
       Dn  = kcond_s*B%Sj(i,j+1,k)/max(B%yc(i,j+1,k)-B%yc(i,j,k), 1.d-30)
       Ds  = kcond_s*B%Sj(i,j,k)/max(B%yc(i,j,k)-B%yc(i,j-1,k), 1.d-30)
       Dt  = kcond_s*B%Sk(i,j,k+1)/max(B%zc(i,j,k+1)-B%zc(i,j,k), 1.d-30)
       Db  = kcond_s*B%Sk(i,j,k)/max(B%zc(i,j,k)-B%zc(i,j,k-1), 1.d-30)
       aE=De; aW=Dwe; aN=Dn; aS=Ds; aT=Dt; aB=Db
       ap  = aE+aW+aN+aS+aT+aB
       src = 0.d0
       if(hv > 0.d0) then
!        h_v*(Tf - Ts): implicit -h_v*Vol on the diagonal, +h_v*Vol*Tf on the RHS
         ap  = ap  + hv*vol
         src = src + hv*vol*B%U(5,i,j,k)
       endif
       Tsold = B%Ts(i,j,k)
       if(ap > 1.d-20) then
         Tsnew = (aE*B%Ts(i+1,j,k)+aW*B%Ts(i-1,j,k)+aN*B%Ts(i,j+1,k) &
               + aS*B%Ts(i,j-1,k)+aT*B%Ts(i,j,k+1)+aB*B%Ts(i,j,k-1)+src)/ap
       else
         Tsnew = Tsold
       endif
       B%Ts(i,j,k) = Tsold + alpha*(Tsnew-Tsold)
     enddo; enddo; enddo
   enddo
  end subroutine porous_energy_Ts

!==============================================================================
! Low-speed flow solver: SIMPLE algorithm (steady, incompressible)
!==============================================================================
  subroutine porous_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_Var
   use const_var
   use porous_work
   implicit none
   integer:: nMesh, mBlock
   real(PRE_EC):: Sfac, Sfac1
   Type (Block_TYPE),pointer:: B
   integer:: nx,ny,nz, iter, i,j,k
   real(PRE_EC):: res_p, res0

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

!  (re)allocate work arrays if size changed
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
   Fi=0.d0; Fj=0.d0; Fk=0.d0; pp=0.d0
   apu=1.d0; apv=1.d0; apw=1.d0
   su_nb=0.d0; sv_nb=0.d0; sw_nb=0.d0
   du=0.d0; dv=0.d0; dw=0.d0

   res0 = 0.d0
   do iter=1, Porous_Max_Iter
!    boundary conditions (ghost cells + boundary face fluxes)
     call porous_boundary_block(nMesh, mBlock)

!    momentum prediction
     call porous_momentum(nMesh, mBlock, 1)
     call porous_momentum(nMesh, mBlock, 2)
     call porous_momentum(nMesh, mBlock, 3)

!    face fluxes (interior Rhie-Chow; boundary set by BC)
     call porous_face_flux(nMesh, mBlock)

!    pressure correction + velocity/pressure update
     call porous_pressure_correction(nMesh, mBlock, res_p)

!    energy (temperature)
     call porous_energy(nMesh, mBlock)
!    solid-frame energy (LTNE)
      call porous_energy_Ts(nMesh, mBlock)

!    catch genuine (10^100) spikes from SIMPLE transients; wide clipping
!    leaves the well-tuned solution untouched but prevents NaN cascade.
     call porous_clip_fields(nMesh, mBlock)

!    DEBUG: locate NaN source
     do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
       if(B%U(5,i,j,k) /= B%U(5,i,j,k)) then
         print*, 'DEBUG T NaN at iter', iter, ' cell', i,j,k, &
           ' T(i-1)', B%U(5,i-1,j,k), ' T(i+1)', B%U(5,i+1,j,k), &
           ' T(j-1)', B%U(5,i,j-1,k), ' T(j+1)', B%U(5,i,j+1,k), &
           ' p', B%p(i,j,k), ' Fi(i,j)', Fi(i,j,k)
         stop
       endif
     enddo; enddo; enddo

!    DEBUG: temperature convergence trend
     if(my_id == 0 .and. mod(iter,50) == 0) then
       print*, 'DEBUG T iter', iter, 'min', minval(B%U(5,1:nx-1,1:ny-1,1:nz-1)), &
               'max', maxval(B%U(5,1:nx-1,1:ny-1,1:nz-1))
     endif

     if(iter == 1) res0 = res_p
     if(my_id == 0 .and. iter <= 3) then
       print*, '  DBG iter', iter, ' maxU2=', maxval(abs(B%U(2,1:nx,1:ny,1:nz))), &
               ' maxFi=', maxval(abs(Fi(1:nx,1:ny,1:nz))), &
               ' maxFj=', maxval(abs(Fj(1:nx,1:ny,1:nz))), &
               ' Fi(21,1,1)=', Fi(21,1,1), ' Fi(21,21,1)=', Fi(21,21,1), &
               ' Fi(41,21,1)=', Fi(41,21,1)
       print*, '  DBG div max=', maxval(abs( Fi(2:nx,1:ny,1:nz)-Fi(1:nx-1,1:ny,1:nz) &
             + Fj(1:nx,2:ny,1:nz)-Fj(1:nx,1:ny-1,1:nz) )), &
               ' res_p=', res_p
     endif
     if(my_id == 0 .and. (iter <= 5 .or. mod(iter,100) == 0)) then
       print*, "  SIMPLE iter", iter, " res_p=", res_p, &
               " u(1,1,1)=", B%U(2,1,1,1), " p(1,1,1)=", B%p(1,1,1)
     endif
!    minimum sweep count avoids premature exit from the quiescent start
!    (initial u=0,p=0 gives an identically-zero mass residual on sweep 1)
     if(iter > 10 .and. res_p < Porous_Tol) exit
   enddo

   if(my_id == 0) then
     print*, "Porous solver Block", B%Block_no, " SIMPLE iterations:", iter, &
             " res_p:", res_p
   endif
  end subroutine porous_solver_one_block

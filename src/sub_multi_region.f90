!----------------------------------------------------------------------
! Multi-region solver stubs
! Supports: fluid (0), solid (1), low-speed (2), porous (3)
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Dispatcher: select solver based on block type
!----------------------------------------------------------------------
  subroutine solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_Var
   implicit none
   integer:: nMesh, mBlock
   real(PRE_EC):: Sfac, Sfac1
   Type (Block_TYPE),pointer:: B
   B=>Mesh(nMesh)%Block(mBlock)
   select case(B%Block_type)
   case(BLOCK_FLUID)
     call Residual_one_block(nMesh, mBlock, Sfac, Sfac1)
     call Uupdate_one_block(nMesh, mBlock)
   case(BLOCK_SOLID)
     call solid_solver_one_block(nMesh, mBlock)
   case(BLOCK_LOWSPEED)
     if(LS_Algorithm .eq. 3) then
       call lowspeed_ac_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
     else
       call lowspeed_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
     endif
   case(BLOCK_POROUS)
     if(LS_Algorithm .eq. 3) then
       call porous_ac_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
     else
       call porous_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
     endif
   case default
     call Residual_one_block(nMesh, mBlock, Sfac, Sfac1)
   end select
  end subroutine solver_one_block

!----------------------------------------------------------------------
! Helper: set ghost cells for a solid block from physical BCs
! Handles bc=2 (wall: Tw/Qw), bc=3 (symmetry: adiabatic)
! bc=11 (fluid-solid), bc=12 (solid-solid), bc=13 (solid-porous)
! are skipped here and handled by coupling routines + buffer exchange.
!----------------------------------------------------------------------
  subroutine set_solid_ghost_BC(nMesh, mBlock)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer:: i,j,k,n,ii,ksub,face_s, i1,j1,k1, i2,j2,k2
   integer:: ib,ie,jb,je,kb,ke, nx,ny,nz
   real(PRE_EC):: Tw_val, Qw_val, dx_b, dy_b, dz_b
   logical:: found_bc

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

!  Loop over all subfaces
   do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     if(Bc%bc .lt. 0) cycle   ! skip internal interfaces (handled by buffer exchange)

!    Skip interface connections (handled by coupling routines).
!    Check bc_msg2: if it has bc<0 (interface connection), skip this subface.
     if(associated(B%bc_msg2)) then
       if(is_interface_bc(B%bc_msg2(ksub)%bc)) cycle
     endif

     face_s = Bc%face
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke

!    Debug: print subface info (first call only)
!    if(my_id == 0) print*, '  DBG: ksub=', ksub, ' bc=', Bc%bc, ' face=', face_s, ' ib=', ib, ' ie=', ie

!    Look up thermal BC type from solid_bc data
!    solid_bc_face_no stores the face number (1=i-,2=j-,3=k-,4=i+,5=j+,6=k+),
!    compare against Bc%face (the face number of the current subface).
     found_bc = .false.
     Tw_val = 0.d0; Qw_val = 0.d0
     do ii=1, B%solid_bc_nface
       if(B%solid_bc_face_no(ii) == Bc%face) then
         found_bc = .true.
         if(B%solid_bc_Tw(ii) > 0.d0) then
!          Isothermal: Tw > 0 (K) -> use dimensional value (Ts is dimensional)
           Tw_val = B%solid_bc_Tw(ii)
         else
!          Heat flux: Tw < 0, use Qw (W/m2) -> dimensional
           Qw_val = B%solid_bc_Qw(ii)
         endif
         exit
       endif
     enddo

!    Determine BC type and apply to ghost cells
!    Note: bc3d.inp uses node indices. Cell-center arrays (Ts) use cell indices.
!    For face 1 (i-): node ib=1 �1�71�1�771�1�71�1�777�1�71�1�771�1�71�1�777 first cell center at i=ib, ghost cells at i=ib-1..ib-LAP
!    For face 4 (i+): node ie=nx �1�71�1�771�1�71�1�777�1�71�1�771�1�71�1�777 last cell center at i=ie-1, ghost cells at i=ie..ie+LAP-1
     select case(face_s)
     case(1)  ! i- face: ghost i=ib-1..ib-LAP, interior i=ib..ib+LAP-1
       do n=1, LAP
         i1 = ib - n                ! ghost cell: 0, -1, -2, -3
         i2 = ib + n - 1            ! interior cell: 1, 2, 3, 4
         do k=kb,ke; do j=jb,je
           select case(Bc%bc)
           case(2)  ! Wall
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i1,j,k) = 2.d0*Tw_val - B%Ts(i2,j,k)
             elseif(found_bc .and. Qw_val /= 0.d0) then
!              Use Euclidean cell-center distance for first ghost layer (FVM accurate).
!              Single-component abs(xc) is WRONG for non-Cartesian grids (e.g. O-grid
!              annulus where radial direction is y-only at theta=pi/2).
              if(n == 1) then
                dx_b = sqrt( (B%xc(i2,j,k)-B%xc(i1,j,k))**2 &
                            +(B%yc(i2,j,k)-B%yc(i1,j,k))**2 &
                            +(B%zc(i2,j,k)-B%zc(i1,j,k))**2 ) * Lscale
              else
                dx_b = sqrt( (B%x(ib+1,j,k)-B%x(ib,j,k))**2 &
                            +(B%y(ib+1,j,k)-B%y(ib,j,k))**2 &
                            +(B%z(ib+1,j,k)-B%z(ib,j,k))**2 ) * Lscale * n
              endif
              B%Ts(i1,j,k) = B%Ts(i2,j,k) + Qw_val * dx_b / B%solid_k
            else
              B%Ts(i1,j,k) = B%Ts(i2,j,k)
            endif
           case(3)  ! Symmetry -> adiabatic
             B%Ts(i1,j,k) = B%Ts(i2,j,k)
           case default
             B%Ts(i1,j,k) = B%Ts(i2,j,k)
           end select
         enddo; enddo
       enddo
     case(4)  ! i+ face: ghost i=ie..ie+LAP-1, interior i=ie-1..ie-LAP
       do n=1, LAP
         i1 = ie + n - 1            ! ghost cell: 51, 52, 53, 54
         i2 = ie - n                ! interior cell: 50, 49, 48, 47
         do k=kb,ke; do j=jb,je
           select case(Bc%bc)
           case(2)
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i1,j,k) = 2.d0*Tw_val - B%Ts(i2,j,k)
             elseif(found_bc .and. Qw_val /= 0.d0) then
              if(n == 1) then
                dx_b = sqrt( (B%xc(i1,j,k)-B%xc(i2,j,k))**2 &
                            +(B%yc(i1,j,k)-B%yc(i2,j,k))**2 &
                            +(B%zc(i1,j,k)-B%zc(i2,j,k))**2 ) * Lscale
              else
                dx_b = sqrt( (B%x(ie,j,k)-B%x(ie-1,j,k))**2 &
                            +(B%y(ie,j,k)-B%y(ie-1,j,k))**2 &
                            +(B%z(ie,j,k)-B%z(ie-1,j,k))**2 ) * Lscale * n
              endif
              B%Ts(i1,j,k) = B%Ts(i2,j,k) + Qw_val * dx_b / B%solid_k
             else
               B%Ts(i1,j,k) = B%Ts(i2,j,k)
             endif
           case(3)
             B%Ts(i1,j,k) = B%Ts(i2,j,k)
           case default
             B%Ts(i1,j,k) = B%Ts(i2,j,k)
           end select
         enddo; enddo
       enddo
     case(2)  ! j- face: ghost j=jb-1..jb-LAP, interior j=jb..jb+LAP-1
       do n=1, LAP
         j1 = jb - n
         j2 = jb + n - 1
         do k=kb,ke; do i=ib,ie
           select case(Bc%bc)
           case(2)
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i,j1,k) = 2.d0*Tw_val - B%Ts(i,j2,k)
             elseif(found_bc .and. Qw_val /= 0.d0) then
              if(n == 1) then
                dy_b = sqrt( (B%xc(i,j2,k)-B%xc(i,j1,k))**2 &
                            +(B%yc(i,j2,k)-B%yc(i,j1,k))**2 &
                            +(B%zc(i,j2,k)-B%zc(i,j1,k))**2 ) * Lscale
              else
                dy_b = sqrt( (B%x(i,jb+1,k)-B%x(i,jb,k))**2 &
                            +(B%y(i,jb+1,k)-B%y(i,jb,k))**2 &
                            +(B%z(i,jb+1,k)-B%z(i,jb,k))**2 ) * Lscale * n
              endif
              B%Ts(i,j1,k) = B%Ts(i,j2,k) + Qw_val * dy_b / B%solid_k
             else
               B%Ts(i,j1,k) = B%Ts(i,j2,k)
             endif
           case(3)
             B%Ts(i,j1,k) = B%Ts(i,j2,k)
           case default
             B%Ts(i,j1,k) = B%Ts(i,j2,k)
           end select
         enddo; enddo
       enddo
     case(5)  ! j+ face: ghost j=je..je+LAP-1, interior j=je-1..je-LAP
       do n=1, LAP
         j1 = je + n - 1
         j2 = je - n
         do k=kb,ke; do i=ib,ie
           select case(Bc%bc)
           case(2)
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i,j1,k) = 2.d0*Tw_val - B%Ts(i,j2,k)
             elseif(found_bc .and. Qw_val /= 0.d0) then
              if(n == 1) then
                dy_b = sqrt( (B%xc(i,j1,k)-B%xc(i,j2,k))**2 &
                            +(B%yc(i,j1,k)-B%yc(i,j2,k))**2 &
                            +(B%zc(i,j1,k)-B%zc(i,j2,k))**2 ) * Lscale
              else
                dy_b = sqrt( (B%x(i,je,k)-B%x(i,je-1,k))**2 &
                            +(B%y(i,je,k)-B%y(i,je-1,k))**2 &
                            +(B%z(i,je,k)-B%z(i,je-1,k))**2 ) * Lscale * n
              endif
              B%Ts(i,j1,k) = B%Ts(i,j2,k) + Qw_val * dy_b / B%solid_k
             else
               B%Ts(i,j1,k) = B%Ts(i,j2,k)
             endif
           case(3)
             B%Ts(i,j1,k) = B%Ts(i,j2,k)
           case default
             B%Ts(i,j1,k) = B%Ts(i,j2,k)
           end select
         enddo; enddo
       enddo
     case(3)  ! k- face: ghost k=kb-1..kb-LAP, interior k=kb..kb+LAP-1
       do n=1, LAP
         k1 = kb - n
         k2 = kb + n - 1
         do j=jb,je; do i=ib,ie
           select case(Bc%bc)
           case(2)
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i,j,k1) = 2.d0*Tw_val - B%Ts(i,j,k2)
             elseif(found_bc .and. Qw_val /= 0.d0) then
              if(n == 1) then
                dz_b = sqrt( (B%xc(i,j,k2)-B%xc(i,j,k1))**2 &
                            +(B%yc(i,j,k2)-B%yc(i,j,k1))**2 &
                            +(B%zc(i,j,k2)-B%zc(i,j,k1))**2 ) * Lscale
              else
                dz_b = sqrt( (B%x(i,j,kb+1)-B%x(i,j,kb))**2 &
                            +(B%y(i,j,kb+1)-B%y(i,j,kb))**2 &
                            +(B%z(i,j,kb+1)-B%z(i,j,kb))**2 ) * Lscale * n
              endif
              B%Ts(i,j,k1) = B%Ts(i,j,k2) + Qw_val * dz_b / B%solid_k
             else
               B%Ts(i,j,k1) = B%Ts(i,j,k2)
             endif
           case(3)
             B%Ts(i,j,k1) = B%Ts(i,j,k2)
           case default
             B%Ts(i,j,k1) = B%Ts(i,j,k2)
           end select
         enddo; enddo
       enddo
     case(6)  ! k+ face: ghost k=ke..ke+LAP-1, interior k=ke-1..ke-LAP
       do n=1, LAP
         k1 = ke + n - 1
         k2 = ke - n
         do j=jb,je; do i=ib,ie
           select case(Bc%bc)
           case(2)
             if(found_bc .and. Tw_val > 0.d0) then
               B%Ts(i,j,k1) = 2.d0*Tw_val - B%Ts(i,j,k2)
             elseif(found_bc .and. Qw_val /= 0.d0) then
              if(n == 1) then
                dz_b = sqrt( (B%xc(i,j,k1)-B%xc(i,j,k2))**2 &
                            +(B%yc(i,j,k1)-B%yc(i,j,k2))**2 &
                            +(B%zc(i,j,k1)-B%zc(i,j,k2))**2 ) * Lscale
              else
                dz_b = sqrt( (B%x(i,j,ke)-B%x(i,j,ke-1))**2 &
                            +(B%y(i,j,ke)-B%y(i,j,ke-1))**2 &
                            +(B%z(i,j,ke)-B%z(i,j,ke-1))**2 ) * Lscale * n
              endif
              B%Ts(i,j,k1) = B%Ts(i,j,k2) + Qw_val * dz_b / B%solid_k
             else
               B%Ts(i,j,k1) = B%Ts(i,j,k2)
             endif
           case(3)
             B%Ts(i,j,k1) = B%Ts(i,j,k2)
           case default
             B%Ts(i,j,k1) = B%Ts(i,j,k2)
           end select
         enddo; enddo
       enddo
     end select
   enddo
  end subroutine set_solid_ghost_BC

!----------------------------------------------------------------------
! Solid heat conduction solver
! Gauss-Seidel implicit iteration for the heat equation:
!   Steady:   Laplacian(T) = 0
!   Unsteady: rho*Cp*dT/dt = k*Laplacian(T)
! Implicit Euler: T_new = (T_old + Fo*sum_neighbors) / (1 + 6*Fo)
!   where Fo = k*dt/(rho*Cp*dx^2)  (Fourier number)
!----------------------------------------------------------------------
  subroutine solid_solver_one_block(nMesh, mBlock)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i,j,k,ii, nx,ny,nz
   real(PRE_EC):: T_old, res, T_new
   real(PRE_EC):: dx_p, dx_m, dy_p, dy_m, dz_p, dz_m
   real(PRE_EC):: ax_p, ax_m, ay_p, ay_m, az_p, az_m
   real(PRE_EC):: coef_diff, coef_tr, T_neighbor_sum
   logical:: is_unsteady
!  NOTE: the solid Gauss-Seidel controls (SOR factor / sweep limits / tolerance)
!  are now read from control.ec namelist "$solid_ec" (globals Solid_GS_Omega,
!  Solid_Max_Iter, Solid_Min_Iter, Solid_Tol).  Their default values are exactly
!  the constants that used to be hard-coded here, so results are unchanged.
   real(PRE_EC),parameter:: EPS_DIST = 1.d-15   ! min distance to avoid div-by-zero

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

!  Determine if unsteady: RK3 or Dual-time stepping
   is_unsteady = (Time_Method == Time_RK3 .or. Time_Method == Time_Dual_LU_SGS)

!  1. Set physical BC ghost cells (wall, symmetry)
   call set_solid_ghost_BC(nMesh, mBlock)

!  2. Gauss-Seidel iteration with FVM non-uniform grid support
!  For cell (i,j,k), the steady-state heat flux balance:
!    k * sum[A_face * (T_neighbor - T_i) / d_face] = 0
!  where A_face = Si/Sj/Sk (face area), d_face = distance between cell centers.
!  Solving for T_i gives a distance-weighted average of neighbors.
   do ii = 1, Solid_Max_Iter
!    Re-set physical BC ghost cells each iteration
     call set_solid_ghost_BC(nMesh, mBlock)
!    Exchange Ts buffer for periodic/internal interfaces
     call update_Ts_buffer_onemesh(nMesh)
     call apply_solid_qw_override(nMesh, mBlock)   ! wall-flux split mode: re-impose interface q_w ghost

     res = 0.d0
     do k = 1, nz-1
     do j = 1, ny-1
     do i = 1, nx-1
       T_old = B%Ts(i,j,k)

!      Cell-center distances (FVM: Euclidean distance between adjacent cell centers)
!      Must use full 3D distance, not just x/y/z component, for non-Cartesian grids.
!      Apply non-orthogonal correction: project cell-center connection vector onto
!      face normal to get the effective normal distance used in Fourier's law.
!      This removes cross-diffusion error on skewed/non-orthogonal grids.
       block
         real(PRE_EC):: vx, vy, vz, vmag, cphi
         vx = B%xc(i+1,j,k)-B%xc(i,j,k)
         vy = B%yc(i+1,j,k)-B%yc(i,j,k)
         vz = B%zc(i+1,j,k)-B%zc(i,j,k)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%ni1(i+1,j,k) + vy*B%ni2(i+1,j,k) + vz*B%ni3(i+1,j,k)) / vmag
           dx_p = vmag * max(cphi, 0.1d0)   ! clamp cphi to avoid degenerate
         else
           dx_p = 0.d0
         endif
         vx = B%xc(i,j,k)-B%xc(i-1,j,k)
         vy = B%yc(i,j,k)-B%yc(i-1,j,k)
         vz = B%zc(i,j,k)-B%zc(i-1,j,k)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%ni1(i,j,k) + vy*B%ni2(i,j,k) + vz*B%ni3(i,j,k)) / vmag
           dx_m = vmag * max(cphi, 0.1d0)
         else
           dx_m = 0.d0
         endif
         vx = B%xc(i,j+1,k)-B%xc(i,j,k)
         vy = B%yc(i,j+1,k)-B%yc(i,j,k)
         vz = B%zc(i,j+1,k)-B%zc(i,j,k)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%nj1(i,j+1,k) + vy*B%nj2(i,j+1,k) + vz*B%nj3(i,j+1,k)) / vmag
           dy_p = vmag * max(cphi, 0.1d0)
         else
           dy_p = 0.d0
         endif
         vx = B%xc(i,j,k)-B%xc(i,j-1,k)
         vy = B%yc(i,j,k)-B%yc(i,j-1,k)
         vz = B%zc(i,j,k)-B%zc(i,j-1,k)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%nj1(i,j,k) + vy*B%nj2(i,j,k) + vz*B%nj3(i,j,k)) / vmag
           dy_m = vmag * max(cphi, 0.1d0)
         else
           dy_m = 0.d0
         endif
         vx = B%xc(i,j,k+1)-B%xc(i,j,k)
         vy = B%yc(i,j,k+1)-B%yc(i,j,k)
         vz = B%zc(i,j,k+1)-B%zc(i,j,k)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%nk1(i,j,k+1) + vy*B%nk2(i,j,k+1) + vz*B%nk3(i,j,k+1)) / vmag
           dz_p = vmag * max(cphi, 0.1d0)
         else
           dz_p = 0.d0
         endif
         vx = B%xc(i,j,k)-B%xc(i,j,k-1)
         vy = B%yc(i,j,k)-B%yc(i,j,k-1)
         vz = B%zc(i,j,k)-B%zc(i,j,k-1)
         vmag = sqrt(vx*vx + vy*vy + vz*vz)
         if(vmag > EPS_DIST) then
           cphi = abs(vx*B%nk1(i,j,k) + vy*B%nk2(i,j,k) + vz*B%nk3(i,j,k)) / vmag
           dz_m = vmag * max(cphi, 0.1d0)
         else
           dz_m = 0.d0
         endif
       end block

!      FVM coefficients: face_area / distance
!      Skip degenerate directions (e.g., 2D problems with 1 cell in k)
       if(abs(dx_p) > EPS_DIST) then
         ax_p = B%Si(i+1,j,k) / dx_p
       else
         ax_p = 0.d0
       endif
       if(abs(dx_m) > EPS_DIST) then
         ax_m = B%Si(i,j,k) / dx_m
       else
         ax_m = 0.d0
       endif
       if(abs(dy_p) > EPS_DIST) then
         ay_p = B%Sj(i,j+1,k) / dy_p
       else
         ay_p = 0.d0
       endif
       if(abs(dy_m) > EPS_DIST) then
         ay_m = B%Sj(i,j,k) / dy_m
       else
         ay_m = 0.d0
       endif
       if(abs(dz_p) > EPS_DIST) then
         az_p = B%Sk(i,j,k+1) / dz_p
       else
         az_p = 0.d0
       endif
       if(abs(dz_m) > EPS_DIST) then
         az_m = B%Sk(i,j,k) / dz_m
       else
         az_m = 0.d0
       endif

!      Sum of diffusion coefficients
       coef_diff = ax_p + ax_m + ay_p + ay_m + az_p + az_m

!      Weighted sum of neighbor temperatures
       T_neighbor_sum = ax_p*B%Ts(i+1,j,k) + ax_m*B%Ts(i-1,j,k) + &
                        ay_p*B%Ts(i,j+1,k) + ay_m*B%Ts(i,j-1,k) + &
                        az_p*B%Ts(i,j,k+1) + az_m*B%Ts(i,j,k-1)

       if(is_unsteady .and. dt_global > 0.d0) then
!        Transient term: rho*Cp*Vol*Lscale^2/dt
!        (Lscale converts mesh length to physical length)
         coef_tr = B%solid_rho * B%solid_Cp * B%Vol(i,j,k) * Lscale*Lscale / dt_global
!        Implicit Euler: (coef_tr * T_old + k * T_neighbor_sum) / (coef_tr + k * coef_diff)
         T_new = (coef_tr * T_old + B%solid_k * T_neighbor_sum) / &
                 (coef_tr + B%solid_k * coef_diff)
       else
!        Steady state: weighted average of neighbors
         T_new = T_neighbor_sum / (coef_diff + 1.d-30)
       endif

!      SOR: Ts_new = (1-omega)*Ts_old + omega*T_new
       B%Ts(i,j,k) = (1.d0-Solid_GS_Omega)*T_old + Solid_GS_Omega*T_new

       res = max(res, abs(B%Ts(i,j,k) - T_old))
     enddo; enddo; enddo

!    Debug: print progress for first few iterations
     if(ii <= 5 .and. my_id == 0) then
       print*, "  GS iter", ii, " res=", res, " Ts(1,1,1)=", B%Ts(1,1,1), &
               " Ts(nx-1,ny-1,1)=", B%Ts(nx-1,ny-1,1)
     endif

     if(mod(ii,100) == 0 .and. ii >= Solid_Min_Iter) then
       if(res < Solid_Tol) exit
     endif
   enddo

   if(my_id == 0 .or. B%Block_no == 3) then
     print*, "Solid solver Block", B%Block_no, " GS iterations:", ii, " final res:", res
     print*, "  DBG: Ts(1,1,1)=", B%Ts(1,1,1), " Ts(nx/2,ny/2,1)=", B%Ts(nx/2,ny/2,1), &
             " Ts(nx-1,ny-1,1)=", B%Ts(nx-1,ny-1,1)
   endif
  end subroutine solid_solver_one_block

!==============================================================================
! Solid GS interface heat-flux override (wall-flux CHT split mode).
! When stg_ow_block==mBlock, re-impose the interface Ts ghost from fp_qw each GS
! sweep so the steady conduction solve keeps the conjugate heat flux boundary
! (otherwise set_solid_ghost_BC would overwrite the interface ghost with the
! physical wall BC).  Supports j- (2) and j+ (5) solid faces.
!==============================================================================
  subroutine apply_solid_qw_override(nMesh, mBlock)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh, mBlock
   Type (Block_TYPE),pointer:: B
   integer:: i, k, n1, js, jg
   real(PRE_EC):: qw, dx

   if(stg_ow_block .le. 0) return
   if(stg_ow_block .ne. mBlock) return
   B => Mesh(nMesh)%Block(mBlock)
   if(stg_ow_face .eq. 2) then
     js = stg_ow_jb
   else if(stg_ow_face .eq. 5) then
     js = stg_ow_je - 1
   else
     return
   endif
   do k = stg_ow_kb, stg_ow_ke-1
   do i = stg_ow_ib, stg_ow_ie-1
     qw = fp_qw(i-stg_ow_ib+1, k-stg_ow_kb+1)
     do n1 = 0, LAP-1
       if(stg_ow_face .eq. 2) then
         jg = stg_ow_jb - 1 - n1
       else
         jg = stg_ow_je + n1
       endif
       dx = sqrt( (B%xc(i,js,k)-B%xc(i,jg,k))**2 &
                + (B%yc(i,js,k)-B%yc(i,jg,k))**2 &
                + (B%zc(i,js,k)-B%zc(i,jg,k))**2 ) * Lscale
       B%Ts(i,jg,k) = B%Ts(i,js,k) + qw*dx/max(B%solid_k,1.d-30)
     enddo
   enddo; enddo
  end subroutine apply_solid_qw_override


!----------------------------------------------------------------------
! Fluid-solid interface coupling
! For each bc=11 interface, compute interface temperature from
! heat flux balance:
!   T_i = (k_f*T_f/dx_f + k_s*T_s/dx_s) / (k_f/dx_f + k_s/dx_s)
! Then set ghost cells on both sides:
!   Solid ghost: T_ghost_s = 2*T_i - T_s_interior
!   Fluid ghost: isothermal wall at T_i (set via boundary condition)
!----------------------------------------------------------------------
  subroutine couple_fluid_solid_interfaces(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: mBlock, ksub, nb, mb, n, ii, face_s
   integer:: ib,ie,jb,je,kb,ke, i1,j1,k1, i2,j2,k2
   integer:: ib1,ie1,jb1,je1,kb1,ke1, face1
   real(PRE_EC):: T_f, T_s, k_f, k_s, dx_f, dx_s, T_i, T_ghost
   real(PRE_EC):: d1, uu1, v1, w1, p1, T1, d2, uu2, v2, w2, p2, T2
   integer:: i, j, k, NVAR1

   NVAR1 = Mesh(nMesh)%NVAR

   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%bc_msg2)) cycle

     do ksub=1, B%subface
       Bc2 => B%bc_msg2(ksub)
!      Skip physical boundaries (bc >= 0); only process interface connections (bc < 0)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1             ! Neighbor block number (global)
       if(nb <= 0) cycle
       if(B_proc(nb) .ne. my_id) then
!        Neighbour on another process: cross-process conjugate (CHT) coupling.
!        Both processes owning the pair call cht_mpi_fluid_solid_face, each with
!        its own block and its own interface entry Bc2.
         if((B%Block_type == BLOCK_SOLID .and. &
             (Block_Type_List(nb) == BLOCK_FLUID .or. Block_Type_List(nb) == BLOCK_LOWSPEED)) .or. &
            ((B%Block_type == BLOCK_FLUID .or. B%Block_type == BLOCK_LOWSPEED) .and. &
             Block_Type_List(nb) == BLOCK_SOLID)) then
           call cht_mpi_fluid_solid_face(nMesh, B, nb, Bc2)
         endif
         cycle
       endif
       mb = B_n(nb)
       Bn => Mesh(nMesh)%Block(mb)

!      Same-class fluid-fluid interfaces are handled by the buffer exchange
!      mechanism, not by this routine: skip them here.
       if(B%Block_type == BLOCK_FLUID .and. Bn%Block_type == BLOCK_FLUID) cycle
       if(B%Block_type == BLOCK_LOWSPEED .and. Bn%Block_type == BLOCK_LOWSPEED) cycle

!      Determine if this is a fluid-solid interface (compressible BLOCK_FLUID
!      or incompressible BLOCK_LOWSPEED against BLOCK_SOLID).  Each pair is
!      processed ONCE from the solid side (B = solid block) so that the
!      interface entry Bc2 (face/face1/ranges) has the solid face as "face".
       if(B%Block_type == BLOCK_SOLID .and. &
          .not.(Bn%Block_type == BLOCK_FLUID .or. Bn%Block_type == BLOCK_LOWSPEED)) cycle
       if(Bn%Block_type == BLOCK_SOLID) then
         if(.not.(B%Block_type == BLOCK_FLUID .or. B%Block_type == BLOCK_LOWSPEED)) cycle
         cycle   ! fluid-side entry: handled once from the solid block above
       endif

!      High-speed (compressible) <-> low-speed fluid pair: handled once from
!      the compressible side (B = FLUID, Bn = LOWSPEED)
       if(B%Block_type == BLOCK_FLUID .and. Bn%Block_type == BLOCK_LOWSPEED) then
         call couple_highlow_fluid_face(nMesh, B, Bn, Bc2)
         cycle
       else if(B%Block_type == BLOCK_LOWSPEED .and. Bn%Block_type == BLOCK_FLUID) then
         cycle   ! handled from the compressible side above
       else if(B%Block_type == BLOCK_FLUID .and. Bn%Block_type == BLOCK_POROUS) then
!        Interface 19: coolant velocity inlet (LS_Inlet_Type=1) activates the
!        transpiration (blowing) coupling; otherwise the phase-A shear-type
!        ghost exchange is used (G=0 reference).
         if(LS_Inlet_Type == 1) then
           call couple_compressible_porous_blowing_face(nMesh, B, Bn, Bc2)
         else
           call couple_highlow_fluid_face(nMesh, B, Bn, Bc2)
         endif
         cycle
       else if(B%Block_type == BLOCK_POROUS .and. Bn%Block_type == BLOCK_FLUID) then
         cycle   ! handled from the compressible side above
       endif

       face_s = Bc2%face
       ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
       face1 = Bc2%face1
       ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1

!      Here B is the solid block, Bn the fluid block.
       if(Bn%Block_type == BLOCK_FLUID) then
!        Compressible fluid: per-face-cell interface temperature (unit
!        conversion to physical K done inside)
         call couple_compressible_fluid_solid_face(nMesh, B, Bn, Bc2)
       else
!        Incompressible (low-speed) fluid: T is stored directly in U(5) [K]
!        and the thermal conductivity is the low-speed LS_k (same units as
!        the solid k).  The fluid ghost is a no-slip isothermal wall at T_i.
         call couple_lowspeed_fluid_solid_face(nMesh, B, Bn, Bc2)
       endif
     enddo
   enddo

   contains

!    Compute interface temperature from heat flux balance
     subroutine compute_interface_T(nMesh, Bs, Bf, ksub, Bc2, T_i)
      implicit none
      Type (Block_TYPE),pointer:: Bs, Bf
      integer:: nMesh, ksub
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      real(PRE_EC):: T_i
      real(PRE_EC):: T_s, T_f, k_s, k_f
      real(PRE_EC):: dx_s, dx_f
      integer:: is, js, ks, iflu, jflu, kflu, face_s, n

!     Solid side: interior cell adjacent to interface
      face_s = Bc2%face  ! face of solid block at interface
      n = 1              ! innermost ghost layer
      select case(face_s)
      case(1); is = Bc2%ib + LAP; js = Bc2%jb; ks = Bc2%kb
      case(4); is = Bc2%ie - LAP; js = Bc2%je; ks = Bc2%ke
      case(2); js = Bc2%jb + LAP; is = Bc2%ib; ks = Bc2%kb
      case(5); js = Bc2%je - LAP; is = Bc2%ie; ks = Bc2%ke
      case(3); ks = Bc2%kb + LAP; is = Bc2%ib; js = Bc2%jb
      case(6); ks = Bc2%ke - LAP; is = Bc2%ie; js = Bc2%je
      end select
      T_s = Bs%Ts(is,js,ks)
      k_s = Bs%solid_k

!     Fluid side: interior cell adjacent to interface (from U)
      face_s = Bc2%face1  ! face of fluid block at interface
      select case(face_s)
      case(1); iflu = Bc2%ib1 + LAP; jflu = Bc2%jb1; kflu = Bc2%kb1
      case(4); iflu = Bc2%ie1 - LAP; jflu = Bc2%je1; kflu = Bc2%ke1
      case(2); jflu = Bc2%jb1 + LAP; iflu = Bc2%ib1; kflu = Bc2%kb1
      case(5); jflu = Bc2%je1 - LAP; iflu = Bc2%ie1; kflu = Bc2%ke1
      case(3); kflu = Bc2%kb1 + LAP; iflu = Bc2%ib1; jflu = Bc2%jb1
      case(6); kflu = Bc2%ke1 - LAP; iflu = Bc2%ie1; jflu = Bc2%je1
      end select
      d1 = Bf%U(1,iflu,jflu,kflu)
      uu1 = Bf%U(2,iflu,jflu,kflu) / d1
      v1 = Bf%U(3,iflu,jflu,kflu) / d1
      w1 = Bf%U(4,iflu,jflu,kflu) / d1
      p1 = (Bf%U(5,iflu,jflu,kflu) - 0.5d0*d1*(uu1*uu1+v1*v1+w1*w1)) * (gamma-1.d0)
      T_f = gamma * Ma * Ma * p1 / max(d1, 1.d-20)
!     Fluid thermal conductivity (from laminar mu and Pr)
      k_f = Bf%mu(iflu,jflu,kflu) * gamma * Ma * Ma / (PrL * (gamma-1.d0))

!     Cell sizes at interface
      select case(face_s)
      case(1,4); dx_s = abs(Bs%x(is+1,js,ks) - Bs%x(is,js,ks))
                dx_f = abs(Bf%x(iflu+1,jflu,kflu) - Bf%x(iflu,jflu,kflu))
      case(2,5); dx_s = abs(Bs%y(is,js+1,ks) - Bs%y(is,js,ks))
                dx_f = abs(Bf%y(iflu,jflu+1,kflu) - Bf%y(iflu,jflu,kflu))
      case(3,6); dx_s = abs(Bs%z(is,js,ks+1) - Bs%z(is,js,ks))
                dx_f = abs(Bf%z(iflu,jflu,kflu+1) - Bf%z(iflu,jflu,kflu))
      end select

!     Interface temperature: heat flux balance
      T_i = (k_f*T_f/max(dx_f,1.d-20) + k_s*T_s/max(dx_s,1.d-20)) &
          / (k_f/max(dx_f,1.d-20) + k_s/max(dx_s,1.d-20))
     end subroutine compute_interface_T

!    Set solid ghost cell from interface temperature
!    face_s: 1=i-, 2=j-, 3=k-, 4=i+, 5=j+, 6=k+
!    ib,ie,jb,je,kb,ke: subface range in bc_msg format
     subroutine set_solid_ghost_from_interface_face(nMesh, mBlock_s, face_s, &
         ib,ie,jb,je,kb,ke, T_i, Bs)
      implicit none
      integer:: nMesh, mBlock_s, face_s
      integer:: ib,ie,jb,je,kb,ke
      real(PRE_EC):: T_i
      Type (Block_TYPE),pointer:: Bs
      integer:: n, i1,j1,k1, i2,j2,k2
      integer:: i,j,k

      select case(face_s)
      case(1)
        do n=1, LAP
          i1 = ib - n; i2 = ib + n - 1
          do k=kb,ke; do j=jb,je
            Bs%Ts(i1,j,k) = 2.d0*T_i - Bs%Ts(i2,j,k)
          enddo; enddo
        enddo
      case(4)
        do n=1, LAP
          i1 = ie + n - 1; i2 = ie - n
          do k=kb,ke; do j=jb,je
            Bs%Ts(i1,j,k) = 2.d0*T_i - Bs%Ts(i2,j,k)
          enddo; enddo
        enddo
      case(2)
        do n=1, LAP
          j1 = jb - n; j2 = jb + n - 1
          do k=kb,ke; do i=ib,ie
            Bs%Ts(i,j1,k) = 2.d0*T_i - Bs%Ts(i,j2,k)
          enddo; enddo
        enddo
      case(5)
        do n=1, LAP
          j1 = je + n - 1; j2 = je - n
          do k=kb,ke; do i=ib,ie
            Bs%Ts(i,j1,k) = 2.d0*T_i - Bs%Ts(i,j2,k)
          enddo; enddo
        enddo
      case(3)
        do n=1, LAP
          k1 = kb - n; k2 = kb + n - 1
          do j=jb,je; do i=ib,ie
            Bs%Ts(i,j,k1) = 2.d0*T_i - Bs%Ts(i,j,k2)
          enddo; enddo
        enddo
      case(6)
        do n=1, LAP
          k1 = ke + n - 1; k2 = ke - n
          do j=jb,je; do i=ib,ie
            Bs%Ts(i,j,k1) = 2.d0*T_i - Bs%Ts(i,j,k2)
          enddo; enddo
        enddo
      end select
     end subroutine set_solid_ghost_from_interface_face

!    Set fluid ghost cell as isothermal wall at T_i (face version)
!    face_s: 1=i-, 2=j-, 3=k-, 4=i+, 5=j+, 6=k+
!    ib,ie,jb,je,kb,ke: subface range in bc_msg format
     subroutine set_fluid_ghost_wall_face(nMesh, mBlock_f, &
         face_s, ib, ie, jb, je, kb, ke, T_i)
      implicit none
      integer:: nMesh, mBlock_f, face_s
      integer:: ib,ie,jb,je,kb,ke
      real(PRE_EC):: T_i
      Type (Block_TYPE),pointer:: Bf
      integer:: n, i1,j1,k1, i2,j2,k2
      integer:: i,j,k

      Bf => Mesh(nMesh)%Block(mBlock_f)

!     Use innermost ghost layer (n=1) for isothermal wall
      n = 1
      select case(face_s)
      case(1)
        i1 = ib + n - 1; i2 = ib + LAP + n - 1
        do k=kb,ke-1; do j=jb,je-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i2,j,k), Bf%U(:,i1,j,k), Ma, gamma, T_i, &
                                   Bf%mu(i2,j,k), Bf%dw(i2,j,k), Re)
        enddo; enddo
      case(4)
        i1 = ie - n + 1; i2 = ie - LAP - n + 1
        do k=kb,ke-1; do j=jb,je-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i2,j,k), Bf%U(:,i1,j,k), Ma, gamma, T_i, &
                                   Bf%mu(i2,j,k), Bf%dw(i2,j,k), Re)
        enddo; enddo
      case(2)
        j1 = jb + n - 1; j2 = jb + LAP + n - 1
        do k=kb,ke-1; do i=ib,ie-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i,j2,k), Bf%U(:,i,j1,k), Ma, gamma, T_i, &
                                   Bf%mu(i,j2,k), Bf%dw(i,j2,k), Re)
        enddo; enddo
      case(5)
        j1 = je - n + 1; j2 = je - LAP - n + 1
        do k=kb,ke-1; do i=ib,ie-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i,j2,k), Bf%U(:,i,j1,k), Ma, gamma, T_i, &
                                   Bf%mu(i,j2,k), Bf%dw(i,j2,k), Re)
        enddo; enddo
      case(3)
        k1 = kb + n - 1; k2 = kb + LAP + n - 1
        do j=jb,je-1; do i=ib,ie-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i,j,k2), Bf%U(:,i,j,k1), Ma, gamma, T_i, &
                                   Bf%mu(i,j,k2), Bf%dw(i,j,k2), Re)
        enddo; enddo
      case(6)
        k1 = ke - n + 1; k2 = ke - LAP - n + 1
        do j=jb,je-1; do i=ib,ie-1
          call wall_bound_with_Tw(NVAR1, Bf%U(:,i,j,k2), Bf%U(:,i,j,k1), Ma, gamma, T_i, &
                                   Bf%mu(i,j,k2), Bf%dw(i,j,k2), Re)
        enddo; enddo
      end select
     end subroutine set_fluid_ghost_wall_face

!   Incompressible (BLOCK_LOWSPEED) fluid - solid conjugate interface.
!   Bs = solid block, Bf = low-speed fluid block, Bc2 = interface entry as
!   seen from the solid block (Bc2%face = solid face, Bc2%face1 = fluid face).
!   For every face cell the interface temperature follows from the flux
!   balance  k_f*(T_f-T_i)/dx_f = k_s*(T_i-T_s)/dx_s, then ghost cells are
!   set on both sides:
!     solid ghost : Ts = 2*T_i - Ts(interior)
!     fluid ghost : no-slip + p zero-gradient + T = 2*T_i - T(interior)
!   Requires conformal interface grids (node indices 1:1), as the fluid block
!   machinery does.  One ghost layer is sufficient for the low-speed solver.
!   Incompressible (BLOCK_LOWSPEED) fluid - solid conjugate interface.
!   Bs = solid block, Bf = low-speed fluid block, Bc2 = interface entry as seen
!   from the solid block (Bc2%face = solid face, Bc2%face1 = fluid face).
!   Any of the 6x6 face combinations is supported; the tangential index
!   correspondence follows the connection descriptors L1,L2,L3.  Both sides
!   are in physical units (K, W/mK).
!     solid ghost : Ts = 2*T_i - Ts(interior)
!     fluid ghost : no-slip + p zero-gradient + T = 2*T_i - T(interior)
     subroutine couple_lowspeed_fluid_solid_face(nMesh, Bs, Bf, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: Bs, Bf
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1, nd_s, nd_f, d, m, td1, td2, c1, c2, o1, o2, n1, lo_s, lo_f
      integer:: kbar(3), kear(3), kbar1(3), kear1(3)
      integer:: mo(3), sg(3), Lv(3)
      integer:: cs(3), cf(3), gs(3), gf(3), fpt(3)
      integer:: fnode_s, fnode_f
      real(PRE_EC):: T_s, T_f, k_s, k_f, dx_s, dx_f, T_i

      face_s = Bc2%face
      face1  = Bc2%face1
      kbar  = (/Bc2%ib,  Bc2%jb,  Bc2%kb /)
      kear  = (/Bc2%ie,  Bc2%je,  Bc2%ke /)
      kbar1 = (/Bc2%ib1, Bc2%jb1, Bc2%kb1/)
      kear1 = (/Bc2%ie1, Bc2%je1, Bc2%ke1/)
      Lv    = (/Bc2%L1,  Bc2%L2,  Bc2%L3 /)
      k_s = Bs%solid_k
      k_f = LS_k
      if(k_s <= 0.d0 .or. k_f <= 0.d0) return

      nd_s = mod(face_s-1,3)+1
      do d=1,3
        mo(d) = abs(Lv(d)); sg(d) = sign(1,Lv(d))
      enddo
      nd_f = mo(nd_s)
      if(nd_f .lt. 1 .or. nd_f .gt. 3) then
        if(my_id == 0) print*, 'couple_lowspeed: bad connection table'
        return
      endif

      if(face_s <= 3) then
        cs(nd_s) = kbar(nd_s);   fnode_s = kbar(nd_s);  lo_s = 1
      else
        cs(nd_s) = kear(nd_s)-1; fnode_s = kear(nd_s);  lo_s = 0
      endif
      if(face1 <= 3) then
        cf(nd_f) = kbar1(nd_f);   fnode_f = kbar1(nd_f);  lo_f = 1
      else
        cf(nd_f) = kear1(nd_f)-1; fnode_f = kear1(nd_f);  lo_f = 0
      endif

      td1 = 0; td2 = 0
      do d=1,3
        if(d == nd_s) cycle
        if(td1 == 0) then; td1 = d; else; td2 = d; endif
      enddo

      do c2 = kbar(td2), kear(td2)-1
        cs(td2) = c2
        m = mo(td2); o2 = c2 - kbar(td2)
        if(sg(td2) > 0) then; cf(m) = kbar1(m) + o2
        else; cf(m) = (kear1(m)-1) - o2; endif
        do c1 = kbar(td1), kear(td1)-1
          cs(td1) = c1
          m = mo(td1); o1 = c1 - kbar(td1)
          if(sg(td1) > 0) then; cf(m) = kbar1(m) + o1
          else; cf(m) = (kear1(m)-1) - o1; endif

          T_s = Bs%Ts(cs(1),cs(2),cs(3))
          T_f = Bf%U(5,cf(1),cf(2),cf(3))
          fpt = cs; fpt(nd_s) = fnode_s
          dx_s = sqrt( (Bs%xc(cs(1),cs(2),cs(3))-Bs%x(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bs%yc(cs(1),cs(2),cs(3))-Bs%y(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bs%zc(cs(1),cs(2),cs(3))-Bs%z(fpt(1),fpt(2),fpt(3)))**2 )
          fpt = cf; fpt(nd_f) = fnode_f
          dx_f = sqrt( (Bf%xc(cf(1),cf(2),cf(3))-Bf%x(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bf%yc(cf(1),cf(2),cf(3))-Bf%y(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bf%zc(cf(1),cf(2),cf(3))-Bf%z(fpt(1),fpt(2),fpt(3)))**2 )
          dx_s = max(dx_s,1.d-20); dx_f = max(dx_f,1.d-20)

          T_i = (k_f*T_f/dx_f + k_s*T_s/dx_s) / (k_f/dx_f + k_s/dx_s)

          do n1 = 1, LAP
            gs = cs
            if(lo_s == 1) then; gs(nd_s) = cs(nd_s) - n1
            else; gs(nd_s) = cs(nd_s) + n1; endif
            Bs%Ts(gs(1),gs(2),gs(3)) = 2.d0*T_i - T_s
          enddo

          do n1 = 1, LAP
            gf = cf
            if(lo_f == 1) then; gf(nd_f) = cf(nd_f) - n1
            else; gf(nd_f) = cf(nd_f) + n1; endif
            Bf%U(1,gf(1),gf(2),gf(3)) = LS_rho
            Bf%U(2,gf(1),gf(2),gf(3)) = -Bf%U(2,cf(1),cf(2),cf(3))
            Bf%U(3,gf(1),gf(2),gf(3)) = -Bf%U(3,cf(1),cf(2),cf(3))
            Bf%U(4,gf(1),gf(2),gf(3)) = -Bf%U(4,cf(1),cf(2),cf(3))
            Bf%U(5,gf(1),gf(2),gf(3)) = 2.d0*T_i - T_f
            Bf%p(gf(1),gf(2),gf(3))   = Bf%p(cf(1),cf(2),cf(3))
          enddo
        enddo
      enddo
     end subroutine couple_lowspeed_fluid_solid_face

!   Compressible (BLOCK_FLUID) fluid - solid conjugate interface.
!   Bs = solid block, Bf = compressible fluid block, Bc2 = interface entry as
!   seen from the solid block.
!
!   Units: the compressible solver is non-dimensional (T* = T/T_inf with
!   T_inf the physical free-stream temperature; mu* = mu_dim/(rho_inf*U_inf*L)
!   so the freestream mu* = 1/Re).  The solid block is dimensional (K, W/mK,
!   m).  Conversion factors are derived from Ma/Re/gamma/T_inf/PrL with the
!   standard-air convention rho_inf = 1 kg/m^3, R = 287 J/(kg K):
!     U_inf = Ma*sqrt(gamma*R*T_inf), mu_ref = rho_inf*U_inf*Lscale/Re,
!     k_ref  = mu_ref*cp/PrL  (cp = gamma*R/(gamma-1))
!   so that the non-dimensional fluid conductivity equals mu* (constant Pr).
!
!   Any of the 6x6 face combinations is supported: the normal faces are read
!   from the connection entry and the tangential index correspondence is taken
!   from the connection descriptors L1,L2,L3 (see Convert_bc).  Conformal
!   interface grids with index-aligned (possibly reversed) directions required.
     subroutine couple_compressible_fluid_solid_face(nMesh, Bs, Bf, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: Bs, Bf
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1, nd_s, nd_f, d, m, td1, td2, c1, c2, o1, o2
      integer:: n1, lo_s, lo_f
      integer:: kbar(3), kear(3), kbar1(3), kear1(3)
      integer:: mo(3), sg(3), Lv(3)
      integer:: cs(3), cf(3), gs(3), gf(3), fpt(3)
      integer:: fnode_s, fnode_f
      real(PRE_EC):: T_s, T_f, k_s_star, k_f, dx_s, dx_f, T_i_star, T_i_K
      real(PRE_EC):: d1, uu1, v1, w1, p1
      real(PRE_EC),parameter:: R_AIR = 287.0d0, RHO_REF = 1.0d0
      real(PRE_EC):: a_ref, mu_ref, cp_ref, k_ref

!     General compressible (BLOCK_FLUID) - solid (BLOCK_SOLID) conjugate
!     interface.  Bs = solid block, Bf = compressible fluid block, Bc2 = the
!     interface entry as seen from the solid block.  Any of the 6x6 face
!     combinations is supported; the tangential index correspondence is taken
!     from the connection descriptors L1,L2,L3 (as produced by convert_inp).
!     Units: fluid non-dimensional (T*=T/T_inf, k*=mu*, lengths by Lscale);
!     solid dimensional (K, W/mK).  The interface balance is done in the
!     fluid non-dimensional system (common Lscale cancels in the weights).
      face_s = Bc2%face
      face1  = Bc2%face1
      kbar  = (/Bc2%ib,  Bc2%jb,  Bc2%kb /)
      kear  = (/Bc2%ie,  Bc2%je,  Bc2%ke /)
      kbar1 = (/Bc2%ib1, Bc2%jb1, Bc2%kb1/)
      kear1 = (/Bc2%ie1, Bc2%je1, Bc2%ke1/)
      Lv    = (/Bc2%L1,  Bc2%L2,  Bc2%L3 /)

      nd_s = mod(face_s-1,3)+1
      do d=1,3
        mo(d) = abs(Lv(d)); sg(d) = sign(1,Lv(d))
      enddo
      if(mo(nd_s) .lt. 1 .or. mo(nd_s) .gt. 3) then
        if(my_id == 0) print*, 'couple_compressible: bad connection table; face', face_s, face1
        return
      endif
      nd_f = mo(nd_s)

      k_s_star = Bs%solid_k
      if(k_s_star <= 0.d0) return
      a_ref  = sqrt(gamma*R_AIR*T_inf)
      mu_ref = RHO_REF*Ma*a_ref*max(Lscale,1.d-30)/Re
      cp_ref = gamma*R_AIR/(gamma-1.d0)
      k_ref  = mu_ref*cp_ref/max(PrL,1.d-30)
      k_s_star = k_s_star / k_ref

!     solid/fluid interior cell (in the normal direction) and face node
      if(face_s <= 3) then
        cs(nd_s) = kbar(nd_s);   fnode_s = kbar(nd_s);  lo_s = 1
      else
        cs(nd_s) = kear(nd_s)-1; fnode_s = kear(nd_s);  lo_s = 0
      endif
      if(face1 <= 3) then
        cf(nd_f) = kbar1(nd_f);   fnode_f = kbar1(nd_f);  lo_f = 1
      else
        cf(nd_f) = kear1(nd_f)-1; fnode_f = kear1(nd_f);  lo_f = 0
      endif

!     the two tangential directions of the solid interface
      td1 = 0; td2 = 0
      do d=1,3
        if(d == nd_s) cycle
        if(td1 == 0) then; td1 = d; else; td2 = d; endif
      enddo

      do c2 = kbar(td2), kear(td2)-1
        cs(td2) = c2
        m = mo(td2); o2 = c2 - kbar(td2)
        if(sg(td2) > 0) then; cf(m) = kbar1(m) + o2
        else; cf(m) = (kear1(m)-1) - o2; endif
        do c1 = kbar(td1), kear(td1)-1
          cs(td1) = c1
          m = mo(td1); o1 = c1 - kbar(td1)
          if(sg(td1) > 0) then; cf(m) = kbar1(m) + o1
          else; cf(m) = (kear1(m)-1) - o1; endif

          T_s = Bs%Ts(cs(1),cs(2),cs(3))
          d1  = Bf%U(1,cf(1),cf(2),cf(3))
          uu1 = Bf%U(2,cf(1),cf(2),cf(3))/d1
          v1  = Bf%U(3,cf(1),cf(2),cf(3))/d1
          w1  = Bf%U(4,cf(1),cf(2),cf(3))/d1
          p1  = (Bf%U(5,cf(1),cf(2),cf(3)) - 0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gamma-1.d0)
          T_f = gamma*Ma*Ma*p1/max(d1,1.d-20)          ! non-dimensional T*
          k_f = max(Bf%mu(cf(1),cf(2),cf(3)), 1.d-30)  ! non-dimensional k = mu*

!         normal distance cell-centre -> interface face node
          fpt = cs; fpt(nd_s) = fnode_s
          dx_s = sqrt( (Bs%xc(cs(1),cs(2),cs(3))-Bs%x(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bs%yc(cs(1),cs(2),cs(3))-Bs%y(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bs%zc(cs(1),cs(2),cs(3))-Bs%z(fpt(1),fpt(2),fpt(3)))**2 )
          fpt = cf; fpt(nd_f) = fnode_f
          dx_f = sqrt( (Bf%xc(cf(1),cf(2),cf(3))-Bf%x(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bf%yc(cf(1),cf(2),cf(3))-Bf%y(fpt(1),fpt(2),fpt(3)))**2 &
                     + (Bf%zc(cf(1),cf(2),cf(3))-Bf%z(fpt(1),fpt(2),fpt(3)))**2 )
          dx_s = max(dx_s, 1.d-20); dx_f = max(dx_f, 1.d-20)

!         interface temperature (non-dimensional heat flux balance)
          T_i_star = (k_f*T_f/dx_f + k_s_star*(T_s/T_inf)/dx_s) &
                   / (k_f/dx_f + k_s_star/dx_s)

!         solid ghost layers: isothermal at T_i (physical K)
          T_i_K = T_i_star*T_inf
          do n1 = 1, LAP
            gs = cs
            if(lo_s == 1) then; gs(nd_s) = cs(nd_s) - n1
            else; gs(nd_s) = cs(nd_s) + n1; endif
            Bs%Ts(gs(1),gs(2),gs(3)) = 2.d0*T_i_K - T_s
          enddo

!         fluid ghost layers: no-slip isothermal wall at T_i*
          do n1 = 1, LAP
            gf = cf
            if(lo_f == 1) then; gf(nd_f) = cf(nd_f) - n1
            else; gf(nd_f) = cf(nd_f) + n1; endif
            call wall_bound_with_Tw(NVAR1, Bf%U(:,cf(1),cf(2),cf(3)), &
                 Bf%U(:,gf(1),gf(2),gf(3)), Ma, gamma, T_i_star, &
                 Bf%mu(cf(1),cf(2),cf(3)), Bf%dw(cf(1),cf(2),cf(3)), Re)
          enddo
        enddo
      enddo
     end subroutine couple_compressible_fluid_solid_face

!===============================================================================
! Cross-process compressible-fluid <-> solid conjugate coupling.
! Called on BOTH processes that own the two blocks of an interface 11 pair when
! the neighbour is not local.  B is the local block (solid on one side, fluid on
! the other), Bc2 is B's own interface entry, nb the neighbour's global number.
! Two passes:
!   pass 1 (fluid -> solid): fluid sends (T_f*, k_f*, dx_f) per interface cell;
!   pass 2 (solid -> fluid): solid sends T_i* back; each side fills its own ghost.
! The cell correspondence is obtained from the connection descriptors L1,L2,L3
! (the tangential index of the solid region is mapped to the fluid region exactly
! as in Umessage_send_mpi, packing in the receiver's index space).
!===============================================================================
!===============================================================================
! Cross-process solid <-> fluid conjugate coupling.  The fluid block may be
! compressible (BLOCK_FLUID) or incompressible (BLOCK_LOWSPEED, AC/SIMPLE).
! Called on BOTH processes that own the interface pair when the neighbour is
! remote; B is the local block, Bc2 its own interface entry, nb the neighbour's
! global block number.  Two passes with the payload in PHYSICAL units
! (T [K], k [W/mK], dx):
!   pass 1 (fluid -> solid): fluid sends (T_f, k_f, dx_f) per interface cell;
!   pass 2 (solid -> fluid): solid computes T_i, sets its own ghost and sends
!                            T_i [K] back; the fluid fills its own ghost.
! The cell correspondence uses the connection descriptors L1,L2,L3 (packed in
! the receiver's index space, as in Umessage_send_mpi).
!===============================================================================
     subroutine cht_mpi_fluid_solid_face(nMesh, B, nb, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh, nb
      Type (Block_TYPE),pointer:: B
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: kbar(3),kear(3),kbar1(3),kear1(3),Lv(3),Pv(3),inv(3)
      integer:: sb_lo(3),sb_hi(3),fnb_lo(3),fnb_hi(3)
      integer:: fs, ff, nd_s, nd_f, td1, td2, d, q, n1, idx, ncell, nvar1
      integer:: c1, c2, ka, ksd, sl, lo_s, lo_f, fnode_s, fnode_f
      integer:: cs(3), cf(3), gs(3), gf(3), fpt(3)
      integer:: solid_blk, tag1, tag2, ierr
      integer:: status(MPI_status_size)
      real(PRE_EC):: T_s, T_f, kf, ks, dx_s, dx_f
      real(PRE_EC):: d1,uu1,v1,w1,p1,Ti
      real(PRE_EC),parameter:: R_AIR=287.0d0, RHO_REF=1.0d0
      real(PRE_EC):: a_ref,mu_ref,cp_ref,k_ref
      real(PRE_EC),allocatable:: sbuf(:), rbuf(:)

      nvar1 = Mesh(nMesh)%NVAR
      kbar  = (/Bc2%ib,  Bc2%jb,  Bc2%kb /)
      kear  = (/Bc2%ie,  Bc2%je,  Bc2%ke /)
      kbar1 = (/Bc2%ib1, Bc2%jb1, Bc2%kb1/)
      kear1 = (/Bc2%ie1, Bc2%je1, Bc2%ke1/)
      Lv    = (/Bc2%L1,  Bc2%L2,  Bc2%L3 /)
      do d=1,3; Pv(d)=sign(1,Lv(d)); enddo
      do q=1,3; inv(q)=0; enddo
      do d=1,3
        q=abs(Lv(d)); if(q>=1 .and. q<=3) inv(q)=d
      enddo

      if(B%Block_type == BLOCK_SOLID) then
        sb_lo=kbar;   sb_hi=kear;    fs=Bc2%face
        fnb_lo=kbar1; fnb_hi=kear1;  ff=Bc2%face1
        solid_blk = B%Block_no
      else
        sb_lo=kbar1;  sb_hi=kear1;   fs=Bc2%face1
        fnb_lo=kbar;  fnb_hi=kear;   ff=Bc2%face
        solid_blk = nb
      endif
      nd_s = mod(fs-1,3)+1
      nd_f = mod(ff-1,3)+1
      td1=0; td2=0
      do d=1,3
        if(d==nd_s) cycle
        if(td1==0) then; td1=d; else; td2=d; endif
      enddo
      ncell = (sb_hi(td2)-sb_lo(td2)) * (sb_hi(td1)-sb_lo(td1))
      if(ncell <= 0) return
      tag1 = 20000 + solid_blk
      tag2 = 22000 + solid_blk
      if(fs <= 3) then; fnode_s = sb_lo(nd_s); lo_s=1; else; fnode_s = sb_hi(nd_s); lo_s=0; endif
      if(ff <= 3) then; fnode_f = fnb_lo(nd_f); lo_f=1; else; fnode_f = fnb_hi(nd_f); lo_f=0; endif

      a_ref  = sqrt(gamma*R_AIR*T_inf)
      mu_ref = RHO_REF*Ma*a_ref*max(Lscale,1.d-30)/Re
      cp_ref = gamma*R_AIR/(gamma-1.d0)
      k_ref  = mu_ref*cp_ref/max(PrL,1.d-30)

!     ---------------------------------------------------------------- fluid
      if(B%Block_type /= BLOCK_SOLID) then
        allocate(sbuf(3*ncell))
        idx=0
        do c2=sb_lo(td2),sb_hi(td2)-1
        do c1=sb_lo(td1),sb_hi(td1)-1
          cs(td1)=c1; cs(td2)=c2; cs(nd_s)=sb_lo(nd_s)
          if(ff <= 3) then; cf(nd_f)=fnb_lo(nd_f); else; cf(nd_f)=fnb_hi(nd_f)-1; endif
          do q=1,3
            if(q==nd_s) cycle
            d=inv(q)
            if(Pv(d)>0) then; ksd=kbar(d); else; ksd=kear(d); endif
            ka = cs(q) - sb_lo(q)
            sl = ksd + ka*Pv(d)
            if(Pv(d)>0) then; cf(d)=sl; else; cf(d)=sl-1; endif
          enddo
          if(B%Block_type == BLOCK_LOWSPEED) then
            T_f = B%U(5,cf(1),cf(2),cf(3))          ! [K]
            kf  = LS_k                              ! [W/(m.K)]
          else
            d1=B%U(1,cf(1),cf(2),cf(3))
            uu1=B%U(2,cf(1),cf(2),cf(3))/d1
            v1=B%U(3,cf(1),cf(2),cf(3))/d1
            w1=B%U(4,cf(1),cf(2),cf(3))/d1
            p1=(B%U(5,cf(1),cf(2),cf(3))-0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gamma-1.d0)
            T_f = gamma*Ma*Ma*p1/max(d1,1.d-20)*T_inf                 ! [K]
            kf  = max(B%mu(cf(1),cf(2),cf(3)),1.d-30)*k_ref           ! [W/(m.K)]
          endif
          fpt=cf; fpt(nd_f)=fnode_f
          dx_f=sqrt( (B%xc(cf(1),cf(2),cf(3))-B%x(fpt(1),fpt(2),fpt(3)))**2 &
                   + (B%yc(cf(1),cf(2),cf(3))-B%y(fpt(1),fpt(2),fpt(3)))**2 &
                   + (B%zc(cf(1),cf(2),cf(3))-B%z(fpt(1),fpt(2),fpt(3)))**2 )
          dx_f=max(dx_f,1.d-20)
          sbuf(3*idx+1)=T_f; sbuf(3*idx+2)=kf; sbuf(3*idx+3)=dx_f
          idx=idx+1
        enddo
        enddo
        call MPI_Send(sbuf, 3*ncell, OCFD_DATA_TYPE, B_proc(nb), tag1, MPI_COMM_WORLD, ierr)
        allocate(rbuf(ncell))
        call MPI_Recv(rbuf, ncell, OCFD_DATA_TYPE, B_proc(nb), tag2, MPI_COMM_WORLD, status, ierr)
        idx=0
        do c2=sb_lo(td2),sb_hi(td2)-1
        do c1=sb_lo(td1),sb_hi(td1)-1
          cs(td1)=c1; cs(td2)=c2; cs(nd_s)=sb_lo(nd_s)
          if(ff <= 3) then; cf(nd_f)=fnb_lo(nd_f); else; cf(nd_f)=fnb_hi(nd_f)-1; endif
          do q=1,3
            if(q==nd_s) cycle
            d=inv(q)
            if(Pv(d)>0) then; ksd=kbar(d); else; ksd=kear(d); endif
            ka = cs(q) - sb_lo(q)
            sl = ksd + ka*Pv(d)
            if(Pv(d)>0) then; cf(d)=sl; else; cf(d)=sl-1; endif
          enddo
          T_f = rbuf(idx+1)
          do n1=1,LAP
            gf=cf
            if(lo_f==1) then; gf(nd_f)=cf(nd_f)-n1; else; gf(nd_f)=cf(nd_f)+n1; endif
            if(B%Block_type == BLOCK_LOWSPEED) then
              B%U(1,gf(1),gf(2),gf(3)) = LS_rho
              B%U(2,gf(1),gf(2),gf(3)) = -B%U(2,cf(1),cf(2),cf(3))
              B%U(3,gf(1),gf(2),gf(3)) = -B%U(3,cf(1),cf(2),cf(3))
              B%U(4,gf(1),gf(2),gf(3)) = -B%U(4,cf(1),cf(2),cf(3))
              B%U(5,gf(1),gf(2),gf(3)) = 2.d0*T_f - B%U(5,cf(1),cf(2),cf(3))
              B%p(gf(1),gf(2),gf(3))   = B%p(cf(1),cf(2),cf(3))
            else
              call wall_bound_with_Tw(nvar1, B%U(:,cf(1),cf(2),cf(3)), &
                   B%U(:,gf(1),gf(2),gf(3)), Ma, gamma, T_f/T_inf, &
                   B%mu(cf(1),cf(2),cf(3)), B%dw(cf(1),cf(2),cf(3)), Re)
            endif
          enddo
          idx=idx+1
        enddo
        enddo
        deallocate(sbuf,rbuf)
        return
      endif

!     ---------------------------------------------------------------- solid
      allocate(rbuf(3*ncell))
      call MPI_Recv(rbuf, 3*ncell, OCFD_DATA_TYPE, B_proc(nb), tag1, MPI_COMM_WORLD, status, ierr)
      allocate(sbuf(ncell))
      ks = B%solid_k
      idx=0
      do c2=sb_lo(td2),sb_hi(td2)-1
      do c1=sb_lo(td1),sb_hi(td1)-1
        cs(td1)=c1; cs(td2)=c2; cs(nd_s)=sb_lo(nd_s)
        T_s = B%Ts(cs(1),cs(2),cs(3))
        T_f = rbuf(3*idx+1); kf = rbuf(3*idx+2); dx_f = rbuf(3*idx+3)
        fpt=cs; fpt(nd_s)=fnode_s
        dx_s=sqrt( (B%xc(cs(1),cs(2),cs(3))-B%x(fpt(1),fpt(2),fpt(3)))**2 &
                 + (B%yc(cs(1),cs(2),cs(3))-B%y(fpt(1),fpt(2),fpt(3)))**2 &
                 + (B%zc(cs(1),cs(2),cs(3))-B%z(fpt(1),fpt(2),fpt(3)))**2 )
        dx_s=max(dx_s,1.d-20)
        Ti=(kf*T_f/dx_f + ks*T_s/dx_s)/(kf/dx_f + ks/dx_s)
        do n1=1,LAP
          gs=cs
          if(lo_s==1) then; gs(nd_s)=cs(nd_s)-n1; else; gs(nd_s)=cs(nd_s)+n1; endif
          B%Ts(gs(1),gs(2),gs(3))=2.d0*Ti - T_s
        enddo
        sbuf(idx+1)=Ti
        idx=idx+1
      enddo
      enddo
      call MPI_Send(sbuf, ncell, OCFD_DATA_TYPE, B_proc(nb), tag2, MPI_COMM_WORLD, ierr)
      deallocate(rbuf,sbuf)
     end subroutine cht_mpi_fluid_solid_face

!   High-speed (BLOCK_FLUID, compressible) <-> low-speed (BLOCK_LOWSPEED)
!   conjugate interface.  Each pair is processed ONCE from the compressible
!   block entry (B = compressible, Bn = low-speed).  The compressible
!   solution is non-dimensional (rho* = rho/rho_inf, u* = u/U_inf,
!   T* = T/T_inf, p* = p/(rho_inf*U_inf^2)); the low-speed solution is
!   physical SI (rho = LS_rho, u m/s, T K).  Units are unified with the
!   standard-air convention used elsewhere:
!     rho_inf = 1 kg/m^3, R = 287 J/(kgK), a_inf = sqrt(gamma*R*T_inf),
!     U_inf = Ma*a_inf,  p* = rho* T* / (gamma*Ma^2).
!   Ghost buffers are exchanged directly (ghost = converted value of the
!   neighbour's first interior cell), as requested.
!
!   Implemented for the current topology: compressible face_s = 2 (j-)
!   against low-speed face1 = 5 (j+).  Conformal aligned grids required.
     subroutine couple_highlow_fluid_face(nMesh, B, Bn, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: B, Bn
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1, i, k, jc_h, jc_l, jg_l, n1, n2, ic_h
      real(PRE_EC):: rho, uu, vv, TT, rho_s, u_s, v_s, T_s, p_nd, E_s
      real(PRE_EC),parameter:: R_AIR = 287.0d0
      real(PRE_EC),parameter:: MU_SI0 = 1.716d-5, T_SI0 = 273.15d0, S_SI = 110.4d0
      real(PRE_EC):: RHO_REF, mu_inf_p
      real(PRE_EC),parameter:: U_MAX_LS = 30.d0   ! m/s tangential cap (low-speed side)
      real(PRE_EC),parameter:: V_MAX_LS = 1.0d0   ! m/s normal cap (low-speed side)
      real(PRE_EC):: a_ref, U_ref
      real(PRE_EC):: U1,U2,U3,U5, r1, u1s, v1s, p1n, T1s
      real(PRE_EC):: p_anch
      logical:: has_pfix
      integer:: ksub2
      TYPE (BC_MSG_TYPE),pointer:: Bcl

      face_s = Bc2%face
      face1  = Bc2%face1
 !    Accepted topological pairs (both have the PHYSICAL interface normal along
 !    y, so the Cartesian component roles used below - U(2)=tangential x,
 !    U(3)=normal y - are identical and only the index bookkeeping differs):
 !      (2,5) gas j-  <-> low-speed j+     (legacy high_low_fluid grids)
 !      (4,2) gas i+  <-> low-speed j-     (cases/fluid_solid: the gas block2
 !                                          i=max face touches block1 j=1)
      if(.not. ((face_s == 2 .and. face1 == 5) .or. &
                (face_s == 4 .and. face1 == 2))) then
        print*, 'couple_highlow: supports comp j-/low j+ (2/5) or comp i+/low j-'// &
                ' (4/2); got', face_s, face1
        return
      endif
      a_ref = sqrt(gamma*R_AIR*T_inf)
      mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
      RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
      U_ref = Ma*a_ref

      if(face_s == 2) then
 !      ---- legacy: gas j- face against low-speed j+ face ------------------
 !      compressible block B (j- face): interior row jc_h = jb, ghost row jb-1
        jc_h = Bc2%jb
 !      low-speed block Bn (j+ face): interior cell row je1-1, ghosts je1..
        jc_l = Bc2%je1 - 1
        jg_l = Bc2%je1
        do k = Bc2%kb, Bc2%ke-1
        do i = Bc2%ib, Bc2%ie-1
 !        --- compressible ghost <- low-speed interior ---
          rho = Bn%U(1,i,jc_l,k)          ! physical density [kg/m3]
          uu  = Bn%U(2,i,jc_l,k)          ! velocity [m/s] (LS stores velocity in U(2..4))
          vv  = Bn%U(3,i,jc_l,k)
          TT  = Bn%U(5,i,jc_l,k)          ! temperature [K]
          rho_s = rho/RHO_REF
          u_s   = uu/U_ref
          v_s   = vv/U_ref
          T_s   = TT/T_inf
          p_nd  = rho_s*T_s/(gamma*Ma*Ma)
          E_s   = p_nd/(gamma-1.d0) + 0.5d0*rho_s*(u_s*u_s+v_s*v_s)
          B%U(1,i,jc_h-1,k) = rho_s
          B%U(2,i,jc_h-1,k) = rho_s*u_s
          B%U(3,i,jc_h-1,k) = rho_s*v_s
          B%U(4,i,jc_h-1,k) = 0.d0
          B%U(5,i,jc_h-1,k) = E_s
 !        --- low-speed ghost <- compressible interior ---
          U1 = B%U(1,i,jc_h,k)
          U2 = B%U(2,i,jc_h,k)
          U3 = B%U(3,i,jc_h,k)
          U5 = B%U(5,i,jc_h,k)
          r1  = max(U1,1.d-20)
          u1s = U2/r1
          v1s = U3/r1
          p1n = (U5 - 0.5d0*U1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          T1s = gamma*Ma*Ma*p1n/max(r1,1.d-20)
          rho = r1*RHO_REF
          uu  = u1s*U_ref
          vv  = v1s*U_ref
 !        The low-speed (incompressible) solver cannot sustain supersonic
 !        drag from the interface; clamp the exchanged velocity so the
 !        low-speed model stays in its valid (Ma << 1) regime.
          uu  = max(min(uu,  U_MAX_LS), -U_MAX_LS)
          vv  = max(min(vv,  V_MAX_LS), -V_MAX_LS)
          TT  = T1s*T_inf
          do n1 = jg_l, jg_l+LAP-1
            Bn%U(1,i,n1,k) = rho
 !          low-speed block stores VELOCITY in U(2..4) (U(1)=density); the old
 !          rho*uu form was only harmless while rho=1 kg/m^3
            Bn%U(2,i,n1,k) = uu
            Bn%U(3,i,n1,k) = vv
            Bn%U(4,i,n1,k) = 0.d0
            Bn%U(5,i,n1,k) = TT
            Bn%p(i,n1,k)   = Bn%p(i,jc_l,k)
          enddo
        enddo; enddo
      else
 !      ---- cases/fluid_solid: gas i+ face (i=ie) against low-speed j- face --
 !      compressible block B: interior column ic_h = ie-1, ghost columns
 !      ie, ie+1, ..., ie+LAP-1
        ic_h = Bc2%ie - 1
 !      low-speed block Bn (j- face): interior cell row jb1 (=1), ghost rows
 !      jb1-1, jb1-2, ... (the coolant enters through the far j=max face)
        jc_l = Bc2%jb1
        jg_l = Bc2%jb1 - 1
 !      tangential map: gas j cell <-> low-speed i cell (reversed; this mesh
 !      has the L=(2,-1,3) connection descriptor for the gas entry):
 !        ls_i = ie1-1 - (gas_j - jb)   ->   gas_j = jb + (ie1-1) - ls_i
        do k = Bc2%kb1, Bc2%ke1-1
        do i = Bc2%ib1, Bc2%ie1-1        ! low-speed tangential cell index
          i2 = Bc2%jb + (Bc2%ie1-1) - i  ! corresponding gas tangential cell
 !        --- compressible ghost (i=ie..) <- low-speed interior ---
          rho = Bn%U(1,i,jc_l,k)
          uu  = Bn%U(2,i,jc_l,k)
          vv  = Bn%U(3,i,jc_l,k)
          TT  = Bn%U(5,i,jc_l,k)
          rho_s = rho/RHO_REF
          u_s   = uu/U_ref
          v_s   = vv/U_ref
          T_s   = TT/T_inf
          p_nd  = rho_s*T_s/(gamma*Ma*Ma)
          E_s   = p_nd/(gamma-1.d0) + 0.5d0*rho_s*(u_s*u_s+v_s*v_s)
          do n1 = 0, LAP-1
            B%U(1,ic_h+1+n1,i2,k) = rho_s
            B%U(2,ic_h+1+n1,i2,k) = rho_s*u_s
            B%U(3,ic_h+1+n1,i2,k) = rho_s*v_s
            B%U(4,ic_h+1+n1,i2,k) = 0.d0
            B%U(5,ic_h+1+n1,i2,k) = E_s
          enddo
 !        --- low-speed ghost (j=jb1-1, jb1-2, ...) <- compressible interior ---
          U1 = B%U(1,ic_h,i2,k)
          U2 = B%U(2,ic_h,i2,k)
          U3 = B%U(3,ic_h,i2,k)
          U5 = B%U(5,ic_h,i2,k)
          r1  = max(U1,1.d-20)
          u1s = U2/r1
          v1s = U3/r1
          p1n = (U5 - 0.5d0*U1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          T1s = gamma*Ma*Ma*p1n/max(r1,1.d-20)
          rho = r1*RHO_REF
          uu  = u1s*U_ref
          vv  = v1s*U_ref
          uu  = max(min(uu,  U_MAX_LS), -U_MAX_LS)
          vv  = max(min(vv,  V_MAX_LS), -V_MAX_LS)
          TT  = T1s*T_inf
          do n1 = 0, LAP-1
            Bn%U(1,i,jg_l-n1,k) = rho
            Bn%U(2,i,jg_l-n1,k) = uu
            Bn%U(3,i,jg_l-n1,k) = vv
            Bn%U(4,i,jg_l-n1,k) = 0.d0
            Bn%U(5,i,jg_l-n1,k) = TT
            Bn%p(i,jg_l-n1,k)   = Bn%p(i,jc_l,k)
          enddo
        enddo; enddo
      endif

!     Anchor the low-speed pressure only when the block has NO pressure-
!     Dirichlet face (velocity inlet + walls + exchange face): with a
!     pressure inlet (LS_Inlet_Type=3) or pressure outlet the level is fixed
!     by that face, and the old whole-field shift (subtract p(1,1,1)) would
!     zero the interior pressure while the inlet ghost keeps p_face=LS_P_in,
!     producing a spurious O(LS_P_in) gradient at the inlet.
      has_pfix = .false.
      do ksub2 = 1, Bn%subface
        Bcl => Bn%bc_msg(ksub2)
        if(is_interface_bc(Bcl%bc)) cycle
        if(Bcl%bc == BC_Outflow .or. Bcl%bc == BC_LS_Outlet) then
          has_pfix = .true.; exit
        endif
        if((Bcl%bc == BC_Inflow .or. Bcl%bc == BC_LS_Inlet) .and. &
           LS_Inlet_Type == 3) then
          has_pfix = .true.; exit
        endif
      enddo
      if(.not. has_pfix) then
        p_anch = Bn%p(1,1,1)
        if(abs(p_anch) > 1.d-12) then
          do k = 1, Bn%nz-1
          do j = 1, Bn%ny-1
          do i = 1, Bn%nx-1
            Bn%p(i,j,k) = Bn%p(i,j,k) - p_anch
          enddo; enddo; enddo
        endif
      endif
     end subroutine couple_highlow_fluid_face
!----------------------------------------------------------------------
! Compressible (BLOCK_FLUID) <-> porous (BLOCK_POROUS) transpiration
! interface (code 19), blowing variant (Phase B).
!   - porous hot face (porous j+ face) is a coolant outlet: ghost layers get
!     zero-gradient state and the wall pressure from the compressible side
!     (SI), closing the porous SIMPLE pressure problem.
!   - compressible face (fluid j- face) is a wall with normal mass injection:
!     coolant at G = LS_rho*LS_U_in leaves the porous wall at temperature
!     T_w = Tf(porous hot cell) and enters the flow with v_w = G/rho_w,
!     rho_w = p_w/(R*T_w).  Tangential velocity mirrored (wall at rest).
!   Explicit staggered coupling per outer step (frame heat balance in phase C).
!   Conformal aligned grids: compressible j- (face_s=2) vs porous j+
!   (face1=5).  Indices/dimensions of ghost arrays as in couple_highlow.
     subroutine couple_compressible_porous_blowing_face(nMesh, B, Bn, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: B, Bn
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1, i, k, n1, jc_h, jc_l, jg_l, jg, ic_h, i2, ig
      real(PRE_EC):: G, T_w, v_w, rho_w, p_phys
      real(PRE_EC):: r1, u1s, v1s, p1n, v_w_s, v_int_s, vg_s, ug_s, T_s, rho_s, E_s
      real(PRE_EC):: RHO_REF, mu_inf_p, a_ref, U_ref
     real(PRE_EC),parameter:: R_AIR = 287.0d0
     logical:: pair42
      real(PRE_EC),parameter:: MU_SI0 = 1.716d-5, T_SI0 = 273.15d0, S_SI = 110.4d0

      face_s = Bc2%face
      face1  = Bc2%face1
!     Accepted topologies (the porous block is always the coolant side):
!       (2,5) gas j-  <-> porous j+   (porous_fluid_phaseB, high_low_fluid_800K)
!       (4,2) gas i+  <-> porous j-   (cases/fluid_solid: the gas block2 i=max
!                                      face touches block1 j=1; coolant enters
!                                      through the far j=max face)
      if(.not. ((face_s == 2 .and. face1 == 5) .or. &
                (face_s == 4 .and. face1 == 2))) then
        if(my_id == 0) print*, 'couple_compressible_porous_blowing: unsupported'// &
            ' face pair', face_s, face1, '(skipped)'
        return
      endif
      pair42 = (face_s == 4 .and. face1 == 2)
      a_ref = sqrt(gamma*R_AIR*T_inf)
      mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
      RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
      U_ref = Ma*a_ref
      G = LS_rho*max(LS_V_in,0.d0)          ! coolant mass flux kg/(m2 s), +y

      if(pair42) then
!      --- gas i+ face against porous j- face -------------------------------
!      gas: interior column ic_h = ie-1, ghost columns ig = ie..ie+LAP-1
!      (the gas normal direction is its i index, i.e. the physical y axis, so
!      U(3) carries the wall-normal velocity and U(2) the tangential one)
        ic_h = Bc2%ie - 1
!      porous: interior row jc_l = jb1 (=1), ghost rows jb1-1, jb1-2, ...
        jc_l = Bc2%jb1
        jg_l = Bc2%jb1 - 1
!      tangential map (L=(2,-1,3) for the gas entry): gas_j = jb + (ie1-1) - i
        do k = Bc2%kb1, Bc2%ke1-1
        do i = Bc2%ib1, Bc2%ie1-1
          i2 = Bc2%jb + (Bc2%ie1-1) - i
          r1  = max(B%U(1,ic_h,i2,k), 1.d-20)
          u1s = B%U(2,ic_h,i2,k)/r1
          v1s = B%U(3,ic_h,i2,k)/r1
          p1n = (B%U(5,ic_h,i2,k) - 0.5d0*r1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          p_phys = p1n*RHO_REF*U_ref*U_ref
          T_w = Bn%U(5,i,jc_l,k)            ! coolant exit temperature [K]
!         porous ghost layers: coolant outlet (zero-gradient) + wall pressure
          do n1 = 0, LAP-1
            jg = jg_l - n1
            Bn%U(1,i,jg,k)=Bn%U(1,i,jc_l,k); Bn%U(2,i,jg,k)=Bn%U(2,i,jc_l,k)
            Bn%U(3,i,jg,k)=Bn%U(3,i,jc_l,k); Bn%U(4,i,jg,k)=Bn%U(4,i,jc_l,k)
            Bn%U(5,i,jg,k)=Bn%U(5,i,jc_l,k); Bn%p(i,jg,k) = p_phys
          enddo
!         compressible ghost layers: blowing wall (temperature T_w)
          rho_w = p_phys/(R_AIR*max(T_w,1.d0))
          v_w   = G/max(rho_w,1.d-30)
          v_w_s = v_w/U_ref
          do n1 = 1, LAP
            ig = ic_h + n1
            vg_s = 2.d0*v_w_s - v1s
            ug_s = -u1s                     ! mirror tangential (wall at rest)
            T_s  = T_w/T_inf
            rho_s = p1n*gamma*Ma*Ma/max(T_s,1.d-30)
            E_s  = p1n/(gamma-1.d0) + 0.5d0*rho_s*(ug_s*ug_s+vg_s*vg_s)
            B%U(1,ig,i2,k)=rho_s; B%U(2,ig,i2,k)=rho_s*ug_s
            B%U(3,ig,i2,k)=rho_s*vg_s; B%U(4,ig,i2,k)=0.d0
            B%U(5,ig,i2,k)=E_s
          enddo
        enddo; enddo
      else
!      --- legacy: gas j- face against porous j+ face -----------------------
        jc_h = Bc2%jb          ! compressible first interior cell row (j- face)
        jc_l = Bc2%je1 - 1     ! porous first interior cell row (j+ face)
        jg_l = Bc2%je1         ! porous first ghost cell row (j+ face)

        do k = Bc2%kb, Bc2%ke-1
        do i = Bc2%ib, Bc2%ie-1
          r1  = max(B%U(1,i,jc_h,k), 1.d-20)
          u1s = B%U(2,i,jc_h,k)/r1
          v1s = B%U(3,i,jc_h,k)/r1
          p1n = (B%U(5,i,jc_h,k) - 0.5d0*B%U(1,i,jc_h,k)*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          p_phys = p1n*RHO_REF*U_ref*U_ref
          T_w = Bn%U(5,i,jc_l,k)              ! coolant exit temperature [K]
!         porous ghost layers: coolant outlet (zero-gradient) + wall pressure
          do n1 = 0, LAP-1
            jg = jg_l + n1
            Bn%U(1,i,jg,k)=Bn%U(1,i,jc_l,k); Bn%U(2,i,jg,k)=Bn%U(2,i,jc_l,k)
            Bn%U(3,i,jg,k)=Bn%U(3,i,jc_l,k); Bn%U(4,i,jg,k)=Bn%U(4,i,jc_l,k)
            Bn%U(5,i,jg,k)=Bn%U(5,i,jc_l,k); Bn%p(i,jg,k) = p_phys
          enddo
!         compressible ghost layers: blowing wall (temperature T_w)
          rho_w = p_phys/(R_AIR*max(T_w,1.d0))
          v_w   = G/max(rho_w,1.d-30)
          v_w_s = v_w/U_ref
          do n1 = 1, LAP
            jg = jc_h - n1
            v_int_s = v1s
            vg_s = 2.d0*v_w_s - v_int_s
            ug_s = -u1s                       ! mirror tangential (wall at rest)
            T_s  = T_w/T_inf
            rho_s = p1n*gamma*Ma*Ma/max(T_s,1.d-30)
            E_s  = p1n/(gamma-1.d0) + 0.5d0*rho_s*(ug_s*ug_s+vg_s*vg_s)
            B%U(1,i,jg,k)=rho_s; B%U(2,i,jg,k)=rho_s*ug_s
            B%U(3,i,jg,k)=rho_s*vg_s; B%U(4,i,jg,k)=0.d0
            B%U(5,i,jg,k)=E_s
          enddo
        enddo; enddo
      endif
     end subroutine couple_compressible_porous_blowing_face



  end subroutine couple_fluid_solid_interfaces

!----------------------------------------------------------------------
! Isothermal wall boundary condition with specified wall temperature
! (local version of wall_bound with Tw parameter instead of global Twall)
!----------------------------------------------------------------------
  subroutine wall_bound_with_Tw(NVAR, U1, Ug1, Ma, gamma, Tw, mu1, dw, Re)
   use precision_EC
   implicit none
   integer:: NVAR
   real(PRE_EC):: U1(NVAR), Ug1(NVAR)
   real(PRE_EC):: d1,uu1,v1,w1,T1,p1,d2,uu2,v2,w2,T2,p2,Ma,gamma,Tw,mu1,dw,wt,Re
   real(PRE_EC),parameter:: beta1_SST=0.075d0

   if(Tw .gt. 0.d0) then
     d1=U1(1)
     uu1=U1(2)/d1
     v1=U1(3)/d1
     w1=U1(4)/d1
     p1=(U1(5)-0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gamma-1.d0)
     T1=gamma*Ma*Ma*p1/d1
     p2=p1
     T2=2.d0*Tw-T1
     if(T2 .lt. 0.5d0*T1) T2=0.5d0*T1
     uu2=-uu1
     v2=-v1
     w2=-w1
     d2=gamma*Ma*Ma*p2/T2
     Ug1(1)=d2
     Ug1(2)=d2*uu2
     Ug1(3)=d2*v2
     Ug1(4)=d2*w2
     Ug1(5)=p2/(gamma-1.d0)+0.5d0*d2*(uu2*uu2+v2*v2+w2*w2)
   else
!    Adiabatic wall
     Ug1(1)=U1(1)
     Ug1(2:4)=-U1(2:4)
     Ug1(5)=U1(5)
   endif

!  Turbulence model ghost cell values
   if(NVAR .ge. 6) then
     Ug1(6)=U1(6)
   endif
   if(NVAR .ge. 7) then
     Ug1(7)=U1(7)
   endif
  end subroutine wall_bound_with_Tw

!----------------------------------------------------------------------
! Solid-solid interface coupling
! For each bc=12 interface, directly copy Ts from the source block's
! interior to the destination block's ghost cells.
!----------------------------------------------------------------------
  subroutine couple_solid_solid_interfaces(nMesh)
   use Global_Var
   implicit none
   integer:: nMesh
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: mBlock, ksub, nb, mb, face_s, n
   integer:: ib,ie,jb,je,kb,ke, i1,j1,k1, i2,j2,k2
   integer:: i,j,k

   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%bc_msg2)) cycle

     do ksub=1, B%subface
       Bc2 => B%bc_msg2(ksub)
!      Skip physical boundaries (bc >= 0); only process interface connections (bc < 0)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1             ! Neighbor block number (global)
       if(nb <= 0) cycle
       if(B_proc(nb) .ne. my_id) cycle   ! neighbour on another process
       mb = B_n(nb)
       Bn => Mesh(nMesh)%Block(mb)

!      Determine if this is a solid-solid interface
       if(.not. (B%Block_type == BLOCK_SOLID .and. Bn%Block_type == BLOCK_SOLID)) cycle

!      Copy Ts from neighbor's interior to this block's ghost cells
       face_s = Bc2%face
       ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke

       select case(face_s)
       case(1)
         do n=1, LAP
           i1 = ib + n - 1; i2 = ib + LAP + n - 1
           do k=kb,ke; do j=jb,je
             B%Ts(i1,j,k) = Bn%Ts(i2,j,k)
           enddo; enddo
         enddo
       case(4)
         do n=1, LAP
           i1 = ie - n + 1; i2 = ie - LAP - n + 1
           do k=kb,ke; do j=jb,je
             B%Ts(i1,j,k) = Bn%Ts(i2,j,k)
           enddo; enddo
         enddo
       case(2)
         do n=1, LAP
           j1 = jb + n - 1; j2 = jb + LAP + n - 1
           do k=kb,ke; do i=ib,ie
             B%Ts(i,j1,k) = Bn%Ts(i,j2,k)
           enddo; enddo
         enddo
       case(5)
         do n=1, LAP
           j1 = je - n + 1; j2 = je - LAP - n + 1
           do k=kb,ke; do i=ib,ie
             B%Ts(i,j1,k) = Bn%Ts(i,j2,k)
           enddo; enddo
         enddo
       case(3)
         do n=1, LAP
           k1 = kb + n - 1; k2 = kb + LAP + n - 1
           do j=jb,je; do i=ib,ie
             B%Ts(i,j,k1) = Bn%Ts(i,j,k2)
           enddo; enddo
         enddo
       case(6)
         do n=1, LAP
           k1 = ke - n + 1; k2 = ke - LAP - n + 1
           do j=jb,je; do i=ib,ie
             B%Ts(i,j,k1) = Bn%Ts(i,j,k2)
           enddo; enddo
         enddo
       end select
     enddo
   enddo
  end subroutine couple_solid_solid_interfaces

!----------------------------------------------------------------------
! Low-speed (BLOCK_LOWSPEED) <-> porous (BLOCK_POROUS) interface coupling
! (explicit interface code BC_Interface_LowPorous = 16).
!
! Both solvers store the SAME incompressible primitive layout on the cells:
!   U(1)   = rho  (density, constant = LS_rho)
!   U(2:4) = rho*(u,v,w)  (momentum, SI units)
!   U(5)   = Tf           (fluid temperature, K)
!   B%p    = pressure     (Pa, cell-centred)
! The porous block additionally carries the solid-frame temperature B%Ts,
! which is NOT exchanged here (isothermal low-porous test cases keep the
! inter-phase heat exchange off, hv=0; the porous ghost Ts keeps its
! zero-gradient initial value on the interface face).
!
! For every interface face the LAP ghost cells of block B are filled with
! the neighbour block's first interior cell values so that the tangential
! velocity, pressure and temperature are continuous across the face.
! The normal velocity on a conformal flat interface is ~0 in the fully
! developed Beavers-Joseph test, so the interface face mass flux stays 0 on
! both sides and only the tangential momentum (viscous shear) is coupled
! through the ghosts.
!
! Supported: conformal, axis-aligned, index-aligned block faces (1:1
! tangential node ranges, opposite low/high normal orientation).  A pair is
! processed once from each side (like couple_solid_solid_interfaces).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Low-speed (BLOCK_LOWSPEED) <-> porous (BLOCK_POROUS) interface coupling
! (explicit interface code BC_Interface_LowPorous = 16).
! Both solvers store the SAME cell layout (U(1)=rho=LS_rho, U(2:4)=velocity,
! U(5)=T, p=B%p).  For every interface face the LAP ghost layers of block B
! are filled from the neighbour's first interior cell so that tangential
! velocity, pressure and temperature are continuous.  Arbitrary index-aligned
! (possibly rotated / reversed) face pairs are handled through the connection
! descriptors L1,L2,L3; remote neighbours are exchanged with MPI.
!----------------------------------------------------------------------
  subroutine couple_lowspeed_porous_interfaces(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: mBlock, ksub, nb, mb, nbt

   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%bc_msg2)) cycle
     do ksub=1, B%subface
       Bc2 => B%bc_msg2(ksub)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1
       if(nb <= 0) cycle
       nbt = Block_Type_List(nb)
       if(.not. ((B%Block_type == BLOCK_LOWSPEED .and. nbt == BLOCK_POROUS) .or. &
                 (B%Block_type == BLOCK_POROUS  .and. nbt == BLOCK_LOWSPEED))) cycle
       if(B_proc(nb) .ne. my_id) then
         call ls_por_iface_mpi(nMesh, B, nb, Bc2)
       else
         mb = B_n(nb)
         if(mb <= 0) cycle
         Bn => Mesh(nMesh)%Block(mb)
         call ls_por_iface_local(nMesh, B, Bn, Bc2)
       endif
     enddo
   enddo

   contains

!  fill block B's interface ghost layers from the (same-process) neighbour Bn
   subroutine ls_por_iface_local(nMesh, B, Bn, Bc2)
     use Global_Var
     implicit none
     integer:: nMesh
     Type (Block_TYPE),pointer:: B, Bn
     TYPE (BC_MSG_TYPE),pointer:: Bc2
     integer:: face_s, face1, nd_n, d, m, td1, td2, c1, c2, o1, o2, n1, lo_s
     integer:: kbar(3),kear(3),kbar1(3),kear1(3),mo(3),sg(3),Lv(3)
     integer:: cs(3), cf(3), gs(3)
     face_s = Bc2%face
     face1  = Bc2%face1
     kbar  = (/Bc2%ib,  Bc2%jb,  Bc2%kb /)
     kear  = (/Bc2%ie,  Bc2%je,  Bc2%ke /)
     kbar1 = (/Bc2%ib1, Bc2%jb1, Bc2%kb1/)
     kear1 = (/Bc2%ie1, Bc2%je1, Bc2%ke1/)
     Lv    = (/Bc2%L1,  Bc2%L2,  Bc2%L3 /)
     do d=1,3; mo(d)=abs(Lv(d)); sg(d)=sign(1,Lv(d)); enddo
     nd_n = mo(mod(face_s-1,3)+1)
     if(nd_n .lt. 1 .or. nd_n .gt. 3) return
     if(face_s <= 3) then; cs(mod(face_s-1,3)+1) = kbar(mod(face_s-1,3)+1); lo_s=1
     else; cs(mod(face_s-1,3)+1) = kear(mod(face_s-1,3)+1)-1; lo_s=0; endif
     if(face1 <= 3) then; cf(nd_n) = kbar1(nd_n); else; cf(nd_n) = kear1(nd_n)-1; endif
     td1=0; td2=0
     do d=1,3
       if(d == mod(face_s-1,3)+1) cycle
       if(td1 == 0) then; td1=d; else; td2=d; endif
     enddo
     do c2=kbar(td2),kear(td2)-1
       cs(td2)=c2; m=mo(td2); o2=c2-kbar(td2)
       if(sg(td2) > 0) then; cf(m)=kbar1(m)+o2; else; cf(m)=(kear1(m)-1)-o2; endif
       do c1=kbar(td1),kear(td1)-1
         cs(td1)=c1; m=mo(td1); o1=c1-kbar(td1)
         if(sg(td1) > 0) then; cf(m)=kbar1(m)+o1; else; cf(m)=(kear1(m)-1)-o1; endif
         do n1=1,LAP
           gs=cs
           if(lo_s == 1) then; gs(mod(face_s-1,3)+1)=cs(mod(face_s-1,3)+1)-n1
           else; gs(mod(face_s-1,3)+1)=cs(mod(face_s-1,3)+1)+n1; endif
           B%U(1,gs(1),gs(2),gs(3)) = Bn%U(1,cf(1),cf(2),cf(3))
           B%U(2,gs(1),gs(2),gs(3)) = Bn%U(2,cf(1),cf(2),cf(3))
           B%U(3,gs(1),gs(2),gs(3)) = Bn%U(3,cf(1),cf(2),cf(3))
           B%U(4,gs(1),gs(2),gs(3)) = Bn%U(4,cf(1),cf(2),cf(3))
           B%U(5,gs(1),gs(2),gs(3)) = Bn%U(5,cf(1),cf(2),cf(3))
           B%p(gs(1),gs(2),gs(3))   = Bn%p(cf(1),cf(2),cf(3))
         enddo
       enddo
     enddo
   end subroutine ls_por_iface_local

!  cross-process exchange of the interface cell values (U(1:5), p)
   subroutine ls_por_iface_mpi(nMesh, B, nb, Bc2)
     use Global_Var
     implicit none
     integer:: nMesh, nb
     Type (Block_TYPE),pointer:: B
     TYPE (BC_MSG_TYPE),pointer:: Bc2
     integer:: kbar(3),kear(3),kbar1(3),kear1(3),Lv(3),Pv(3),inv(3)
     integer:: own_nd, nbr_nd, td1,td2, nt1,nt2, d, q, n1, idx, ncell
     integer:: c1,c2, ka, ksd, sl, cs(3), nv(3), gs(3)
     integer:: tag, ierr, status(MPI_status_size)
     real(PRE_EC),allocatable:: sbuf(:), rbuf(:)
     kbar  = (/Bc2%ib,  Bc2%jb,  Bc2%kb /)
     kear  = (/Bc2%ie,  Bc2%je,  Bc2%ke /)
     kbar1 = (/Bc2%ib1, Bc2%jb1, Bc2%kb1/)
     kear1 = (/Bc2%ie1, Bc2%je1, Bc2%ke1/)
     Lv    = (/Bc2%L1,  Bc2%L2,  Bc2%L3 /)
     do d=1,3; Pv(d)=sign(1,Lv(d)); enddo
     do q=1,3; inv(q)=0; enddo
     do d=1,3; q=abs(Lv(d)); if(q>=1 .and. q<=3) inv(q)=d; enddo
     own_nd = mod(Bc2%face-1,3)+1
     nbr_nd = mod(Bc2%face1-1,3)+1
     td1=0; td2=0
     do d=1,3
       if(d == own_nd) cycle
       if(td1 == 0) then; td1=d; else; td2=d; endif
     enddo
     nt1=0; nt2=0
     do d=1,3
       if(d == nbr_nd) cycle
       if(nt1 == 0) then; nt1=d; else; nt2=d; endif
     enddo
     ncell = (kear(td2)-kbar(td2))*(kear(td1)-kbar(td1))
     if(ncell <= 0) return
     tag = 30000 + min(B%Block_no, nb)
     allocate(sbuf(6*ncell)); allocate(rbuf(6*ncell))
!    pack our interface cell values in the receiver's index space
     idx=0
     do c2=kbar1(nt2),kear1(nt2)-1
       nv(nt2)=c2
       do c1=kbar1(nt1),kear1(nt1)-1
         nv(nt1)=c1; nv(nbr_nd)=kbar1(nbr_nd)
         do q=1,3
           if(q == nbr_nd) cycle
           d=inv(q)
           if(Pv(d) > 0) then; ksd=kbar(d); else; ksd=kear(d); endif
           ka = nv(q) - kbar1(q)
           sl = ksd + ka*Pv(d)
           if(Pv(d) > 0) then; cs(d)=sl; else; cs(d)=sl-1; endif
         enddo
         if(Bc2%face <= 3) then; cs(own_nd)=kbar(own_nd); else; cs(own_nd)=kear(own_nd)-1; endif
         sbuf(6*idx+1)=B%U(1,cs(1),cs(2),cs(3))
         sbuf(6*idx+2)=B%U(2,cs(1),cs(2),cs(3))
         sbuf(6*idx+3)=B%U(3,cs(1),cs(2),cs(3))
         sbuf(6*idx+4)=B%U(4,cs(1),cs(2),cs(3))
         sbuf(6*idx+5)=B%U(5,cs(1),cs(2),cs(3))
         sbuf(6*idx+6)=B%p(cs(1),cs(2),cs(3))
         idx=idx+1
       enddo
     enddo
     call MPI_Sendrecv(sbuf, 6*ncell, OCFD_DATA_TYPE, B_proc(nb), tag, &
                       rbuf, 6*ncell, OCFD_DATA_TYPE, B_proc(nb), tag, &
                       MPI_COMM_WORLD, status, ierr)
!    fill our ghost layers (received array indexed in our own region)
     idx=0
     do c2=kbar(td2),kear(td2)-1
       do c1=kbar(td1),kear(td1)-1
         do n1=1,LAP
           gs=0; gs(td1)=c1; gs(td2)=c2
           if(Bc2%face <= 3) then; gs(own_nd)=kbar(own_nd)-n1; else; gs(own_nd)=kear(own_nd)-1+n1; endif
           B%U(1,gs(1),gs(2),gs(3))=rbuf(6*idx+1)
           B%U(2,gs(1),gs(2),gs(3))=rbuf(6*idx+2)
           B%U(3,gs(1),gs(2),gs(3))=rbuf(6*idx+3)
           B%U(4,gs(1),gs(2),gs(3))=rbuf(6*idx+4)
           B%U(5,gs(1),gs(2),gs(3))=rbuf(6*idx+5)
           B%p(gs(1),gs(2),gs(3))  =rbuf(6*idx+6)
         enddo
         idx=idx+1
       enddo
     enddo
     deallocate(sbuf,rbuf)
   end subroutine ls_por_iface_mpi

  end subroutine couple_lowspeed_porous_interfaces


!----------------------------------------------------------------------
! (porous_solver_one_block is implemented in sub_porous.f90)
!----------------------------------------------------------------------
! Output solid-frame temperature field (Ts) for solid and porous blocks
! Writes one file per block: Ts_block_NNN.dat (formatted, PLOT3D-like)
!----------------------------------------------------------------------
  subroutine output_Ts
   use Global_Var
   implicit none
   integer:: mBlock, nMesh, i, j, k, nx, ny, nz
   Type (Block_TYPE),pointer:: B
   character(len=32):: fname

   nMesh = 1
   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type /= BLOCK_SOLID .and. B%Block_type /= BLOCK_POROUS) cycle
     nx = B%nx; ny = B%ny; nz = B%nz
     write(fname, '("Ts_block_",I0,".dat")') mBlock
     open(201, file=fname, status='replace')
     write(201,*) 'TITLE = "Solid temperature field"'
     write(201,*) 'VARIABLES = "x" "y" "z" "Ts"'
     write(201,*) 'ZONE I=', nx-1, ' J=', ny-1, ' K=', nz-1, ' DATAPACKING=POINT'
     do k=1, nz-1
       do j=1, ny-1
         do i=1, nx-1
           write(201, '(4ES20.10)') B%xc(i,j,k), B%yc(i,j,k), B%zc(i,j,k), B%Ts(i,j,k)
         enddo
       enddo
     enddo
     close(201)
     if(my_id == 0) print*, 'Output Ts to ', trim(fname)
   enddo
  end subroutine output_Ts

!==============================================================================
! Staggered segmented coupling for FLUID<->POROUS interface code 19
! (Iflag_Couple_Scheme=1, thermal wall / transpiration-cooled porous wall).
!
! Outer loop (steady, segregated / staggered):
!   for it = 1 .. Niter_Couple_Outer
!     (1) GAS CHUNK: advance only the compressible (BLOCK_FLUID) blocks for
!         Kstep_Couple_Comp time steps.  The code-19 subfaces of the gas block
!         (the "coupling" region, in this topology the whole j- wall) are
!         imposed as a no-slip isothermal wall at the per-face temperature
!         fp_Tw(i,k) (first chunk: Twall_Couple_Init, e.g. 300 K).  Ordinary
!         code-2 wall segments stay adiabatic (standard Twall<0 handling).
!         The porous block is NOT advanced in this phase.
!     (2) EXTRACT: evaluate the per-face wall heat flux fp_qw(i,k) [W/m2,
!         >0 into the wall] and wall static pressure fp_pw(i,k) [Pa].
!     (3) POROUS CHUNK: advance only the porous block to convergence.  At its
!         hot face (same code-19 face, porous j+): fluid outlet ghost pressure
!         = fp_pw, velocity/T zero-gradient; solid-frame Ts ghost enforces the
!         incoming heat flux fp_qw (Ts_g = Ts_i + qw*dx/ks_eff).  Coolant
!         supply at the porous underside keeps physical BCs from control.ec.
!     (4) RETURN: new hot-face temperature fp_Tw = (Ts_i + Ts_g)/2 -> gas
!         chunk; repeat until max|dT_w| < Tol_Couple_Tw.
!
! Supports one conformal pair: compressible face_s=2 (j-) vs porous
! face1=5 (j+), aligned indices (porous_fluid_phaseB topology).
!==============================================================================
  subroutine run_staggered_fluid_porous(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   integer:: mf, mp, mBlock, mBlock2, ksub, nb, mb, it, step, pc, NVAR1
   integer:: nstep, nhalve
   integer:: face_s, face1, ib,ie,jb,je,kb,ke, ib1,ie1,jb1,je1,kb1,ke1
   Type (Block_TYPE),pointer:: Bf, Bp, B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: i, k, jc_h, jc_l, jg, n1, ic_h, jg_l, gi
  real(PRE_EC):: Sfac, Sfac1, twmax, qwmax, xw, vcmax
   real(PRE_EC):: a_ref, U_ref, RHO_REF, mu_inf_p, mu_ref, cp_ref
   real(PRE_EC),parameter:: R_AIR=287.d0
   real(PRE_EC),parameter:: MU_SI0=1.716d-5, T_SI0=273.15d0, S_SI=110.4d0
   real(PRE_EC):: r1, u1s, v1s, p1n, T1_nd, T1_K, mu_SI, k_gas, dxp
   real(PRE_EC):: mu1c, dwc
   logical:: found, converged, pair42
   real(PRE_EC):: Utmp(7)

!  ---- locate compressible and porous blocks --------------------------------
   mf=0; mp=0
   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type == BLOCK_FLUID .and. mf == 0) mf = mBlock
     if(B%Block_type == BLOCK_POROUS .and. mp == 0) mp = mBlock
   enddo
   if(mf == 0 .or. mp == 0) then
     print*, 'run_staggered_fluid_porous: need one BLOCK_FLUID and one BLOCK_POROUS block, got', mf, mp
     return
   endif
   Bp => Mesh(nMesh)%Block(mp)
   NVAR1 = Mesh(nMesh)%NVAR

!  ---- find the FLUID code-19 subface (coupling region).  In multi-block gas
!       topologies (e.g. high_low_fluid_800K, 4 blocks: one porous + three
!       compressible) the code-19 subface may live on any one of the
!       BLOCK_FLUID blocks, not necessarily the first one -- scan them all.
   found = .false.
   do mb=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mb)
     if(B%Block_type /= BLOCK_FLUID) cycle
     do ksub=1, B%subface
       Bc2 => B%bc_msg2(ksub)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1
       if(nb <= 0) cycle
       mBlock2 = 0
       do mBlock=1, Mesh(nMesh)%Num_Block
         if(Mesh(nMesh)%Block(mBlock)%Block_no == nb) then
           Bn => Mesh(nMesh)%Block(mBlock); mBlock2 = mBlock; exit
         endif
       enddo
       if(mBlock2 == 0) cycle
       if(Bn%Block_type /= BLOCK_POROUS) cycle
      face_s = Bc2%face; face1 = Bc2%face1
!     Accepted topologies:
!       (2,5) gas j-  <-> porous j+   (porous_fluid_phaseB / high_low_fluid_800K)
!       (4,2) gas i+  <-> porous j-   (cases/fluid_solid: gas block2 i=max
!                                      touches block1 j=1)
      if(.not. ((face_s == 2 .and. face1 == 5) .or. &
                (face_s == 4 .and. face1 == 2))) then
        if(my_id == 0) print*, 'run_staggered: supports gas j-(2)/porous j+(5)'// &
            ' or gas i+(4)/porous j-(2); got', face_s, face1
        return
      endif
      pair42 = (face_s == 4 .and. face1 == 2)
      ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
      ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1
      Bf => Mesh(nMesh)%Block(mb); mf = mb
      found = .true.
      exit
    enddo
    if(found) exit
  enddo
  if(.not. found) then
    print*, 'run_staggered_fluid_porous: no FLUID-POROUS (code 19) interface found'
    return
  endif

!  ---- interface index descriptors -------------------------------------------
!   pair42 = .false. : gas j- face (interior row jc_h=jb, ghost rows jb-1..)
!                      against porous j+ face (interior row jc_l=je1-1,
!                      ghost rows jg_l.. = je1..)
!   pair42 = .true.  : gas i+ face (interior column ic_h=ie-1, ghost columns
!                      ie..) against porous j- face (interior row jc_l=jb1=1,
!                      ghost rows jg_l.. = jb1-1, jb1-2, ...)
!   The fp_* work arrays are always indexed by the POROUS interface cell
!   (i = ib1..ie1-1, k = kb1..ke1-1); the corresponding gas tangential cell is
!     gas_j = jb + (ie1-1) - ls_i        (reversed L=(2,-1,3) descriptor)
  if(pair42) then
    ic_h = ie - 1
    jc_l = jb1
    jg_l = jb1 - 1
  else
    ic_h = 0
    jc_h = jb
    jc_l = je1 - 1
    jg_l = je1
  endif

!  ---- per-face work arrays (porous interface cells) -------------------------
  if(.not. allocated(fp_Tw)) then
    allocate(fp_Tw(1:max(Bf%nx,Bp%nx)+1, 1:max(Bf%nz,Bp%nz)+1))
    allocate(fp_Tw_old(1:max(Bf%nx,Bp%nx)+1, 1:max(Bf%nz,Bp%nz)+1))
    allocate(fp_qw(1:max(Bf%nx,Bp%nx)+1, 1:max(Bf%nz,Bp%nz)+1))
    allocate(fp_pw(1:max(Bf%nx,Bp%nx)+1, 1:max(Bf%nz,Bp%nz)+1))
  endif
  fp_Tw = Twall_Couple_Init
  fp_qw = 0.d0
  fp_pw = 0.d0
!  重启：trailer 里若带着本界面（pair 19）的量，直接恢复（零跳变）；
!  Iflag_Couple_Restart=-1 时忽略（=旧行为）。
  if(Iflag_Couple_Restart >= 0 .and. Couple_State_Found == 1 .and. &
     Couple_State_Pair == 19 .and. Couple_State_nf == size(fp_Tw,1) .and. &
     Couple_State_nk == size(fp_Tw,2) .and. allocated(Couple_Tw_save)) then
    fp_Tw = Couple_Tw_save; fp_qw = Couple_qw_save; fp_pw = Couple_pw_save
    fp_Tw_old = fp_Tw
    if(my_id == 0) print*, ' restart: porous-interface T_w/q_w/p_w restored', &
           ' from the restart file (no re-seeding; last metric=', Couple_State_Twmax, ')'
  endif

!  ---- physical references (same convention as couple_highlow) --------------
  a_ref   = sqrt(gamma*R_AIR*T_inf)
  mu_inf_p= MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
  RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
  U_ref   = Ma*a_ref
  mu_ref  = RHO_REF*U_ref*max(Lscale,1.d-30)/Re
  cp_ref  = gamma*R_AIR/(gamma-1.d0)

  !  restart file: keep the saved counters instead of resetting
  if(restart_found == 1) then
    Mesh(nMesh)%tt    = restart_tt_saved
    Mesh(nMesh)%Kstep = restart_Kstep_saved
    if(my_id == 0) print*, ' restart: continue staggered run at Kstep=', &
        restart_Kstep_saved, ' tt=', restart_tt_saved
  !    NOTE: the per-face interface quantities (fp_Tw/fp_u) are re-seeded
  !    from the control.ec initial guess; the outer iterations re-converge
  !    them within a few rounds (state itself is fully restored).
  else
  Mesh(nMesh)%tt = 0.d0
  Mesh(nMesh)%Kstep = 0
  endif
  if(my_id == 0) then
    print*, ' run_staggered_fluid_porous: Kstep_Couple_Comp=', Kstep_Couple_Comp, &
            ' Niter_Couple_Outer=', Niter_Couple_Outer, &
            ' Twall_Couple_Init=', Twall_Couple_Init, ' K'
    print*, '   gas block', mf, ' porous block', mp, ' pair42(gas i+/porous j-)=', pair42
    print*, '   interface cells (porous index) i=', ib1, '..', ie1-1, ' k=', kb1, '..', ke1-1
  endif

  converged = .false.
  outer: do it=1, Niter_Couple_Outer
    fp_Tw_old = fp_Tw
!   adaptive gas-chunk length: keep the full Kstep_Couple_Comp for the first
!   Niter_Couple_Warm outer iterations (warm start), then halve it each outer
!   iteration down to Kstep_Couple_Min for the final refinement stage.
    nstep = Kstep_Couple_Comp
    if(it > Niter_Couple_Warm) then
      nhalve = it - Niter_Couple_Warm
      do i=1, nhalve
        nstep = max(Kstep_Couple_Min, nstep/2)
      enddo
    endif
    nstep = max(1, nstep)
    if(my_id == 0) print*, ' outer iter', it, ': gas chunk steps =', nstep

!   ==================== (1) GAS CHUNK ===========================
    do step=1, nstep
!     impose the isothermal wall on the code-19 face of the gas block
      call fill_gas_wall_ghost(NVAR1)
!     advance compressible blocks only
      call comput_Sfac(Sfac,Sfac1)
      call Set_Un(nMesh)
      do mBlock=1, Mesh(nMesh)%Num_Block
        B => Mesh(nMesh)%Block(mBlock)
        if(B%Block_type == BLOCK_FLUID) call solver_one_block(nMesh, mBlock, Sfac, Sfac1)
      enddo
      if(IFLAG_LIMIT_FLOW == 1) call limit_flow(nMesh)
      call Boundary_condition_onemesh(nMesh)
      call update_buffer_onemesh(nMesh)
      call update_Ts_buffer_onemesh(nMesh)
      Mesh(nMesh)%tt = Mesh(nMesh)%tt + dt_global
      Mesh(nMesh)%Kstep = Mesh(nMesh)%Kstep + 1
      if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_show) == 0) then
        call comput_force
        call output_Res(nMesh)
      endif
      if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_save) == 0) then
        call output_flow
        call output_Ts
        call output_vtk
      endif
      !  restart file + node-centred SI flow field: MUST be called by ALL
      !  ranks (the routine does its own MPI + Kstep-based due check)
      call restart_step_output(nMesh)
    enddo
    if(my_id == 0) print*, ' gas chunk done, Kstep=', Mesh(nMesh)%Kstep, ' tt=', Mesh(nMesh)%tt

!   refresh the wall ghost once more against the final gas interior so the
!   extracted q_w is consistent with this chunk's isothermal wall state
    call fill_gas_wall_ghost(NVAR1)

!   ==================== (2) EXTRACT q_w, p_w =========================
    qwmax = 0.d0
    if(pair42) then
      do k = kb1, ke1-1
      do i = ib1, ie1-1
        gi = jb + (ie1-1) - i
        r1  = max(Bf%U(1,ic_h,gi,k), 1.d-20)
        u1s = Bf%U(2,ic_h,gi,k)/r1
        v1s = Bf%U(3,ic_h,gi,k)/r1
        p1n = (Bf%U(5,ic_h,gi,k) - 0.5d0*r1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
        T1_nd = gamma*Ma*Ma*p1n/max(r1,1.d-20)
        T1_K  = T1_nd*T_inf
        mu_SI = MU_SI0*sqrt((T1_K/T_SI0)**3)*(T_SI0+S_SI)/(T1_K+S_SI)
        k_gas  = mu_SI*cp_ref/max(PrL,1.d-30)
        dxp = sqrt( (Bf%xc(ic_h,gi,k)-Bf%xc(ic_h+1,gi,k))**2 &
                  + (Bf%yc(ic_h,gi,k)-Bf%yc(ic_h+1,gi,k))**2 &
                  + (Bf%zc(ic_h,gi,k)-Bf%zc(ic_h+1,gi,k))**2 )*Lscale
        dxp = max(dxp, 1.d-20)
        Utmp(1:NVAR1) = Bf%U(1:NVAR1,ic_h+1,gi,k)
!       wall heat flux into the wall (>0): q = k_gas*(T_int - T_ghost)/dx
        fp_qw(i,k) = k_gas*(T1_nd - T_nd_from_U(Utmp, NVAR1, Ma, gamma))/dxp * T_inf
        fp_pw(i,k) = p1n*RHO_REF*U_ref*U_ref
        qwmax = max(qwmax, abs(fp_qw(i,k)))
      enddo; enddo
    else
      do k=kb, ke-1
      do i=ib, ie-1
        r1  = max(Bf%U(1,i,jc_h,k), 1.d-20)
        u1s = Bf%U(2,i,jc_h,k)/r1
        v1s = Bf%U(3,i,jc_h,k)/r1
        p1n = (Bf%U(5,i,jc_h,k) - 0.5d0*r1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
        T1_nd = gamma*Ma*Ma*p1n/max(r1,1.d-20)
        T1_K  = T1_nd*T_inf
        mu_SI = MU_SI0*sqrt((T1_K/T_SI0)**3)*(T_SI0+S_SI)/(T1_K+S_SI)
        k_gas  = mu_SI*cp_ref/max(PrL,1.d-30)
        dxp = max((Bf%yc(i,jc_h,k)-Bf%yc(i,jc_h-1,k))*Lscale, 1.d-20)
        Utmp(1:NVAR1) = Bf%U(1:NVAR1,i,jc_h-1,k)
        fp_qw(i,k) = k_gas*(T1_nd - T_nd_from_U(Utmp, NVAR1, Ma, gamma))/dxp * T_inf
        fp_pw(i,k) = p1n*RHO_REF*U_ref*U_ref
        qwmax = max(qwmax, abs(fp_qw(i,k)))
      enddo; enddo
    endif
    if(my_id == 0) print*, ' gas heat-flux max|q_w|=', qwmax, ' W/m2'

!   ==================== (3) POROUS CHUNK =========================
    do pc=1, Porous_Chunk_Iter
!     porous interface face boundary: outlet pressure fp_pw + heat-flux Ts BC
      do k=kb1, ke1-1
      do i=ib1, ie1-1
        if(pair42) then
          do n1=0, LAP-1
            jg = jg_l - n1
            Bp%U(1,i,jg,k) = Bp%U(1,i,jc_l,k)
            Bp%U(2,i,jg,k) = Bp%U(2,i,jc_l,k)
            Bp%U(3,i,jg,k) = Bp%U(3,i,jc_l,k)
            Bp%U(4,i,jg,k) = Bp%U(4,i,jc_l,k)
            Bp%U(5,i,jg,k) = Bp%U(5,i,jc_l,k)
            Bp%p(i,jg,k)   = fp_pw(i,k)
          enddo
          call set_porous_Ts_flux(Bp, i, jc_l, jg_l, k, -1, fp_qw(i,k))
        else
          do n1=0, LAP-1
            jg = jg_l + n1
            Bp%U(1,i,jg,k) = Bp%U(1,i,jc_l,k)
            Bp%U(2,i,jg,k) = Bp%U(2,i,jc_l,k)
            Bp%U(3,i,jg,k) = Bp%U(3,i,jc_l,k)
            Bp%U(4,i,jg,k) = Bp%U(4,i,jc_l,k)
            Bp%U(5,i,jg,k) = Bp%U(5,i,jc_l,k)
            Bp%p(i,jg,k)   = fp_pw(i,k)
          enddo
          call set_porous_Ts_flux(Bp, i, jc_l, jg_l, k, +1, fp_qw(i,k))
        endif
      enddo; enddo
      do mBlock=1, Mesh(nMesh)%Num_Block
        B => Mesh(nMesh)%Block(mBlock)
        if(B%Block_type == BLOCK_POROUS) call solver_one_block(nMesh, mBlock, Sfac, Sfac1)
      enddo
    enddo

!   ==================== (4) RETURN T_w from porous interface face ===
    do k=kb1, ke1-1
    do i=ib1, ie1-1
      if(pair42) then
        fp_Tw(i,k) = 0.5d0*(Bp%Ts(i,jc_l,k) + Bp%Ts(i,jg_l,k))
      else
        fp_Tw(i,k) = 0.5d0*(Bp%Ts(i,jc_l,k) + Bp%Ts(i,jc_l+1,k))
      endif
    enddo; enddo
    !   coolant exit (transpiration) velocity at the porous interface face
        vcmax = 0.d0
        do k=kb1, ke1-1
        do i=ib1, ie1-1
          vcmax = max(vcmax, abs(Bp%U(3,i,jc_l,k)))
        enddo; enddo
    twmax = 0.d0
    do k=kb1, ke1-1
    do i=ib1, ie1-1
      twmax = max(twmax, abs(fp_Tw(i,k)-fp_Tw_old(i,k)))
    enddo; enddo
    if(my_id == 0) then
      print*, ' porous chunk done, outer iter', it, ' max|dT_w|=', twmax, ' K', &
              '  coolant exit |v|_max=', vcmax, ' m/s', &
              '  T_w range [', minval(fp_Tw(ib1:ie1-1,kb1:ke1-1)), ',', &
              maxval(fp_Tw(ib1:ie1-1,kb1:ke1-1)), ']'
      open(203, file='iface_couple.dat', status='replace')
      write(203,'(A)') '# x_w(m)  T_w(K)  q_w(W/m2)  p_w(Pa)'
      do k=kb1, ke1-1
      do i=ib1, ie1-1
        if(pair42) then
          gi = jb + (ie1-1) - i
          xw = 0.5d0*(Bf%x(ie,gi,k)+Bf%x(ie,gi+1,k))
        else
          xw = 0.5d0*(Bf%x(i,jb,k)+Bf%x(i+1,jb,k))
        endif
        write(203,'(4ES16.7)') xw*Lscale, fp_Tw(i,k), fp_qw(i,k), fp_pw(i,k)
      enddo; enddo
      close(203)
    endif
!    登记本轮耦合状态（pair 19；mode=1 交错）
    call set_couple_state(1, merge(1,0,twmax < Tol_Couple_Tw), 19, it, &
                          twmax, Tol_Couple_Tw)
    if(it >= 2 .and. twmax < Tol_Couple_Tw) then
      converged = .true.
      if(my_id == 0) print*, ' Staggered coupling converged at outer iter', it, &
                             ' (max|dT_w| <', Tol_Couple_Tw, ')'
      exit outer
    endif
  enddo outer

  if(.not. converged) then
    if(my_id == 0) print*, ' Staggered coupling reached Niter_Couple_Outer=', &
                           Niter_Couple_Outer, ' (not fully converged)'
  endif

  call output_flow
  call output_Ts
  call output_vtk
  call write_restart
  contains

!   Impose the isothermal wall state on the gas block code-19 face cells.
!   pair42=.false. : j- face, interior row j=jb, ghost rows jb-1, jb-2, ...
!   pair42=.true.  : i+ face, interior column i=ie-1, ghost columns ie, ie+1,...
!   T_w = fp_Tw(porous cell index, k) [K].
   subroutine fill_gas_wall_ghost(nv)
     implicit none
     integer:: nv
     integer:: i2, k2, n1
     real(PRE_EC):: mu1c, dwc
     if(pair42) then
       do k2 = kb, ke-1
       do i2 = jb, je-1
         if(If_viscous == 1) then
           mu1c = Bf%mu(ic_h,i2,k2); dwc = Bf%dw(ic_h,i2,k2)
         else
           mu1c = 1.d0/Re; dwc = 0.d0
         endif
         gi = jb + (ie1-1) - i2
         if(LS_Inlet_Type == 1 .and. Bp%U(3,gi,jc_l,k2) > 0.d0) then
         ! transpiration (blowing) wall: the coolant leaves the porous block
         ! in +y with the velocity Bp%U(3,gi,jc_l,k2) [m/s] and enters the gas
           call wall_bound_blowing(nv, Bf%U(:,ic_h,i2,k2), Bf%U(:,ic_h+1,i2,k2), &
              Ma, gamma, fp_Tw(gi,k2)/T_inf, Bp%U(3,gi,jc_l,k2), U_ref)
         else
         call wall_bound_with_Tw(nv, Bf%U(:,ic_h,i2,k2), Bf%U(:,ic_h+1,i2,k2), &
              Ma, gamma, fp_Tw(gi,k2)/T_inf, mu1c, dwc, Re)
         endif
         do n1=1, LAP-1
           Bf%U(:,ic_h+1+n1,i2,k2) = Bf%U(:,ic_h+1,i2,k2)
         enddo
       enddo; enddo
     else
       do k2 = kb, ke-1
       do i2 = ib, ie-1
         if(If_viscous == 1) then
           mu1c = Bf%mu(i2,jc_h,k2); dwc = Bf%dw(i2,jc_h,k2)
         else
           mu1c = 1.d0/Re; dwc = 0.d0
         endif
         call wall_bound_with_Tw(nv, Bf%U(:,i2,jc_h,k2), Bf%U(:,i2,jc_h-1,k2), &
              Ma, gamma, fp_Tw(i2,k2)/T_inf, mu1c, dwc, Re)
         do n1=2, LAP
           Bf%U(:,i2,jc_h-n1,k2) = Bf%U(:,i2,jc_h-1,k2)   ! deeper ghosts
         enddo
       enddo; enddo
     endif
   end subroutine fill_gas_wall_ghost
!   Isothermal wall with normal mass injection (blowing / transpiration).
!   Tw is non-dimensional (T/T_inf), Vn the wall-normal velocity [m/s] in +y,
!   Uref the velocity scale.  Tangential velocity mirrored (wall at rest),
!   pressure copied, T = 2*Tw - T_int (clamped at 0.5*T_int).
   subroutine wall_bound_blowing(nvar, U1, Ug1, Ma1, gam, Tw, Vn, Uref)
     implicit none
     integer:: nvar
     real(PRE_EC):: U1(nvar), Ug1(nvar), Ma1, gam, Tw, Vn, Uref
     real(PRE_EC):: d1, uu1, v1, w1, p1, T1, p2, T2, u2, v2, w2, d2
     d1 = max(U1(1),1.d-20)
     uu1 = U1(2)/d1; v1 = U1(3)/d1; w1 = U1(4)/d1
     p1 = (U1(5) - 0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gam-1.d0)
     T1 = gam*Ma1*Ma1*p1/d1
     p2 = p1
     T2 = 2.d0*Tw - T1
     if(T2 .lt. 0.5d0*T1) T2 = 0.5d0*T1
     d2 = gam*Ma1*Ma1*p2/T2
     u2 = -uu1
     v2 = 2.d0*(Vn/max(Uref,1.d-30)) - v1
     w2 = -w1
     Ug1(1)=d2
     Ug1(2)=d2*u2
     Ug1(3)=d2*v2
     Ug1(4)=d2*w2
     Ug1(5)=p2/(gam-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
     if(nvar .ge. 6) Ug1(6)=U1(6)
     if(nvar .ge. 7) Ug1(7)=U1(7)
   end subroutine wall_bound_blowing

!   Non-dimensional temperature T/T_inf of a compressible state U(:)
   real(PRE_EC) function T_nd_from_U(U1, nv, Ma1, gam)
     implicit none
     real(PRE_EC):: U1(nv)
     integer:: nv
     real(PRE_EC):: Ma1, gam, d1, uu, vv, ww, p1
     d1 = max(U1(1), 1.d-20)
     uu = U1(2)/d1; vv = U1(3)/d1; ww = U1(4)/d1
     p1 = (U1(5) - 0.5d0*d1*(uu*uu+vv*vv+ww*ww))*(gam-1.d0)
     T_nd_from_U = gam*Ma1*Ma1*p1/max(d1,1.d-20)
   end function T_nd_from_U

!   Solid-frame Ts ghost on the porous interface face from the incoming heat
!   flux qw>0 (W/m2, heat entering the wall): Ts_g = Ts_i + qw*dx/ks_eff.
!   sgn = +1 : ghost rows grow with +j (j+ face);  sgn = -1 : with -j (j- face)
   subroutine set_porous_Ts_flux(B, i, jint, jgh, k, sgn, qw)
     implicit none
     Type (Block_TYPE),pointer:: B
     integer:: i, jint, jgh, k, sgn
     real(PRE_EC):: qw, dx_g, ks_eff
     integer:: n1, jg2
     ks_eff = max((1.d0 - B%porous_eps)*B%solid_k, 1.d-30)
     do n1 = 0, LAP-1
       jg2 = jgh + sgn*n1
       dx_g = sqrt( (B%xc(i,jint,k)-B%xc(i,jg2,k))**2 &
                  + (B%yc(i,jint,k)-B%yc(i,jg2,k))**2 &
                  + (B%zc(i,jint,k)-B%zc(i,jg2,k))**2 ) * Lscale
       B%Ts(i,jg2,k) = B%Ts(i,jint,k) + qw*dx_g/ks_eff
     enddo
   end subroutine set_porous_Ts_flux
  end subroutine run_staggered_fluid_porous

!==============================================================================
! Staggered segmented coupling for FLUID<->SOLID / LOWSPEED(AC)<->SOLID
! (interface codes 11 / 13), conjugate-heat-transfer (CHT) variant.
!
! Same outer driver as interface 19: warm-up fluid/LS chunks of
! Kstep_Couple_Comp steps (halved after Niter_Couple_Warm) alternating with
! fully-converged solid chunks.  Interface ghosts are refreshed by the standard
! per-step conjugate couple (couple_compressible_fluid_solid_face for 11,
! couple_lowspeed_fluid_solid_face for 13) at every chunk step; the solid chunk
! is a full Gauss-Seidel steady solve under the current interface ghost.
! Convergence is monitored on the interface temperature T_w recovered from the
! solid ghost (T_i = (Ts_js + Ts_jg)/2).
!
! paircode = 11 : BLOCK_FLUID(compressible) <-> BLOCK_SOLID
!          = 13 : BLOCK_LOWSPEED(AC)       <-> BLOCK_SOLID
!==============================================================================
  subroutine run_staggered_fluid_solid(nMesh, paircode)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh, paircode
   integer:: ftype, s1, f1, mBlock, mBlock2, ksub, nb, it, step, i1, i, k
   integer:: nstep, nhalve
   integer:: face_s, face1, ib,ie,jb,je,kb,ke, ib1,ie1,jb1,je1,kb1,ke1
   integer:: js, jg, nf, nk
   real(PRE_EC):: Sfac, Sfac1, twmax, T_i, qw, dx_s, tmpu
   logical:: found, converged, wfmode
   Type (Block_TYPE),pointer:: Bs, Bf, B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2

   ftype = BLOCK_FLUID
   if(paircode .eq. 13) ftype = BLOCK_LOWSPEED
   wfmode = (Iflag_Couple_WallFlux == 1)   ! CHT split: isothermal-Tw wall / q_w-flux solid

   s1 = 0; f1 = 0
   do mBlock = 1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type == BLOCK_SOLID .and. s1 == 0) s1 = mBlock
     if(B%Block_type == ftype .and. f1 == 0) f1 = mBlock
   enddo
   if(s1 == 0 .or. f1 == 0) then
     print*, 'run_staggered_fluid_solid: need one BLOCK_SOLID and one fluid block'// &
             ' (type ', ftype, '), got', s1, f1
     return
   endif
   Bs => Mesh(nMesh)%Block(s1)

!  find the solid-side interface subface to the fluid block.  Accept both the
!  explicit interface markers (bc<0 or 11..19) and the link-style entries used
!  by some grids (physical wall code with neighbour fields nb1/face1 set).
   found = .false.
   do ksub = 1, Bs%subface
     Bc2 => Bs%bc_msg2(ksub)
     if(.not. is_interface_bc(Bc2%bc)) then
       if(Bc2%nb1 <= 0 .or. Bc2%face1 <= 0) cycle
     endif
     nb = Bc2%nb1
     if(nb <= 0) cycle
     mBlock2 = 0
     do mBlock = 1, Mesh(nMesh)%Num_Block
       if(Mesh(nMesh)%Block(mBlock)%Block_no == nb) then
         Bn => Mesh(nMesh)%Block(mBlock); mBlock2 = mBlock; exit
       endif
     enddo
     if(mBlock2 == 0) cycle
     if(Bn%Block_type /= ftype) cycle
     face_s = Bc2%face; face1 = Bc2%face1
     if(paircode .eq. 11) then
       if(face_s /= 2 .or. .not.(face1 == 1 .or. face1 == 4)) then
         if(my_id == 0) print*, 'run_staggered_fluid_solid(11): supports solid j- (2) vs'// &
             ' comp i-/i+ (1/4) only; got', face_s, face1
         return
       endif
     else
       if(.not. (face_s == 2 .or. face_s == 5)) then
         if(my_id == 0) print*, 'run_staggered_lowspeed_solid(13): supports solid j-/'// &
             'j+ (2/5) only; got', face_s
         return
       endif
       if(.not. (face1 == 2 .or. face1 == 5)) then
         if(my_id == 0) print*, 'run_staggered_lowspeed_solid(13): fluid face must be'// &
             ' j- or j+ (2/5), got', face1
         return
       endif
     endif
     ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
     ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1
     Bf => Mesh(nMesh)%Block(mBlock2); f1 = mBlock2
     found = .true.
     exit
   enddo
   if(.not. found) then
     print*, 'run_staggered_fluid_solid: no FLUID/LOWSPEED - SOLID (code ', paircode, &
             ') interface found'
     return
   endif

   nf = ie - ib
   nk = ke - kb
   if(.not. allocated(fp_Tw)) then
     allocate(fp_Tw(1:nf,1:nk), fp_Tw_old(1:nf,1:nk), fp_qw(1:nf,1:nk))
     allocate(fp_pw(1:nf,1:nk), fp_u(1:nf,1:nk), fp_pw_old(1:nf,1:nk), fp_u_old(1:nf,1:nk))
   else
     if(size(fp_Tw,1) /= nf .or. size(fp_Tw,2) /= nk) then
       deallocate(fp_Tw, fp_Tw_old, fp_qw, fp_pw, fp_u, fp_pw_old, fp_u_old)
       allocate(fp_Tw(1:nf,1:nk), fp_Tw_old(1:nf,1:nk), fp_qw(1:nf,1:nk))
       allocate(fp_pw(1:nf,1:nk), fp_u(1:nf,1:nk), fp_pw_old(1:nf,1:nk), fp_u_old(1:nf,1:nk))
     endif
   endif
   fp_Tw = Twall_Couple_Init
   fp_Tw_old = fp_Tw; fp_qw = 0.d0; fp_pw = 0.d0; fp_u = 0.d0
   fp_pw_old = 0.d0; fp_u_old = 0.d0
!  重启：若 restart 文件的 trailer 里正好带着本界面的量（同 pair、同尺寸），
!  直接恢复，不再用 control.ec 初值重播（跨重启零跳变）。Iflag_Couple_Restart=-1
!  表示忽略重启里的耦合状态（=旧行为）。
   if(Iflag_Couple_Restart >= 0 .and. Couple_State_Found == 1 .and. &
      Couple_State_Pair == paircode .and. Couple_State_nf == nf .and. &
      Couple_State_nk == nk .and. allocated(Couple_Tw_save)) then
     fp_Tw = Couple_Tw_save;  fp_qw = Couple_qw_save
     fp_u  = Couple_u_save;   fp_pw = Couple_pw_save
     fp_Tw_old = fp_Tw; fp_pw_old = fp_pw; fp_u_old = fp_u
     if(my_id == 0) print*, ' restart: interface T_w/q_w restored from the', &
            ' restart file (no re-seeding; last metric=', Couple_State_Twmax, ')'
   endif

   !  restart file: keep the saved counters instead of resetting
   if(restart_found == 1) then
     Mesh(nMesh)%tt    = restart_tt_saved
     Mesh(nMesh)%Kstep = restart_Kstep_saved
     if(my_id == 0) print*, ' restart: continue staggered run at Kstep=', &
         restart_Kstep_saved, ' tt=', restart_tt_saved
   !    NOTE: the per-face interface quantities (fp_Tw/fp_u) are re-seeded
   !    from the control.ec initial guess; the outer iterations re-converge
   !    them within a few rounds (state itself is fully restored).
   else
   Mesh(nMesh)%tt = 0.d0
   Mesh(nMesh)%Kstep = 0
   endif
   if(my_id == 0) then
     print*, ' run_staggered_fluid_solid(paircode=', paircode, '): solid block', s1, &
             ' fluid block', f1, ' Kstep_Couple_Comp=', Kstep_Couple_Comp, &
             ' Niter_Couple_Outer=', Niter_Couple_Outer, ' Twall_Couple_Init=', Twall_Couple_Init
   endif

   converged = .false.
   outer: do it = 1, Niter_Couple_Outer
     fp_Tw_old = fp_Tw
!    adaptive chunk length: full Kstep_Couple_Comp for the first
!    Niter_Couple_Warm outer iterations, then halve down to Kstep_Couple_Min
     nstep = Kstep_Couple_Comp
     if(it > Niter_Couple_Warm) then
       nhalve = it - Niter_Couple_Warm
       do i1 = 1, nhalve
         nstep = max(Kstep_Couple_Min, nstep/2)
       enddo
     endif
     nstep = max(1, nstep)
!    code 13: the low-speed flow side is solved by the AC solver "to
!    convergence" per outer round (one solver_one_block call); the chunk-step
!    schedule applies to the compressible case (11) only.
     if(paircode .eq. 13) nstep = 1
     if(my_id == 0) print*, ' outer iter', it, ': flow chunk steps =', nstep

!   ==================== (1) FLUID / LS CHUNK ========================
      if(.not. wfmode) call couple_fluid_solid_interfaces(nMesh)   ! seed interface ghosts before the chunk
     do step = 1, nstep
        if(wfmode) call stagger_fill_flow_wall
       call comput_Sfac(Sfac, Sfac1)
       call Set_Un(nMesh)
       do mBlock = 1, Mesh(nMesh)%Num_Block
         B => Mesh(nMesh)%Block(mBlock)
         if(B%Block_type == ftype) call solver_one_block(nMesh, mBlock, Sfac, Sfac1)
       enddo
        if(.not. wfmode) call couple_fluid_solid_interfaces(nMesh)
       call couple_solid_solid_interfaces(nMesh)
       call couple_lowspeed_porous_interfaces(nMesh)
       if(IFLAG_LIMIT_FLOW == 1) call limit_flow(nMesh)
       call Boundary_condition_onemesh(nMesh)
       call update_buffer_onemesh(nMesh)
       call update_Ts_buffer_onemesh(nMesh)
       Mesh(nMesh)%tt = Mesh(nMesh)%tt + dt_global
       Mesh(nMesh)%Kstep = Mesh(nMesh)%Kstep + 1
       if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_show) == 0) then
         call comput_force
         call output_Res(nMesh)
       endif
       if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_save) == 0) then
         call output_flow
         call output_Ts
         call output_vtk
       endif
       !  restart file + node-centred SI flow field: MUST be called by ALL
       !  ranks (the routine does its own MPI + Kstep-based due check)
       call restart_step_output(nMesh)
     enddo
      if(wfmode) then
        call stagger_fill_flow_wall
        call stagger_extract_flow_qw
      endif

     if(my_id == 0) print*, ' flow chunk done, Kstep=', Mesh(nMesh)%Kstep, &
         ' tt=', Mesh(nMesh)%tt

!   ==================== (2) SOLID CHUNK (steady GS) =================
      if(wfmode) then
        stg_ow_block = s1; stg_ow_face = face_s
        stg_ow_ib=ib; stg_ow_ie=ie; stg_ow_jb=jb; stg_ow_je=je; stg_ow_kb=kb; stg_ow_ke=ke
      endif
!    interface ghost is current from the last couple of the flow chunk
     do mBlock = 1, Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
       if(B%Block_type == BLOCK_SOLID) call solid_solver_one_block(nMesh, mBlock)
     enddo
!    refresh ghosts with the updated solid and store the interface T_w
      if(wfmode) then
        stg_ow_block = 0
      else
        call couple_fluid_solid_interfaces(nMesh)
      endif

     if(face_s == 2) then
       js = Bc2%jb
       jg = Bc2%jb - 1
     else
       js = Bc2%je - 1
       jg = Bc2%je
     endif
     twmax = 0.d0
     do k = kb, ke-1
     do i = ib, ie-1
       T_i = 0.5d0*(Bs%Ts(i,js,k) + Bs%Ts(i,jg,k))
       if(face_s == 2) then
         dx_s = Bs%y(i,jb,k) - Bs%yc(i,js,k)
       else
         dx_s = Bs%yc(i,js,k) - Bs%y(i,je,k)
       endif
       dx_s = max(abs(dx_s), 1.d-20)*Lscale
       qw   = Bs%solid_k*(T_i - Bs%Ts(i,js,k))/dx_s
       fp_Tw(i-ib+1,k-kb+1) = T_i
       if(.not. wfmode) fp_qw(i-ib+1,k-kb+1) = qw
       twmax = max(twmax, abs(T_i - fp_Tw_old(i-ib+1,k-kb+1)))
     enddo; enddo
     if(my_id == 0) then
       print*, ' solid chunk done, outer iter', it, ' max|dT_w|=', twmax, ' K', &
               '  T_w range [', minval(fp_Tw), ',', maxval(fp_Tw), ']  max|q_w|=', maxval(abs(fp_qw))
       open(203, file='iface_couple.dat', status='replace')
       write(203,'(A)') '# x_w(m)  T_w(K)  q_w(W/m2)  p_w(Pa)'
       do k = kb, ke-1
       do i = ib, ie-1
         tmpu = 0.5d0*(Bs%x(i,js,k)+Bs%x(i+1,js,k))
         write(203,'(4ES16.7)') tmpu*Lscale, fp_Tw(i-ib+1,k-kb+1), &
                                fp_qw(i-ib+1,k-kb+1), 0.d0
       enddo; enddo
       close(203)
     endif
!    把本轮耦合状态登记进重启文件的 trailer（mode=1 交错；conv=1 表示已满足判据）
     call set_couple_state(1, merge(1,0,twmax < Tol_Couple_Tw), paircode, it, &
                           twmax, Tol_Couple_Tw)
     if(it >= 2 .and. twmax < Tol_Couple_Tw) then
       converged = .true.
       if(my_id == 0) print*, ' Staggered CHT coupling converged at outer iter', it, &
                              ' (max|dT_w| <', Tol_Couple_Tw, ')'
       exit outer
     endif
   enddo outer

   if(.not. converged) then
     if(my_id == 0) print*, ' Staggered CHT coupling reached Niter_Couple_Outer=', &
                            Niter_Couple_Outer, ' (not fully converged)'
   endif

   call output_flow
   call output_Ts
   call output_vtk

  call write_restart
  contains

!   Fill the flow-side interface ghost as an isothermal wall at fp_Tw.
!   code 13: BLOCK_LOWSPEED (mirror velocities, zero-gradient p, U5=2Tw-T);
!   code 11: compressible wall via wall_bound_with_Tw (i- face, tangential j map).
    subroutine stagger_fill_flow_wall
      implicit none
      integer:: i, k, jt, jint, jg, n1
      real(PRE_EC):: Tw
      if(paircode .eq. 13) then
        do k = kb, ke-1
        do i = ib, ie-1
          Tw = fp_Tw(i-ib+1,k-kb+1)
          if(face1 .eq. 2) then
            jint = jb1
          else
            jint = je1-1
          endif
          do n1 = 0, LAP-1
            if(face1 .eq. 2) then
              jg = jb1-1-n1
            else
              jg = je1+n1
            endif
            Bf%U(1,i,jg,k) = LS_rho
            Bf%U(2,i,jg,k) = -Bf%U(2,i,jint,k)
            Bf%U(3,i,jg,k) = -Bf%U(3,i,jint,k)
            Bf%U(4,i,jg,k) = -Bf%U(4,i,jint,k)
            Bf%U(5,i,jg,k) = 2.d0*Tw - Bf%U(5,i,jint,k)
            Bf%p(i,jg,k)   = Bf%p(i,jint,k)
          enddo
        enddo; enddo
      else
        do k = kb, ke-1
        do i = ib, ie-1
          jt = jb1 + (i - ib)
          Tw = fp_Tw(i-ib+1,k-kb+1)
          call wall_bound_with_Tw(Mesh(nMesh)%NVAR, Bf%U(:,ib1,jt,k), Bf%U(:,ib1-1,jt,k), &
               Ma, gamma, Tw/T_inf, Bf%mu(ib1,jt,k), Bf%dw(ib1,jt,k), Re)
          do n1 = 1, LAP-1
            Bf%U(:,ib1-1-n1,jt,k) = Bf%U(:,ib1-1,jt,k)
          enddo
        enddo; enddo
      endif
    end subroutine stagger_fill_flow_wall

!   Extract the interface heat flux into the solid (W/m2, >0 into wall) from the
!   flow-side isothermal-wall ghost after the flow chunk.
    subroutine stagger_extract_flow_qw
      implicit none
      integer:: i, k, jt, jint, jg
      real(PRE_EC):: Tint, Tgh, dxg, q
      real(PRE_EC),parameter:: R_AIR=287.d0, MU_SI0=1.716d-5, T_SI0=273.15d0, S_SI=110.4d0
      real(PRE_EC):: d1, u1s, v1s, p1n, Tnd, Tnd_g, T1K, muSI, kgas, cp_ref
      if(paircode .eq. 13) then
        do k = kb, ke-1
        do i = ib, ie-1
          if(face1 .eq. 2) then
            jint = jb1
          else
            jint = je1-1
          endif
          if(face1 .eq. 2) then
            jg = jb1-1
          else
            jg = je1
          endif
          Tint = Bf%U(5,i,jint,k)
          Tgh  = Bf%U(5,i,jg,k)
          dxg  = sqrt( (Bf%xc(i,jint,k)-Bf%xc(i,jg,k))**2 &
                     + (Bf%yc(i,jint,k)-Bf%yc(i,jg,k))**2 &
                     + (Bf%zc(i,jint,k)-Bf%zc(i,jg,k))**2 )*Lscale
          fp_qw(i-ib+1,k-kb+1) = LS_k*(Tint-Tgh)/max(dxg,1.d-30)
        enddo; enddo
      else
        do k = kb, ke-1
        do i = ib, ie-1
          jt = jb1 + (i - ib)
          d1  = max(Bf%U(1,ib1,jt,k),1.d-20)
          u1s = Bf%U(2,ib1,jt,k)/d1; v1s = Bf%U(3,ib1,jt,k)/d1
          p1n = (Bf%U(5,ib1,jt,k)-0.5d0*d1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          Tnd = gamma*Ma*Ma*p1n/max(d1,1.d-20)
          d1  = max(Bf%U(1,ib1-1,jt,k),1.d-20)
          u1s = Bf%U(2,ib1-1,jt,k)/d1; v1s = Bf%U(3,ib1-1,jt,k)/d1
          p1n = (Bf%U(5,ib1-1,jt,k)-0.5d0*d1*(u1s*u1s+v1s*v1s))*(gamma-1.d0)
          Tnd_g = gamma*Ma*Ma*p1n/max(d1,1.d-20)
          dxg  = sqrt( (Bf%xc(ib1,jt,k)-Bf%xc(ib1-1,jt,k))**2 &
                     + (Bf%yc(ib1,jt,k)-Bf%yc(ib1-1,jt,k))**2 &
                     + (Bf%zc(ib1,jt,k)-Bf%zc(ib1-1,jt,k))**2 )*Lscale
          T1K  = Tnd*T_inf
          muSI = MU_SI0*sqrt((T1K/T_SI0)**3)*(T_SI0+S_SI)/(T1K+S_SI)
          cp_ref = gamma*R_AIR/(gamma-1.d0)
          kgas = muSI*cp_ref/max(PrL,1.d-30)
          fp_qw(i-ib+1,k-kb+1) = kgas*(Tnd-Tnd_g)/max(dxg,1.d-30)*T_inf
        enddo; enddo
      endif
    end subroutine stagger_extract_flow_qw

  end subroutine run_staggered_fluid_solid
!==============================================================================
! Staggered segmented coupling for FLUID(compressible, high speed) <->
! LOWSPEED fluid (interface code 12), two-fluid matching block-Gauss-Seidel
! variant (chosen by the user).
!
! Each outer round: (1) gas chunk - advance only the BLOCK_FLUID blocks for the
! scheduled number of steps while the low-speed blocks stay frozen; (2) LS
! chunk - advance only the BLOCK_LOWSPEED blocks (AC, LS_Algorithm=3, to
! convergence) while the gas stays frozen.  The standard per-step interface
! couple couple_highlow_fluid_face (unit conversion both ways + LS velocity
! clamping / pressure anchoring) is applied after every solver call, so the
! advancing side sees the frozen partner state on its interface ghost.
! Convergence is monitored on the low-speed top-row interface temperature,
! pressure and velocity between outer rounds.
!==============================================================================
  subroutine run_staggered_highlow(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   integer:: mf, ml, mBlock, mBlock2, ksub, nb, it, step, i1, i, k
   integer:: nstep, nhalve, nf, nk, ia1, ia2, ka1, ka2, gi
   integer:: face_s, face1, ib,ie,jb,je,kb,ke, ib1,ie1,jb1,je1,kb1,ke1, jc_l
   real(PRE_EC):: Sfac, Sfac1, dtw, dpw, duw, xw, uu, vv
   logical:: found, converged, pair42
   Type (Block_TYPE),pointer:: Bf, Bp, B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2

   mf = 0; ml = 0
   do mBlock = 1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type == BLOCK_FLUID .and. mf == 0) mf = mBlock
     if(B%Block_type == BLOCK_LOWSPEED .and. ml == 0) ml = mBlock
   enddo
   if(mf == 0 .or. ml == 0) then
     print*, 'run_staggered_highlow: need one BLOCK_FLUID and one BLOCK_LOWSPEED,'// &
             ' got', mf, ml
     return
   endif
   Bf => Mesh(nMesh)%Block(mf)

!  find the code-12 subface on a BLOCK_FLUID block (couple_highlow is applied
!  once from the compressible side); support the j- vs j+ conformal pair.
   found = .false.
   do mBlock = 1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type /= BLOCK_FLUID) cycle
     do ksub = 1, B%subface
       Bc2 => B%bc_msg2(ksub)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1
       if(nb <= 0) cycle
       mBlock2 = 0
       do mBlock2 = 1, Mesh(nMesh)%Num_Block
         if(Mesh(nMesh)%Block(mBlock2)%Block_no == nb) exit
       enddo
       if(mBlock2 > Mesh(nMesh)%Num_Block) cycle
       Bn => Mesh(nMesh)%Block(mBlock2)
       if(Bn%Block_type /= BLOCK_LOWSPEED) cycle
       face_s = Bc2%face; face1 = Bc2%face1
       if(.not. ((face_s == 2 .and. face1 == 5) .or. &
                 (face_s == 4 .and. face1 == 2))) then
         if(my_id == 0) print*, 'run_staggered_highlow: supports comp j- (2) vs LS j+ (5)'// &
             ' (2/5) or comp i+/LS j- (4/2) only; got', face_s, face1
         return
       endif
       ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
       ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1
       Bf => Mesh(nMesh)%Block(mBlock)
       Bp => Mesh(nMesh)%Block(mBlock2)
       mf = mBlock; ml = mBlock2
       found = .true.
       exit
     enddo
     if(found) exit
   enddo
   if(.not. found) then
     print*, 'run_staggered_highlow: no FLUID-LOWSPEED (code 12) interface found'
     return
   endif
   pair42 = (face_s == 4 .and. face1 == 2)
   if(pair42) then
   ! cases/fluid_solid topology: gas i+ face <-> low-speed j- face.  The
   ! low-speed interface cells are i=ib1..ie1-1 (tangential index) at the single
   ! interior row j=jb1; the coolant enters through the far j=max face.
   jc_l = jb1
   ia1 = ib1; ia2 = ie1 - 1
   ka1 = kb1; ka2 = ke1 - 1
   else
   jc_l = je1 - 1             ! low-speed first interior row below its j+ face
   ia1 = ib;  ia2 = ie - 1
   ka1 = kb;  ka2 = ke - 1
   endif

   nf = ia2 - ia1 + 1
   nk = ka2 - ka1 + 1
   if(.not. allocated(fp_Tw)) then
     allocate(fp_Tw(1:nf,1:nk), fp_Tw_old(1:nf,1:nk), fp_qw(1:nf,1:nk))
     allocate(fp_pw(1:nf,1:nk), fp_u(1:nf,1:nk), fp_pw_old(1:nf,1:nk), fp_u_old(1:nf,1:nk))
   else
     if(size(fp_Tw,1) /= nf .or. size(fp_Tw,2) /= nk) then
       deallocate(fp_Tw, fp_Tw_old, fp_qw, fp_pw, fp_u, fp_pw_old, fp_u_old)
       allocate(fp_Tw(1:nf,1:nk), fp_Tw_old(1:nf,1:nk), fp_qw(1:nf,1:nk))
       allocate(fp_pw(1:nf,1:nk), fp_u(1:nf,1:nk), fp_pw_old(1:nf,1:nk), fp_u_old(1:nf,1:nk))
     endif
   endif
   fp_Tw = 0.d0; fp_pw = 0.d0; fp_u = 0.d0; fp_Tw_old = 0.d0
   fp_pw_old = 0.d0; fp_u_old = 0.d0; fp_qw = 0.d0
!  重启：trailer 里若带着本界面（pair 12）的量，直接恢复（零跳变）；
!  Iflag_Couple_Restart=-1 时忽略（=旧行为）。
   if(Iflag_Couple_Restart >= 0 .and. Couple_State_Found == 1 .and. &
      Couple_State_Pair == 12 .and. Couple_State_nf == nf .and. &
      Couple_State_nk == nk .and. allocated(Couple_Tw_save)) then
     fp_Tw = Couple_Tw_save; fp_pw = Couple_pw_save; fp_u = Couple_u_save
     fp_qw = Couple_qw_save
     fp_Tw_old = fp_Tw; fp_pw_old = fp_pw; fp_u_old = fp_u
     if(my_id == 0) print*, ' restart: highlow-interface T_w/p_w/u restored', &
            ' from the restart file (no re-seeding; last metric=', Couple_State_Twmax, ')'
   endif

   !  restart file: keep the saved counters instead of resetting
   if(restart_found == 1) then
     Mesh(nMesh)%tt    = restart_tt_saved
     Mesh(nMesh)%Kstep = restart_Kstep_saved
     if(my_id == 0) print*, ' restart: continue staggered run at Kstep=', &
         restart_Kstep_saved, ' tt=', restart_tt_saved
   !    NOTE: the per-face interface quantities (fp_Tw/fp_u) are re-seeded
   !    from the control.ec initial guess; the outer iterations re-converge
   !    them within a few rounds (state itself is fully restored).
   else
   Mesh(nMesh)%tt = 0.d0
   Mesh(nMesh)%Kstep = 0
   endif
   if(my_id == 0) then
     print*, ' run_staggered_highlow(12): gas block', mf, ' LS block', ml, &
             ' Kstep_Couple_Comp=', Kstep_Couple_Comp, ' Niter_Couple_Outer=', Niter_Couple_Outer
   endif

   converged = .false.
   outer: do it = 1, Niter_Couple_Outer
     fp_Tw_old = fp_Tw; fp_pw_old = fp_pw; fp_u_old = fp_u
     nstep = Kstep_Couple_Comp
     if(it > Niter_Couple_Warm) then
       nhalve = it - Niter_Couple_Warm
       do i1 = 1, nhalve
         nstep = max(Kstep_Couple_Min, nstep/2)
       enddo
     endif
     nstep = max(1, nstep)
     if(my_id == 0) print*, ' outer iter', it, ': gas chunk steps =', nstep

!   ============ (1) GAS CHUNK: compressible only ==================
     do step = 1, nstep
       call comput_Sfac(Sfac, Sfac1)
       call Set_Un(nMesh)
       do mBlock = 1, Mesh(nMesh)%Num_Block
         B => Mesh(nMesh)%Block(mBlock)
         if(B%Block_type == BLOCK_FLUID) call solver_one_block(nMesh, mBlock, Sfac, Sfac1)
       enddo
       call couple_fluid_solid_interfaces(nMesh)
       call couple_solid_solid_interfaces(nMesh)
       call couple_lowspeed_porous_interfaces(nMesh)
       if(IFLAG_LIMIT_FLOW == 1) call limit_flow(nMesh)
       call Boundary_condition_onemesh(nMesh)
       call update_buffer_onemesh(nMesh)
       call update_Ts_buffer_onemesh(nMesh)
       Mesh(nMesh)%tt = Mesh(nMesh)%tt + dt_global
       Mesh(nMesh)%Kstep = Mesh(nMesh)%Kstep + 1
       if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_show) == 0) then
         call comput_force
         call output_Res(nMesh)
       endif
       if(my_id == 0 .and. mod(Mesh(nMesh)%Kstep, Kstep_save) == 0) then
         call output_flow
         call output_Ts
         call output_vtk
       endif
       !  restart file + node-centred SI flow field: MUST be called by ALL
       !  ranks (the routine does its own MPI + Kstep-based due check)
       call restart_step_output(nMesh)
     enddo
     if(my_id == 0) print*, ' gas chunk done, Kstep=', Mesh(nMesh)%Kstep, &
         ' tt=', Mesh(nMesh)%tt

!   ============ (2) LS CHUNK: AC to convergence ====================
     call comput_Sfac(Sfac, Sfac1)
     call Set_Un(nMesh)
     do mBlock = 1, Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
       if(B%Block_type == BLOCK_LOWSPEED) call solver_one_block(nMesh, mBlock, Sfac, Sfac1)
     enddo
     call couple_fluid_solid_interfaces(nMesh)
     call couple_solid_solid_interfaces(nMesh)
     call couple_lowspeed_porous_interfaces(nMesh)
     if(IFLAG_LIMIT_FLOW == 1) call limit_flow(nMesh)
     call Boundary_condition_onemesh(nMesh)
     call update_buffer_onemesh(nMesh)
     call update_Ts_buffer_onemesh(nMesh)
     if(my_id == 0) print*, ' LS (AC) chunk done'

!   ============ (3) RECORD interface metrics & convergence =========
     dtw = 0.d0; dpw = 0.d0; duw = 0.d0
       do k = ka1, ka2
       do i = ia1, ia2
       fp_Tw(i-ia1+1,k-ka1+1) = Bp%U(5,i,jc_l,k)
       fp_pw(i-ia1+1,k-ka1+1) = Bp%p(i,jc_l,k)
       uu = Bp%U(2,i,jc_l,k); vv = Bp%U(3,i,jc_l,k)
       fp_u(i-ia1+1,k-ka1+1) = sqrt(uu*uu + vv*vv)
       dtw = max(dtw, abs(fp_Tw(i-ia1+1,k-ka1+1)-fp_Tw_old(i-ia1+1,k-ka1+1)))
       dpw = max(dpw, abs(fp_pw(i-ia1+1,k-ka1+1)-fp_pw_old(i-ia1+1,k-ka1+1)))
       duw = max(duw, abs(fp_u(i-ia1+1,k-ka1+1)-fp_u_old(i-ia1+1,k-ka1+1)))
       enddo; enddo
     if(my_id == 0) then
       print*, ' highlow interface metrics, outer iter', it, ': max|dT|=', dtw, &
               ' K max|dp|=', dpw, ' Pa max|du|=', duw, ' m/s'
       print*, '   LS top-row T[', minval(fp_Tw), ',', maxval(fp_Tw), '] K  p[', &
               minval(fp_pw), ',', maxval(fp_pw), '] Pa'
       open(204, file='iface_highlow.dat', status='replace')
       write(204,'(A)') '# x_w(m)  T_w(K)  p_w(Pa)  |u|_w(m/s)'
       do k = ka1, ka2
       do i = ia1, ia2
         if(pair42) then
         gi = jb + (ie1-1) - i
         xw = 0.5d0*(Bf%x(ie,gi,k)+Bf%x(ie,gi+1,k))
         else
         xw = 0.5d0*(Bf%x(i,jb,k)+Bf%x(i+1,jb,k))
         endif
         write(204,'(4ES16.7)') xw*Lscale, fp_Tw(i-ia1+1,k-ka1+1), &
         fp_pw(i-ia1+1,k-ka1+1), fp_u(i-ia1+1,k-ka1+1)
       enddo; enddo
       close(204)
     endif
!    登记本轮耦合状态（pair 12；三个判据同时满足才算“已满足”）
     call set_couple_state(1, &
        merge(1,0, dtw < Tol_Couple_Tw .and. dpw < Tol_Couple_p .and. &
                   duw < Tol_Couple_u), 12, it, dtw, Tol_Couple_Tw)
     if(it >= 2 .and. dtw < Tol_Couple_Tw .and. dpw < Tol_Couple_p .and. &
        duw < Tol_Couple_u) then
       converged = .true.
       if(my_id == 0) print*, ' Staggered highlow coupling converged at outer iter', it
       exit outer
     endif
   enddo outer

   if(.not. converged) then
     if(my_id == 0) print*, ' Staggered highlow coupling reached Niter_Couple_Outer=', &
                            Niter_Couple_Outer, ' (not fully converged)'
   endif

   call output_flow
   call output_Ts
   call output_vtk
!  final restart file so an interrupted run can be resumed from the end state
   call write_restart
  end subroutine run_staggered_highlow
!==============================================================================
! Staggered segmented coupling dispatcher.  Detects which cross-region pair
! exists in the mesh (single pair assumed) and routes to the matching driver:
!   19 FLUID<->POROUS    run_staggered_fluid_porous   (wall+q_w thermal)
!   11 FLUID<->SOLID     run_staggered_fluid_solid    (CHT, code 11)
!   13 LOWSPEED<->SOLID  run_staggered_fluid_solid    (CHT/AC, code 13)
!   12 FLUID<->LOWSPEED  run_staggered_highlow        (matching block GS/AC)
!==============================================================================
  subroutine run_staggered_multiregion(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   integer:: mBlock, ksub, nb, mBlock2, npairs
   logical:: has11, has12, has13, has19
   logical:: hasF, hasL, hasS, hasP
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2


   hasF = .false.; hasL = .false.; hasS = .false.; hasP = .false.
   do mBlock = 1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(B%Block_type == BLOCK_FLUID) hasF = .true.
     if(B%Block_type == BLOCK_LOWSPEED) hasL = .true.
     if(B%Block_type == BLOCK_SOLID) hasS = .true.
     if(B%Block_type == BLOCK_POROUS) hasP = .true.
   enddo
   has11 = hasF .and. hasS
   has12 = hasF .and. hasL
   has13 = hasL .and. hasS
   has19 = hasF .and. hasP
   npairs = 0
   if(has11) npairs = npairs+1
   if(has12) npairs = npairs+1
   if(has13) npairs = npairs+1
   if(has19) npairs = npairs+1
   if(my_id == 0) print*, ' run_staggered_multiregion: block types F/L/S/P =', hasF, hasL, hasS, hasP
   if(my_id == 0) print*, ' run_staggered_multiregion: pairs 11/12/13/19 =', &
                          has11, has12, has13, has19
   if(npairs == 0) then
     if(my_id == 0) print*, ' run_staggered_multiregion: no cross-region block combination - abort'
     return
   endif
   if(IF_Debug == 1) call dbg_dump_interfaces(nMesh)
   if(npairs > 1) then
     if(my_id == 0) print*, ' run_staggered_multiregion: multiple cross-region combos present,', &
         ' staggered driver supports a single pair (11/12/13/19) - abort'
     return
   endif

   if(has19) then
     call run_staggered_fluid_porous(nMesh)
   else if(has11) then
     call run_staggered_fluid_solid(nMesh, 11)
   else if(has13) then
     call run_staggered_fluid_solid(nMesh, 13)
   else if(has12) then
     call run_staggered_highlow(nMesh)
   endif
  end subroutine run_staggered_multiregion
!==============================================================================
! Diagnostic dump of every block-to-block interface connection (face/face1,
! node ranges and the L1..L3 connection descriptors) as seen by the coupling
! routines.  Enabled with IF_Debug = 1 in control.ec.  face numbering:
!   1=i-, 2=j-, 3=k-, 4=i+, 5=j+, 6=k+
!==============================================================================
   subroutine dbg_dump_interfaces(nMesh)
    use Global_Var
    use const_var
    implicit none
    integer:: nMesh
    integer:: mBlock, ksub
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc2

    if(my_id .ne. 0) return
    print*, ' ---- interface connection dump (block type / bc / face / nb1 / face1) ----'
    do mBlock = 1, Mesh(nMesh)%Num_Block
      B => Mesh(nMesh)%Block(mBlock)
      if(.not. associated(B%bc_msg2)) cycle
      do ksub = 1, B%subface
        Bc2 => B%bc_msg2(ksub)
        if(Bc2%nb1 <= 0 .and. .not. is_interface_bc(Bc2%bc)) cycle
        print*, '  blk', mBlock, 'type', B%Block_type, 'sub', ksub, &
                ' bc=', Bc2%bc, ' face=', Bc2%face, ' nb1=', Bc2%nb1, &
                ' face1=', Bc2%face1, ' L=', Bc2%L1, Bc2%L2, Bc2%L3
        print*, '     own  rng:', Bc2%ib, Bc2%ie, Bc2%jb, Bc2%je, Bc2%kb, Bc2%ke
        print*, '     nb   rng:', Bc2%ib1, Bc2%ie1, Bc2%jb1, Bc2%je1, Bc2%kb1, Bc2%ke1
      enddo
    enddo
    print*, ' ---------------------------------------------------------------------------'
   end subroutine dbg_dump_interfaces


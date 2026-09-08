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
   real(PRE_EC),parameter:: GS_OMEGA = 1.7d0   ! SOR over-relaxation
   integer,parameter:: MAX_GS_ITER = 20000
   integer,parameter:: MIN_GS_ITER = 5
   real(PRE_EC),parameter:: GS_TOL = 1.d-9
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
   do ii = 1, MAX_GS_ITER
!    Re-set physical BC ghost cells each iteration
     call set_solid_ghost_BC(nMesh, mBlock)
!    Exchange Ts buffer for periodic/internal interfaces
     call update_Ts_buffer_onemesh(nMesh)

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
       B%Ts(i,j,k) = (1.d0-GS_OMEGA)*T_old + GS_OMEGA*T_new

       res = max(res, abs(B%Ts(i,j,k) - T_old))
     enddo; enddo; enddo

!    Debug: print progress for first few iterations
     if(ii <= 5 .and. my_id == 0) then
       print*, "  GS iter", ii, " res=", res, " Ts(1,1,1)=", B%Ts(1,1,1), &
               " Ts(nx-1,ny-1,1)=", B%Ts(nx-1,ny-1,1)
     endif

     if(mod(ii,100) == 0 .and. ii >= MIN_GS_ITER) then
       if(res < GS_TOL) exit
     endif
   enddo

   if(my_id == 0 .or. B%Block_no == 3) then
     print*, "Solid solver Block", B%Block_no, " GS iterations:", ii, " final res:", res
     print*, "  DBG: Ts(1,1,1)=", B%Ts(1,1,1), " Ts(nx/2,ny/2,1)=", B%Ts(nx/2,ny/2,1), &
             " Ts(nx-1,ny-1,1)=", B%Ts(nx-1,ny-1,1)
   endif
  end subroutine solid_solver_one_block

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
       mb = B_n(nb)             ! Local index on this process
       if(mb <= 0) cycle        ! Neighbor not on this process
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
     subroutine couple_lowspeed_fluid_solid_face(nMesh, Bs, Bf, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: Bs, Bf
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1
      integer:: i, k, jcell_s, jcell_f, jg_s, jg_f
      real(PRE_EC):: T_s, T_f, k_s, k_f, dx_s, dx_f, T_i
      real(PRE_EC):: yf_s, yf_f

      face_s = Bc2%face
      face1  = Bc2%face1
      k_s = Bs%solid_k
      k_f = LS_k
      if(k_s <= 0.d0 .or. k_f <= 0.d0) return   ! no conduction on either side

!     j-normal interfaces (2 = j-, 5 = j+); the low-speed case tested here
      if(.not.(face_s == 2 .or. face_s == 5)) then
        print*, 'couple_lowspeed: only j-normal interfaces implemented, face=', face_s
        return
      endif
      if(face_s == 2) then
        jcell_s = Bc2%jb            ! interior cell next to j- face
        jg_s    = Bc2%jb - 1        ! ghost cell
      else
        jcell_s = Bc2%je - 1
        jg_s    = Bc2%je
      endif
      if(face1 == 2) then
        jcell_f = Bc2%jb1
        jg_f    = Bc2%jb1 - 1
      else if(face1 == 5) then
        jcell_f = Bc2%je1 - 1
        jg_f    = Bc2%je1
      else
        print*, 'couple_lowspeed: unsupported fluid face1=', face1
        return
      endif

      do k=Bc2%kb, Bc2%ke-1
      do i=Bc2%ib, Bc2%ie-1
        T_s = Bs%Ts(i,jcell_s,k)
        T_f = Bf%U(5,i,jcell_f,k)
!       distance cell-centre to face node (y direction, conformal grids)
        if(face_s == 2) then
          yf_s = Bs%y(i,Bc2%jb,k)
          dx_s = Bs%yc(i,jcell_s,k) - yf_s
        else
          yf_s = Bs%y(i,Bc2%je,k)
          dx_s = yf_s - Bs%yc(i,jcell_s,k)
        endif
        if(face1 == 2) then
          yf_f = Bf%y(i,Bc2%jb1,k)
          dx_f = Bf%yc(i,jcell_f,k) - yf_f
        else
          yf_f = Bf%y(i,Bc2%je1,k)
          dx_f = yf_f - Bf%yc(i,jcell_f,k)
        endif
        dx_s = max(dx_s, 1.d-20)
        dx_f = max(dx_f, 1.d-20)
        T_i = (k_f*T_f/dx_f + k_s*T_s/dx_s) / (k_f/dx_f + k_s/dx_s)

!       solid ghost (isothermal at T_i)
        Bs%Ts(i,jg_s,k) = 2.d0*T_i - Bs%Ts(i,jcell_s,k)

!       fluid ghost: no-slip isothermal wall at T_i
        Bf%U(1,i,jg_f,k)   = LS_rho
        Bf%U(2,i,jg_f,k)   = -Bf%U(2,i,jcell_f,k)
        Bf%U(3,i,jg_f,k)   = -Bf%U(3,i,jcell_f,k)
        Bf%U(4,i,jg_f,k)   = -Bf%U(4,i,jcell_f,k)
        Bf%U(5,i,jg_f,k)   = 2.d0*T_i - Bf%U(5,i,jcell_f,k)
        Bf%p(i,jg_f,k)     = Bf%p(i,jcell_f,k)
      enddo; enddo
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
!   Implemented for the current grid topology: solid face_s = 2 (j-),
!   fluid face1 = 1 (i-).  Other combinations are rejected with a message.
!   Conformal interface grids with aligned index directions are required.
     subroutine couple_compressible_fluid_solid_face(nMesh, Bs, Bf, Bc2)
      use Global_Var
      implicit none
      integer:: nMesh
      Type (Block_TYPE),pointer:: Bs, Bf
      TYPE (BC_MSG_TYPE),pointer:: Bc2
      integer:: face_s, face1
      integer:: i, j, k, js, jg, if1, jg1
      real(PRE_EC):: T_s, T_f, k_s_star, k_f, dx_s, dx_f, T_i_star, T_i_K
      real(PRE_EC):: p1, d1, uu1, v1, w1
      real(PRE_EC),parameter:: R_AIR = 287.0d0, RHO_REF = 1.0d0
      real(PRE_EC):: a_ref, mu_ref, cp_ref, k_ref

      face_s = Bc2%face
      face1  = Bc2%face1
      if(face_s /= 2 .or. face1 /= 1) then
        print*, 'couple_compressible: supports solid face=2, fluid face=1 only; got', face_s, face1
        return
      endif
      k_s_star = Bs%solid_k
      if(k_s_star <= 0.d0) return
      a_ref = sqrt(gamma*R_AIR*T_inf)
      mu_ref = RHO_REF*Ma*a_ref*max(Lscale,1.d-30)/Re
      cp_ref = gamma*R_AIR/(gamma-1.d0)
      k_ref = mu_ref*cp_ref/max(PrL,1.d-30)
      k_s_star = k_s_star / k_ref

!     solid (Bs): j- face, interior row j=jb, ghost j=jb-1; k aligned
      js = Bc2%jb
      jg = Bc2%jb - 1
!     fluid (Bf): i- face, interior col i=ib1, ghost i=ib1-1
      if1 = Bc2%ib1
      jg1 = Bc2%ib1 - 1
      do k = Bc2%kb, Bc2%ke-1
      do i = Bc2%ib, Bc2%ie-1
        j = Bc2%jb1 + (i - Bc2%ib)     ! fluid tangential index (conformal)
        T_s = Bs%Ts(i,js,k)
        d1  = Bf%U(1,if1,j,k)
        uu1 = Bf%U(2,if1,j,k)/d1
        v1  = Bf%U(3,if1,j,k)/d1
        w1  = Bf%U(4,if1,j,k)/d1
        p1  = (Bf%U(5,if1,j,k) - 0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gamma-1.d0)
        T_f = gamma*Ma*Ma*p1/max(d1,1.d-20)     ! non-dimensional T*
        k_f = max(Bf%mu(if1,j,k), 1.d-30)       ! non-dimensional k = mu*
        dx_s = max(Bs%y(i,Bc2%jb,k) - Bs%yc(i,js,k), 1.d-20)
        dx_f = max(Bf%yc(if1,j,k) - Bf%y(Bc2%ib1,j,k), 1.d-20)
!       interface temperature (non-dimensional balance)
        T_i_star = (k_f*T_f/max(dx_f,1.d-30) + k_s_star*(T_s/T_inf)/max(dx_s,1.d-30)) &
                 / (k_f/max(dx_f,1.d-30) + k_s_star/max(dx_s,1.d-30))
!       solid ghost: isothermal at T_i (physical K)
        T_i_K = T_i_star*T_inf
        Bs%Ts(i,jg,k) = 2.d0*T_i_K - Bs%Ts(i,js,k)
!       fluid ghost: no-slip isothermal wall at T_i* (compressible wall)
        call wall_bound_with_Tw(NVAR1, Bf%U(:,if1,j,k), Bf%U(:,jg1,j,k), &
             Ma, gamma, T_i_star, Bf%mu(if1,j,k), Bf%dw(if1,j,k), Re)
      enddo; enddo
     end subroutine couple_compressible_fluid_solid_face

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
      integer:: face_s, face1, i, k, jc_h, jc_l, jg_l, n1, n2
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
      if(face_s /= 2 .or. face1 /= 5) then
        print*, 'couple_highlow: supports comp face=2, low face=5 only; got', face_s, face1
        return
      endif
      a_ref = sqrt(gamma*R_AIR*T_inf)
      mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
      RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
      U_ref = Ma*a_ref
!     compressible block B (j- face): interior row jc_h = jb, ghost row 0
      jc_h = Bc2%jb
!     low-speed block Bn (j+ face): interior row jc_l = je1-1? ranges are in
!     nodes of the low-speed block: top face j = je1 (node), interior cell je1-1
      jc_l = Bc2%je1 - 1
      jg_l = Bc2%je1        ! first ghost cell row of the low-speed block

!     j- interface of B: ghost row j=0 ; low-speed interior at row jc_l
      do k = Bc2%kb, Bc2%ke-1
      do i = Bc2%ib, Bc2%ie-1
!       --- compressible ghost <- low-speed interior ---
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
!       --- low-speed ghost <- compressible interior ---
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
!       The low-speed (incompressible) solver cannot sustain supersonic
!       drag from the interface; clamp the exchanged velocity so the
!       low-speed model stays in its valid (Ma << 1) regime.
        uu  = max(min(uu,  U_MAX_LS), -U_MAX_LS)
        vv  = max(min(vv,  V_MAX_LS), -V_MAX_LS)
        TT  = T1s*T_inf
        do n1 = jg_l, jg_l+LAP-1
          Bn%U(1,i,n1,k) = rho
!         low-speed block stores VELOCITY in U(2..4) (U(1)=density); the old
!         rho*uu form was only harmless while rho=1 kg/m^3
          Bn%U(2,i,n1,k) = uu
          Bn%U(3,i,n1,k) = vv
          Bn%U(4,i,n1,k) = 0.d0
          Bn%U(5,i,n1,k) = TT
          Bn%p(i,n1,k)   = Bn%p(i,jc_l,k)
        enddo
      enddo; enddo

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
      integer:: face_s, face1, i, k, n1, jc_h, jc_l, jg_l, jg
      real(PRE_EC):: G, T_w, v_w, rho_w, p_phys
      real(PRE_EC):: r1, u1s, v1s, p1n, v_w_s, v_int_s, vg_s, ug_s, T_s, rho_s, E_s
      real(PRE_EC):: RHO_REF, mu_inf_p, a_ref, U_ref
      real(PRE_EC),parameter:: R_AIR = 287.0d0
      real(PRE_EC),parameter:: MU_SI0 = 1.716d-5, T_SI0 = 273.15d0, S_SI = 110.4d0

      face_s = Bc2%face
      face1  = Bc2%face1
      if(face_s /= 2 .or. face1 /= 5) then
        if(my_id == 0) print*, 'couple_compressible_porous_blowing: unsupported face pair', face_s, face1, '(skipped)'
        return
      endif
      a_ref = sqrt(gamma*R_AIR*T_inf)
      mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
      RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
      U_ref = Ma*a_ref
      G = LS_rho*max(LS_V_in,0.d0)          ! coolant mass flux kg/(m2 s), normal +y
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
!       porous ghost layers: coolant outlet (zero-gradient) + wall pressure
        do n1 = 0, LAP-1
          jg = jg_l + n1
          Bn%U(1,i,jg,k)=Bn%U(1,i,jc_l,k); Bn%U(2,i,jg,k)=Bn%U(2,i,jc_l,k)
          Bn%U(3,i,jg,k)=Bn%U(3,i,jc_l,k); Bn%U(4,i,jg,k)=Bn%U(4,i,jc_l,k)
          Bn%U(5,i,jg,k)=Bn%U(5,i,jc_l,k); Bn%p(i,jg,k) = p_phys
        enddo
!       compressible ghost layers: blowing wall (temperature T_w)
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
       mb = B_n(nb)             ! Local index on this process
       if(mb <= 0) cycle        ! Neighbor not on this process
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
  subroutine couple_lowspeed_porous_interfaces(nMesh)
   use Global_Var
   use const_var
   implicit none
   integer:: nMesh
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: mBlock, ksub, nb, mb
   integer:: face_s, face1, n, i, j, k
   integer:: i1, j1, k1          ! normal index of neighbour's first interior cell
   logical:: ok

   do mBlock=1, Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%bc_msg2)) cycle
     do ksub=1, B%subface
       Bc2 => B%bc_msg2(ksub)
       if(.not. is_interface_bc(Bc2%bc)) cycle
       nb = Bc2%nb1
       if(nb <= 0) cycle
       mb = B_n(nb)
       if(mb <= 0) cycle
       Bn => Mesh(nMesh)%Block(mb)

!      this routine couples exactly one LOWSPEED block with one POROUS block
       if(.not. ((B%Block_type == BLOCK_LOWSPEED .and. Bn%Block_type == BLOCK_POROUS) .or. &
                 (B%Block_type == BLOCK_POROUS  .and. Bn%Block_type == BLOCK_LOWSPEED))) cycle

       face_s = Bc2%face     ! face of B at the interface
       face1  = Bc2%face1    ! face of Bn at the interface
!      conformal flat interface: same normal axis, opposite low/high faces
       ok = .false.
       if(face_s == 1 .and. face1 == 4) ok = .true.
       if(face_s == 4 .and. face1 == 1) ok = .true.
       if(face_s == 2 .and. face1 == 5) ok = .true.
       if(face_s == 5 .and. face1 == 2) ok = .true.
       if(face_s == 3 .and. face1 == 6) ok = .true.
       if(face_s == 6 .and. face1 == 3) ok = .true.
       if(.not. ok) then
         if(my_id == 0) print*, 'couple_lowspeed_porous: unsupported face pair', &
                                face_s, face1, ' (block', B%Block_no, nb, ') skipped'
         cycle
       endif

!      first interior cell of Bn adjacent to its interface face
       i1 = Bc2%ib1; j1 = Bc2%jb1; k1 = Bc2%kb1
       select case(face1)
       case(4); i1 = Bc2%ie1 - 1      ! i+ face: last interior cell row ie1-1
       case(2); j1 = Bc2%jb1          ! j- face: first interior cell row jb1
       case(5); j1 = Bc2%je1 - 1      ! j+ face: last interior cell row je1-1
       case(6); k1 = Bc2%ke1 - 1      ! k+ face: last interior cell row ke1-1
       end select

!      fill B's ghost layers (n = 1..LAP) from Bn's first interior cell
       select case(face_s)
       case(1)      ! B i- face: ghosts i = ib - n
         do n=1, LAP
           i = Bc2%ib - n
           do k=Bc2%kb, Bc2%ke-1
           do j=Bc2%jb, Bc2%je-1
             B%U(1,i,j,k)=Bn%U(1,i1,j,k); B%U(2,i,j,k)=Bn%U(2,i1,j,k)
             B%U(3,i,j,k)=Bn%U(3,i1,j,k); B%U(4,i,j,k)=Bn%U(4,i1,j,k)
             B%U(5,i,j,k)=Bn%U(5,i1,j,k); B%p(i,j,k) = Bn%p(i1,j,k)
           enddo; enddo
         enddo
       case(4)      ! B i+ face: ghosts i = ie + n - 1
         do n=1, LAP
           i = Bc2%ie + n - 1
           do k=Bc2%kb, Bc2%ke-1
           do j=Bc2%jb, Bc2%je-1
             B%U(1,i,j,k)=Bn%U(1,i1,j,k); B%U(2,i,j,k)=Bn%U(2,i1,j,k)
             B%U(3,i,j,k)=Bn%U(3,i1,j,k); B%U(4,i,j,k)=Bn%U(4,i1,j,k)
             B%U(5,i,j,k)=Bn%U(5,i1,j,k); B%p(i,j,k) = Bn%p(i1,j,k)
           enddo; enddo
         enddo
       case(2)      ! B j- face: ghosts j = jb - n
         do n=1, LAP
           j = Bc2%jb - n
           do k=Bc2%kb, Bc2%ke-1
           do i=Bc2%ib, Bc2%ie-1
             B%U(1,i,j,k)=Bn%U(1,i,j1,k); B%U(2,i,j,k)=Bn%U(2,i,j1,k)
             B%U(3,i,j,k)=Bn%U(3,i,j1,k); B%U(4,i,j,k)=Bn%U(4,i,j1,k)
             B%U(5,i,j,k)=Bn%U(5,i,j1,k); B%p(i,j,k) = Bn%p(i,j1,k)
           enddo; enddo
         enddo
       case(5)      ! B j+ face: ghosts j = je + n - 1
         do n=1, LAP
           j = Bc2%je + n - 1
           do k=Bc2%kb, Bc2%ke-1
           do i=Bc2%ib, Bc2%ie-1
             B%U(1,i,j,k)=Bn%U(1,i,j1,k); B%U(2,i,j,k)=Bn%U(2,i,j1,k)
             B%U(3,i,j,k)=Bn%U(3,i,j1,k); B%U(4,i,j,k)=Bn%U(4,i,j1,k)
             B%U(5,i,j,k)=Bn%U(5,i,j1,k); B%p(i,j,k) = Bn%p(i,j1,k)
           enddo; enddo
         enddo
       case(3)      ! B k- face: ghosts k = kb - n
         do n=1, LAP
           k = Bc2%kb - n
           do j=Bc2%jb, Bc2%je-1
           do i=Bc2%ib, Bc2%ie-1
             B%U(1,i,j,k)=Bn%U(1,i,j,k1); B%U(2,i,j,k)=Bn%U(2,i,j,k1)
             B%U(3,i,j,k)=Bn%U(3,i,j,k1); B%U(4,i,j,k)=Bn%U(4,i,j,k1)
             B%U(5,i,j,k)=Bn%U(5,i,j,k1); B%p(i,j,k) = Bn%p(i,j,k1)
           enddo; enddo
         enddo
       case(6)      ! B k+ face: ghosts k = ke + n - 1
         do n=1, LAP
           k = Bc2%ke + n - 1
           do j=Bc2%jb, Bc2%je-1
           do i=Bc2%ib, Bc2%ie-1
             B%U(1,i,j,k)=Bn%U(1,i,j,k1); B%U(2,i,j,k)=Bn%U(2,i,j,k1)
             B%U(3,i,j,k)=Bn%U(3,i,j,k1); B%U(4,i,j,k)=Bn%U(4,i,j,k1)
             B%U(5,i,j,k)=Bn%U(5,i,j,k1); B%p(i,j,k) = Bn%p(i,j,k1)
           enddo; enddo
         enddo
       end select
     enddo
   enddo
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
   integer:: mf, mp, mBlock, ksub, nb, mb, it, step, pc, NVAR1
   integer:: nstep, nhalve
   integer:: face_s, face1, ib,ie,jb,je,kb,ke, ib1,ie1,jb1,je1,kb1,ke1
   Type (Block_TYPE),pointer:: Bf, Bp, B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc2
   integer:: i, k, jc_h, jc_l, jg, n1
   real(PRE_EC):: Sfac, Sfac1, twmax, qwmax
   real(PRE_EC):: a_ref, U_ref, RHO_REF, mu_inf_p, mu_ref, cp_ref
   real(PRE_EC),parameter:: R_AIR=287.d0
   real(PRE_EC),parameter:: MU_SI0=1.716d-5, T_SI0=273.15d0, S_SI=110.4d0
   real(PRE_EC):: r1, u1s, v1s, p1n, T1_nd, T1_K, mu_SI, k_gas, dxp
   real(PRE_EC):: mu1c, dwc
   logical:: found, converged
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
   Bf => Mesh(nMesh)%Block(mf)
   Bp => Mesh(nMesh)%Block(mp)
   NVAR1 = Mesh(nMesh)%NVAR

!  ---- find the FLUID code-19 subface (coupling region) ---------------------
   found = .false.
   do ksub=1, Bf%subface
     Bc2 => Bf%bc_msg2(ksub)
     if(.not. is_interface_bc(Bc2%bc)) cycle
     nb = Bc2%nb1
     if(nb <= 0) cycle
     mb = 0
     do mBlock=1, Mesh(nMesh)%Num_Block
       if(Mesh(nMesh)%Block(mBlock)%Block_no == nb) then
         Bn => Mesh(nMesh)%Block(mBlock); mb = mBlock; exit
       endif
     enddo
     if(mb == 0) cycle
     if(Bn%Block_type /= BLOCK_POROUS) cycle
     face_s = Bc2%face; face1 = Bc2%face1
     if(face_s /= 2 .or. face1 /= 5) then
       if(my_id == 0) print*, 'run_staggered: supports gas j- (face=2) vs porous j+ (face1=5) only; got', face_s, face1
       return
     endif
     ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
     ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1
     found = .true.
     exit
   enddo
   if(.not. found) then
     print*, 'run_staggered_fluid_porous: no FLUID-POROUS (code 19) interface found'
     return
   endif


!  ---- per-face work arrays (gas face cells i=ib..ie-1, k=kb..ke-1) ----------
   if(.not. allocated(fp_Tw)) then
     allocate(fp_Tw(1:Bf%nx-1, 1:Bf%nz-1))
     allocate(fp_Tw_old(1:Bf%nx-1, 1:Bf%nz-1))
     allocate(fp_qw(1:Bf%nx-1, 1:Bf%nz-1))
     allocate(fp_pw(1:Bf%nx-1, 1:Bf%nz-1))
   endif
   fp_Tw = Twall_Couple_Init
   fp_qw = 0.d0
   fp_pw = 0.d0

!  ---- physical references (same convention as couple_highlow) --------------
   a_ref   = sqrt(gamma*R_AIR*T_inf)
   mu_inf_p= MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
   RHO_REF = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
   U_ref   = Ma*a_ref
   mu_ref  = RHO_REF*U_ref*max(Lscale,1.d-30)/Re
   cp_ref  = gamma*R_AIR/(gamma-1.d0)

   Mesh(nMesh)%tt = 0.d0
   Mesh(nMesh)%Kstep = 0
   if(my_id == 0) then
     print*, ' run_staggered_fluid_porous: Kstep_Couple_Comp=', Kstep_Couple_Comp, &
             ' Niter_Couple_Outer=', Niter_Couple_Outer, &
             ' Twall_Couple_Init=', Twall_Couple_Init, ' K'
     print*, '   gas block', mf, ' (nx-1 x nz-1 =', Bf%nx-1, 'x', Bf%nz-1, &
             '), porous block', mp, '; interface j- cells i=', ib, '..', ie-1
   endif

   converged = .false.
   outer: do it=1, Niter_Couple_Outer
     fp_Tw_old = fp_Tw
!    adaptive gas-chunk length: keep the full Kstep_Couple_Comp for the first
!    Niter_Couple_Warm outer iterations (warm start), then halve it each outer
!    iteration down to Kstep_Couple_Min for the final refinement stage.
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
       call fill_gas_wall_ghost(Bf, ib, ie, jb, kb, ke, NVAR1)
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
     enddo
     if(my_id == 0) print*, ' gas chunk done, Kstep=', Mesh(nMesh)%Kstep, ' tt=', Mesh(nMesh)%tt

!     refresh the wall ghost once more against the final gas interior so the
!     extracted q_w is consistent with this chunk's isothermal wall state
      call fill_gas_wall_ghost(Bf, ib, ie, jb, kb, ke, NVAR1)


!   ==================== (2) EXTRACT q_w, p_w =========================
     jc_h = jb
     qwmax = 0.d0
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
!      wall heat flux into the wall (>0): q = k_gas*(T_int - T_ghost)/dx
       fp_qw(i,k) = k_gas*(T1_nd - T_nd_from_U(Utmp, NVAR1, Ma, gamma))/dxp * T_inf
       fp_pw(i,k) = p1n*RHO_REF*U_ref*U_ref
       qwmax = max(qwmax, abs(fp_qw(i,k)))
     enddo; enddo
     if(my_id == 0) print*, ' gas heat-flux max|q_w|=', qwmax, ' W/m2'

!   ==================== (3) POROUS CHUNK =========================
     jc_l = je1 - 1     ! porous first interior row below the hot face
     do pc=1, Porous_Chunk_Iter
!     porous hot-face boundary: outlet pressure fp_pw + heat-flux Ts BC
       do k=kb1, ke1-1
       do i=ib1, ie1-1
         do n1=0, LAP-1
           jg = je1 + n1
           Bp%U(1,i,jg,k) = Bp%U(1,i,jc_l,k)
           Bp%U(2,i,jg,k) = Bp%U(2,i,jc_l,k)
           Bp%U(3,i,jg,k) = Bp%U(3,i,jc_l,k)
           Bp%U(4,i,jg,k) = Bp%U(4,i,jc_l,k)
           Bp%U(5,i,jg,k) = Bp%U(5,i,jc_l,k)
           Bp%p(i,jg,k)   = fp_pw(i,k)
         enddo
         call set_porous_Ts_flux(Bp, i, jc_l, je1, k, fp_qw(i,k))
       enddo; enddo
       do mBlock=1, Mesh(nMesh)%Num_Block
         B => Mesh(nMesh)%Block(mBlock)
         if(B%Block_type == BLOCK_POROUS) call porous_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
       enddo
     enddo

!   ==================== (4) RETURN T_w from porous hot face ========
     do k=kb1, ke1-1
     do i=ib1, ie1-1
       fp_Tw(i,k) = 0.5d0*(Bp%Ts(i,jc_l,k) + Bp%Ts(i,jc_l+1,k))
     enddo; enddo
     twmax = 0.d0
     do k=kb1, ke1-1
     do i=ib1, ie1-1
       twmax = max(twmax, abs(fp_Tw(i,k)-fp_Tw_old(i,k)))
     enddo; enddo
     if(my_id == 0) then
       print*, ' porous chunk done, outer iter', it, ' max|dT_w|=', twmax, ' K', &
               '  T_w range [', minval(fp_Tw(ib:ie-1,kb:ke-1)), ',', &
               maxval(fp_Tw(ib:ie-1,kb:ke-1)), ']'
       open(203, file='iface_couple.dat', status='replace')
       write(203,'(A)') '# x_w(m)  T_w(K)  q_w(W/m2)  p_w(Pa)'
       do k=kb, ke-1
       do i=ib, ie-1
         write(203,'(4ES16.7)') 0.5d0*(Bf%x(i,jb,k)+Bf%x(i+1,jb,k)), &
                                fp_Tw(i,k), fp_qw(i,k), fp_pw(i,k)
       enddo; enddo
       close(203)
     endif
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

  contains

!   Impose the isothermal wall state on the gas block code-19 face cells
!   (j- face, rows j=jb.., ghost row jb-1..).  T_w = fp_Tw(i,k) [K].
    subroutine fill_gas_wall_ghost(Bf, ib, ie, jb, kb, ke, nv)
      implicit none
      Type (Block_TYPE),pointer:: Bf
      integer:: ib, ie, jb, kb, ke, nv
      integer:: i, k, n1
      real(PRE_EC):: mu1c, dwc
      do k=kb, ke-1
      do i=ib, ie-1
        if(If_viscous == 1) then
          mu1c = Bf%mu(i,jb,k); dwc = Bf%dw(i,jb,k)
        else
          mu1c = 1.d0/Re; dwc = 0.d0
        endif
        call wall_bound_with_Tw(nv, Bf%U(:,i,jb,k), Bf%U(:,i,jb-1,k), &
             Ma, gamma, fp_Tw(i,k)/T_inf, mu1c, dwc, Re)
        do n1=2, LAP
          Bf%U(:,i,jb-n1,k) = Bf%U(:,i,jb-1,k)   ! deeper ghosts (copy of layer 1)
        enddo
      enddo; enddo
    end subroutine fill_gas_wall_ghost

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

!   Solid-frame Ts ghost on the porous hot face from incoming heat flux
!   qw>0 (W/m2, heat entering the wall): Ts_g = Ts_i + qw*dx/ks_eff
    subroutine set_porous_Ts_flux(B, i, jint, jgh, k, qw)
      implicit none
      Type (Block_TYPE),pointer:: B
      integer:: i, jint, jgh, k
      real(PRE_EC):: qw, dx_g, ks_eff
      integer:: n1
      ks_eff = max((1.d0 - B%porous_eps)*B%solid_k, 1.d-30)
      do n1 = 0, LAP-1
        dx_g = sqrt( (B%xc(i,jint,k)-B%xc(i,jgh+n1,k))**2 &
                   + (B%yc(i,jint,k)-B%yc(i,jgh+n1,k))**2 &
                   + (B%zc(i,jint,k)-B%zc(i,jgh+n1,k))**2 ) * Lscale
        B%Ts(i,jgh+n1,k) = B%Ts(i,jint,k) + qw*dx_g/ks_eff
      enddo
    end subroutine set_porous_Ts_flux
  end subroutine run_staggered_fluid_porous

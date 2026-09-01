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
     call lowspeed_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   case(BLOCK_POROUS)
     call porous_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
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
       if(B%bc_msg2(ksub)%bc < 0) cycle
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
!    For face 1 (i-): node ib=1 ¡ú first cell center at i=ib, ghost cells at i=ib-1..ib-LAP
!    For face 4 (i+): node ie=nx ¡ú last cell center at i=ie-1, ghost cells at i=ie..ie+LAP-1
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
              dx_b = (B%x(ib+1,j,k) - B%x(ib,j,k)) * Lscale
              B%Ts(i1,j,k) = B%Ts(i2,j,k) + Qw_val * dx_b * n / B%solid_k
            else
              B%Ts(i1,j,k) = B%Ts(i2,j,k)
            endif
!            Debug: print ghost cell for first layer
            if(n==1 .and. my_id==0) print*, '  DBG_GHOST: face=', face_s, ' i1=', i1, ' i2=', i2, ' Ts(i1)=', B%Ts(i1,1,1), ' Ts(i2)=', B%Ts(i2,1,1), ' Tw_val=', Tw_val, ' Qw_val=', Qw_val
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
              dx_b = (B%x(ie,j,k) - B%x(ie-1,j,k)) * Lscale
              B%Ts(i1,j,k) = B%Ts(i2,j,k) + Qw_val * dx_b * n / B%solid_k
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
              dy_b = (B%y(i,jb+1,k) - B%y(i,jb,k)) * Lscale
              B%Ts(i,j1,k) = B%Ts(i,j2,k) + Qw_val * dy_b * n / B%solid_k
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
              dy_b = (B%y(i,je,k) - B%y(i,je-1,k)) * Lscale
              B%Ts(i,j1,k) = B%Ts(i,j2,k) + Qw_val * dy_b * n / B%solid_k
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
              dz_b = (B%z(i,j,kb+1) - B%z(i,j,kb)) * Lscale
              B%Ts(i,j,k1) = B%Ts(i,j,k2) + Qw_val * dz_b * n / B%solid_k
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
              dz_b = (B%z(i,j,ke) - B%z(i,j,ke-1)) * Lscale
              B%Ts(i,j,k1) = B%Ts(i,j,k2) + Qw_val * dz_b * n / B%solid_k
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
   real(PRE_EC):: T_old, dx, dy, dz, Fo, alpha, res, T_new
   real(PRE_EC):: rhoCp, hx, hy, hz
   logical:: is_unsteady
   real(PRE_EC),parameter:: GS_OMEGA = 1.7d0   ! SOR over-relaxation
   integer,parameter:: MAX_GS_ITER = 200000
   integer,parameter:: MIN_GS_ITER = 5
   real(PRE_EC),parameter:: GS_TOL = 1.d-12
   real(PRE_EC), parameter:: LARGE_FO = 1.d10

   B=>Mesh(nMesh)%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   rhoCp = B%solid_rho * B%solid_Cp
   alpha = B%solid_k / max(rhoCp, 1.d-20)

!  Determine if unsteady: RK3 or Dual-time stepping
   is_unsteady = (Time_Method == Time_RK3 .or. Time_Method == Time_Dual_LU_SGS)

!  1. Set physical BC ghost cells (wall, symmetry)
   call set_solid_ghost_BC(nMesh, mBlock)

!  2. Gauss-Seidel iteration
!  Use 2nd-order central difference stencil: update cells 2..nx-2 (excluding
!  boundary-adjacent cells that may lack valid neighbors on both sides).
   do ii = 1, MAX_GS_ITER
!    Re-set physical BC ghost cells each iteration
     call set_solid_ghost_BC(nMesh, mBlock)

!    Update all cell centers (1..nx-1, 1..ny-1, 1..nz-1).
!    Neighbors at boundaries (i=0, j=0, k=0 or i=nx, j=ny, k=nz) are
!    ghost cells properly set by set_solid_ghost_BC.
     res = 0.d0
     do k = 1, nz-1
     do j = 1, ny-1
     do i = 1, nx-1
       T_old = B%Ts(i,j,k)

!      Cell-centered spacings
       hx = B%x(i+1,j,k) - B%x(i,j,k)
       hy = B%y(i,j+1,k) - B%y(i,j,k)
       hz = B%z(i,j,k+1) - B%z(i,j,k)

!      Fourier number for implicit Euler
       if(is_unsteady .and. dt_global > 0.d0) then
         Fo = alpha * dt_global / (hx*hx + hy*hy + hz*hz + 1.d-20)
       else
!        Steady: use large Fo -> essentially Laplacian(T)=0
         Fo = LARGE_FO
       endif

!      Gauss-Seidel: average of neighbors
       T_new = (B%Ts(i+1,j,k) + B%Ts(i-1,j,k) + &
                B%Ts(i,j+1,k) + B%Ts(i,j-1,k) + &
                B%Ts(i,j,k+1) + B%Ts(i,j,k-1)) / 6.d0

!      Implicit Euler correction for unsteady
       if(is_unsteady) then
!        T_new = (T_old + 6*Fo*T_new) / (1 + 6*Fo)
         T_new = (T_old + 6.d0*Fo*T_new) / (1.d0 + 6.d0*Fo)
       endif

!      SOR: Ts_new = (1-omega)*Ts_old + omega*T_new
       B%Ts(i,j,k) = (1.d0-GS_OMEGA)*T_old + GS_OMEGA*T_new

!      Track per-iteration change for convergence check
       res = max(res, abs(B%Ts(i,j,k) - T_old))
     enddo; enddo; enddo

!    Check convergence (every 100 iterations)
     if(mod(ii,100) == 0 .and. ii >= MIN_GS_ITER) then
       if(res < GS_TOL) exit
     endif
   enddo

   if(my_id == 0 .or. B%Block_no == 3) then
     print*, "Solid solver Block", B%Block_no, " GS iterations:", ii, " final res:", res
     print*, "  DBG: Ts(1,1,1)=", B%Ts(1,1,1), " Ts(25,1,1)=", B%Ts(25,1,1), " Ts(50,1,1)=", B%Ts(50,1,1)
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
       if(Bc2%bc >= 0) cycle
       nb = Bc2%nb1             ! Neighbor block number (global)
       if(nb <= 0) cycle
       mb = B_n(nb)             ! Local index on this process
       if(mb <= 0) cycle        ! Neighbor not on this process
       Bn => Mesh(nMesh)%Block(mb)

!      Determine if this is a fluid-solid interface
       if(.not. ((B%Block_type == BLOCK_FLUID .and. Bn%Block_type == BLOCK_SOLID) .or. &
                 (B%Block_type == BLOCK_SOLID .and. Bn%Block_type == BLOCK_FLUID))) cycle

       face_s = Bc2%face
       ib=Bc2%ib; ie=Bc2%ie; jb=Bc2%jb; je=Bc2%je; kb=Bc2%kb; ke=Bc2%ke
       face1 = Bc2%face1
       ib1=Bc2%ib1; ie1=Bc2%ie1; jb1=Bc2%jb1; je1=Bc2%je1; kb1=Bc2%kb1; ke1=Bc2%ke1

!      Determine which side is fluid and which is solid
       if(B%Block_type == BLOCK_FLUID .and. Bn%Block_type == BLOCK_SOLID) then
!        B is fluid, Bn is solid
         call compute_interface_T(nMesh, B, Bn, ksub, Bc2, T_i)
!        Set fluid ghost cell (isothermal wall at T_i) for current block (B)
         call set_fluid_ghost_wall_face(nMesh, mBlock, &
             face_s, ib, ie, jb, je, kb, ke, T_i)
!        Set solid ghost cell for neighbor (Bn) using Bc2 face1/ib1/ie1/...
         call set_solid_ghost_from_interface_face(nMesh, mb, &
             face1, ib1, ie1, jb1, je1, kb1, ke1, T_i, Bn)
       elseif(B%Block_type == BLOCK_SOLID .and. Bn%Block_type == BLOCK_FLUID) then
!        B is solid, Bn is fluid
         call compute_interface_T(nMesh, B, Bn, ksub, Bc2, T_i)
!        Set solid ghost cell for current block (B) using Bc face/ib/ie/...
         call set_solid_ghost_from_interface_face(nMesh, mBlock, &
             face_s, ib, ie, jb, je, kb, ke, T_i, B)
!        Set fluid ghost cell (isothermal wall at T_i) for neighbor (Bn) using Bc2 face1/ib1/ie1/...
         call set_fluid_ghost_wall_face(nMesh, mb, &
             face1, ib1, ie1, jb1, je1, kb1, ke1, T_i)
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
       if(Bc2%bc >= 0) cycle
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
! Low-speed flow solver stub
!----------------------------------------------------------------------
  subroutine lowspeed_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_Var
   implicit none
   integer:: nMesh, mBlock
   real(PRE_EC):: Sfac, Sfac1
!  TODO: implement low-speed flow solver
  end subroutine lowspeed_solver_one_block

!----------------------------------------------------------------------
! Porous media solver stub
!----------------------------------------------------------------------
  subroutine porous_solver_one_block(nMesh, mBlock, Sfac, Sfac1)
   use Global_Var
   implicit none
   integer:: nMesh, mBlock
   real(PRE_EC):: Sfac, Sfac1
!  TODO: implement porous media solver
  end subroutine porous_solver_one_block

!----------------------------------------------------------------------
! Output solid temperature field (Ts) for all solid blocks
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
     if(B%Block_type /= BLOCK_SOLID) cycle
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
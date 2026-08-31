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
! Solid heat conduction solver stub
!----------------------------------------------------------------------
  subroutine solid_solver_one_block(nMesh, mBlock)
   use Global_Var
   implicit none
   integer:: nMesh, mBlock
!  TODO: implement solid heat conduction solver
  end subroutine solid_solver_one_block

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
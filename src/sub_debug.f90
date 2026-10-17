   !===============================================================================
!  区域残差统计打印（2026-09-26）：固体 / 低速 / 多孔求解器共用。
!
!  每个块（区域）不再只报一个 max 残差，而是报三个量；口径与可压块既有的
!  Mesh%Res_max / Mesh%Res_rms 完全一致（sub_time_advance.f90::
!  comput_max_Res_onemesh），使四类块用同一种"残差语言"：
!     res_max : 逐胞残差幅值在区域内的最大值
!     res_rms : sqrt( sum(res^2) / N_cell )，即区域平均（均方根）水平
!     (imax,jmax,kmax) : res_max 所在胞，同时给出该胞中心坐标 [m]
!
!  只由 0 号进程调用/打印；参数全为标量（无字符量），调用方不需要显式接口。
!  kk = 迭代/扫掠号：>0 打印过程监视行；<0 打印结束汇总行（|kk| = 总次数）。
!  what 选择残差名称与附加字段（a1/a2，保持既有日志字段 res_q/res_m/res_p 不丢）：
!     1 固体 GS |dTs|        2 低速 AC 流场          3 多孔 AC 流场
!     4 低速 SIMPLE div(F)   5 多孔 SIMPLE div(F)
!     6 低速 AC 能量 |dT|    7 多孔 LTNE 流体 |dTf|   8 多孔 LTNE 骨架 |dTs|
!===============================================================================
  module region_stat_mod
   use precision_EC
   use const_var
   implicit none
!  过程监视行间隔：硬编码的显示常量（不是 namelist 参数，不参与 6 步流程）
   integer,parameter:: RES_PRINT_STEP = 100
  contains

  subroutine region_res_stat(nMesh, mBlock, what, kk, res_max, imax, jmax, kmax, res_rms, a1, a2)
   use Global_var
   implicit none
   integer:: nMesh, mBlock, what, kk, imax, jmax, kmax
   real(PRE_EC):: res_max, res_rms, a1, a2
   Type (Block_TYPE),pointer:: B
   character(len=32):: ctype, cname
   character(len=20):: cknd
   character(len=160):: aux
   character(len=320):: buf

   B => Mesh(nMesh)%Block(mBlock)
   select case(B%Block_type)
   case(BLOCK_FLUID)   ; ctype = 'BLOCK_FLUID'
   case(BLOCK_SOLID)   ; ctype = 'BLOCK_SOLID'
   case(BLOCK_LOWSPEED); ctype = 'BLOCK_LOWSPEED'
   case(BLOCK_POROUS)  ; ctype = 'BLOCK_POROUS'
   case default        ; ctype = 'BLOCK_UNKNOWN'
   end select
   select case(abs(what))
   case(1); cname = 'solid GS |dTs|'
   case(2); cname = 'low-speed AC flow'
   case(3); cname = 'porous AC flow'
   case(4); cname = 'low-speed SIMPLE div(F)'
   case(5); cname = 'porous SIMPLE div(F)'
   case(6); cname = 'low-speed AC energy |dT|'
   case(7); cname = 'porous LTNE fluid |dTf|'
   case(8); cname = 'porous LTNE frame |dTs|'
   case default; cname = 'residual'
   end select

!  附加字段：与既有日志字段兼容，便于老脚本 grep
   aux = ' '
   if(what == 2 .or. what == 3) then
     write(aux,'("  res_q=",E12.5," res_m=",E12.5)') a1, a2
   else if(what == 4 .or. what == 5) then
     write(aux,'("  res_p=",E12.5)') a1
   else if(abs(what) >= 6) then
     write(aux,'("  dT(window)=",E12.5)') a1
   endif

   if(kk > 0) then
     write(cknd,'(" step=",I6)') kk
   else if(kk < 0) then
     write(cknd,'(" final=",I6)') -kk
   else
     cknd = ' final'
   endif

   write(buf,'(1x,A," [",A,"] block",I4," (",A,")",A," : max=",E15.8, &
        "  rms=",E15.8,"  at (i,j,k)=(",3I5,")  xyz=[",3(1x,E13.6),"] m",A)') &
        ' Region res', trim(cname), B%Block_no, trim(ctype), trim(cknd), res_max, &
        res_rms, imax, jmax, kmax, B%xc(imax,jmax,kmax)*Lscale, &
        B%yc(imax,jmax,kmax)*Lscale, B%zc(imax,jmax,kmax)*Lscale, trim(aux)
   write(*,'(A)') trim(buf)
  end subroutine region_res_stat

  end module region_stat_mod


subroutine Output_mesh_debug
   use Global_var
   implicit none
   integer:: m,nx,ny,nz,i,j,k
   Type (Block_TYPE),pointer:: B
   character(len=50):: filename
   write(filename,"('Meshdeb-'I5.5'.dat')") my_id
   open(50,file=filename)
    
!  输出网格文件 tecplot格式；   
   write(50,*) "variables=x,y,z"
   do m=1,Mesh(1)%Num_Block
     B=>Mesh(1)%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
     write(50,*) "zone i=",nx+2, " j= ", ny+2, " k= ", nz+2
	 do k=0,nz+1
	 do j=0,ny+1
	 do i=0,nx+1
	 write(50,"(3f20.10)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k)
	 enddo
	 enddo
	 enddo
	enddo
	close(50)
   end
 

 
   subroutine debug1
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,m,i,j,k
   real*8:: u1,u2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    u1=0.d0
	u2=0.d0
	do m=1,MP%Num_Block
    B=>MP%Block(m)
	nx=B%nx; ny=B%ny; nz=B%nz
	do k=1-LAP,nz+LAP-1
	do j=1-LAP,ny+LAP-1
	do i=1-LAP,nx+LAP-1
	 u1=u1+B%U(1,i,j,k)**2
	 u2=u2+B%U(5,i,j,k)**2
	enddo
	enddo
	enddo
	enddo
	print*, "************ u1,u2=",u1,u2
	end



   subroutine debug2(nx,ny,nz)
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,i,j,k
   real*8:: f1,f2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    f1=0.d0
	f2=0.d0
	do k=1,nz
	do j=1,ny
	do i=1,nx
	 f1=f1+Flux(1,i,j,k)**2
     f2=f2+Flux(5,i,j,k)**2
	enddo
	enddo
	enddo

	print*, "************ f1,f2=",f1,f2
	end

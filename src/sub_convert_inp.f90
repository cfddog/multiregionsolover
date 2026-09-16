!----Boundary message (bc3d.inp, Gridgen general format)----------------------------------------------------------
! 读取边界格式文件(bc3d.inp)，并转化为OpenCFD-EC的内建.inc格式 (见《OpenCFD-EC理论手册》 2.4.3节)
! OpenCFD-EC内建的存储格式比bc3d.inp多了一些冗余信息 （有些类似BXCFD的.in格式），例如多了子面号f_no,
! 面类型face 以及连接的子面号f_no1,连接的面类型face1
! 以及连接次序L1, L2, L3  (例如L1=1表示该维与连接块的第1维正连接， L1=-1表示与连接块的第1为反向连接).
! 这些冗余信息为块-块之间的通信（尤其是MPI并行通信）提供了便利，有利于简化代码
 
 module  Type_def1
   use const_var
    TYPE BC_MSG_TYPE              ! 边界链接信息
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                     ! 子面号，连接顺序描述符
   END TYPE BC_MSG_TYPE


    TYPE Block_TYPE1                          ! 数据结构：仅包含Bc_msg 
	 integer::  nx,ny,nz                      ! 网格数nx,ny,nz
	 integer::  subface                       ! 子面数
 	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
    END TYPE Block_TYPE1   
 End module  Type_def1

!------------------------------------------------------
  subroutine convert_inp_inc 
   use Type_def1
   implicit none
  
   integer:: NB,m,ksub,nx,ny,nz,k,j,k1,ksub1
   integer:: kb(3),ke(3),kb1(3),ke1(3),s(3),p(3),Lp(3)
   TYPE(Block_TYPE1),Pointer,dimension(:):: Block
   Type (Block_TYPE1),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
     
   Interface   
     subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
	   integer,dimension(3):: kb,ke,kb1,ke1
	 end subroutine
    end Interface
   

   print*, "Convert bc3d.inp to bc3d.inc ..."
   open(88,file="bc3d.inp")
   read(88,*)
   read(88,*) NB
   allocate(Block(NB))
   do m=1,NB
    B => Block(m)
    read(88,*) B%nx,B%ny,B%nz
	read(88,*)
    read(88,*) B%subface   !number of the subface in the Block m
     
	allocate(B%bc_msg(B%subface))   ! 边界描述

    do ksub=1, B%subface
      Bc => B%bc_msg(ksub)
      Bc%f_no=ksub                        ! 子面号
	  read(88,*)  kb(1),ke(1),kb(2),ke(2),kb(3),ke(3),Bc%bc
!     A face left without a boundary code gets NO flux in every solver, which
!     silently produces a wrong (often diverging) solution -- warn loudly.
      if(Bc%bc .eq. 0) then
        print*, ' WARNING: bc3d.inp face with boundary code 0 (unset) !!!'
        print*, '   block=', m, ' subface=', ksub, &
                ' ib,ie,jb,je,kb,ke=', kb(1),ke(1),kb(2),ke(2),kb(3),ke(3)
        print*, '   -> no flux is applied there; give it a real code', &
                ' (2=wall 3=symmetry 5=inlet 6=outlet -1=inner link)'
      endif
	 
	  if(Bc%bc .lt. 0) then
 !  --------有连接的情况 (内边界)--------------------------------------------------------	    
	   read(88,*) kb1(1),ke1(1),kb1(2),ke1(2),kb1(3),ke1(3),Bc%nb1
      else
 !---------无连接情况 (物理边界)----------------------------------------
            kb1(:)=0; ke1(:)=0; Bc%nb1=0
	  endif
	 call Convert_bc(Bc,kb,ke,kb1,ke1)
  enddo
  enddo
   
   close(88)

!  搜索连接块的块号 f_no1  (便于MPI并行通信是使用)
   do m=1,NB
     B => Block(m)
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .lt. 0) then
         Bc%f_no1=0
		 B1=>Block(Bc%nb1)         ! 指向连接块
         
		 do ksub1=1,B1%subface
		 Bc1=>B1%bc_msg(ksub1)
		 if(Bc%ib1==Bc1%ib .and. Bc%ie1==Bc1%ie .and. Bc%jb1==Bc1%jb .and. Bc%je1==Bc1%je   &
		    .and. Bc%kb1==Bc1%kb .and. Bc%ke1==Bc1%ke  .and. Bc1%nb1==m) then
		    Bc%f_no1=ksub1
		  exit
		 endif 	 
         enddo

         if(Bc%f_no1 ==0) then
		  print*, " Error in find linked block number !!!"
		  print*, "Block, subface=",m, ksub
		  stop
		 endif
	   endif
      enddo
	enddo

   open(99,file="bc3d.inc")
    write(99,*) " Inp-liked file of OpenCFD-EC"
    write(99,*) NB
	do m=1,NB
     B => Block(m)
     write(99,*) B%nx, B%ny, B%nz
	 write(99,*) "Block ", m
	 write(99,*) B%subface
	  do ksub=1,B%subface
        Bc => B%bc_msg(ksub)
	   write(99,"(9I6)") Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
	   write(99,"(12I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
      enddo
	enddo
	close(99)

   print*, "Convert bc3d.inp to bc3d.inc OK"


  end  

!-----------------------------------------------------------------------------------




!   将Gridgen格式 转换为 OpenCFD-EC 的边界连接格式
!     计算面号、 连接次序 (L1,L2,L3)等
      subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
	   integer,dimension(3):: kb,ke,kb1,ke1,s,p,LP
       integer:: k,j,k1

	   if(Bc%bc .ge. 0) then
	     Bc%ib1=0; Bc%ie1=0; Bc%jb1=0; Bc%je1=0; Bc%kb1=0; Bc%ke1=0; Bc%nb1=0
         Bc%L1=0; Bc%L2=0; Bc%L3=0; Bc%face1=0; Bc%f_no1=0
       endif
     

!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
       do k=1,3
  	     if(kb(k) .eq. ke(k) ) then 
	       s(k)=0                           ! 连接维
	     else if (kb(k) .gt. 0) then 
	       s(k)=1                           ! 正
	     else
	       s(k)=-1                          ! 负
	     endif
       enddo

!    边界子面的大小     
	 Bc%ib=min(abs(kb(1)),abs(ke(1))) ;  Bc%ie=max(abs(kb(1)),abs(ke(1)))
     Bc%jb=min(abs(kb(2)),abs(ke(2))) ;  Bc%je=max(abs(kb(2)),abs(ke(2)))
     Bc%kb=min(abs(kb(3)),abs(ke(3))) ;  Bc%ke=max(abs(kb(3)),abs(ke(3))) 


!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
      if(s(1) .eq. 0) then
	     if (Bc%ib .eq. 1) then
	      Bc%face=1               ! i-
	     else
	      Bc%face=4               ! i+
	    endif
      else if(s(2) .eq. 0) then
	    if(Bc%jb .eq. 1) then
	      Bc%face=2                 ! j-
	    else
	      Bc%face=5                 ! j+
	    endif
      else
 	   if(Bc%kb  .eq. 1) then
	    Bc%face=3                 ! k-
	   else
	    Bc%face=6                 ! k+
	   endif
     endif 

!---------------------------------------------------------------------------
!------内边界的情况，建立连接描述
  if( Bc%bc .lt. 0) then            ! 内边界
!      计算连接顺序描述符L1,L2,L3
!      计算各维之间的连接关系      
     do k=1,3  
	   if(kb1(k) .eq. ke1(k) ) then 
	       p(k)=0                      ! 
       else if (kb1(k) .gt. 0) then
	       p(k)=1                      ! .inp 文件的 正数
       else
	       p(k)=-1
       endif
     enddo
 	   

!    对应连接子面的大小     
	 Bc%ib1=min(abs(kb1(1)),abs(ke1(1))) ;  Bc%ie1=max(abs(kb1(1)),abs(ke1(1)))
     Bc%jb1=min(abs(kb1(2)),abs(ke1(2))) ;  Bc%je1=max(abs(kb1(2)),abs(ke1(2)))
     Bc%kb1=min(abs(kb1(3)),abs(ke1(3))) ;  Bc%ke1=max(abs(kb1(3)),abs(ke1(3))) 
    
 	  
!   判断该面连接面的类型 (i-, i+, j-,j+, k-,k+)     
      if(p(1) .eq. 0) then
	     if (Bc%ib1 .eq. 1) then
	      Bc%face1=1               ! i-
	     else
	      Bc%face1=4               ! i+
	    endif
      else if(p(2) .eq. 0) then
	    if(Bc%jb1 .eq. 1) then
	      Bc%face1=2                 ! j-
	    else
	      Bc%face1=5                 ! j+
	    endif
      else
 	   if(Bc%kb1 .eq. 1) then
	    Bc%face1=3
	   else
	    Bc%face1=6
	   endif
     endif  

!  计算“连接对” 描述符  bc%L1, bc%L2, bc%L3 
 	   do k=1,3
	     do j=1,3
	       if(s(k) .eq. p(j)) Lp(k)=j          ! .inp文件的连接格式： 正对正、 负对负、 0对0； 
	     enddo
	   enddo    
	   
!    计算连接次序 （正为顺序；负为拟序）	  
	  do k=1,3
	   if(s(k) .ne. 0) then
	     k1=Lp(k)
	     if( (ke(k)-kb(k))*(ke1(k1)-kb1(k1)) .lt. 0) Lp(k)=-Lp(k)      ! 逆序连接
	   else
         k1=Lp(k)
!		 if( (mod(Bc%face,2)-mod(Bc%face1,2))==0) Lp(k)=-Lp(k)         ! 逆序连接 （单面） ! Bug 2012-7-13
! 正-正连接 Lp为负 (例， i+ 面连接到 j+面， 则为逆序连接)
		 if( (Bc%face-1)/3 .eq. (Bc%face1-1)/3 ) Lp(k)=-Lp(k)        ! (Bc%face=1,2,3为 +面，4,5,6为-面)  ! 逆序连接 （单面）

	   endif
	 enddo 
     
	 Bc%L1=Lp(1); Bc%L2=Lp(2); Bc%L3=Lp(3)   ! 连接次序描述符 （详见《理论手册》）
   endif
  end

!----------------------------------------------------------------------
! Same as convert_inp_inc but reads from bc3d_interface.inp and writes to bc3d_interface.inc
! Keeps ALL block-to-block connections (does NOT convert any to physical boundaries)
!----------------------------------------------------------------------
  subroutine convert_inp_inc_interface
   use Type_def1
   implicit none
   integer:: NB,m,ksub,nx,ny,nz,k,j,k1,ksub1
   integer:: kb(3),ke(3),kb1(3),ke1(3),s(3),p(3),Lp(3)
   TYPE(Block_TYPE1),Pointer,dimension(:):: Block
   Type (Block_TYPE1),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
   Interface
     subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
       integer,dimension(3):: kb,ke,kb1,ke1
     end subroutine
    end Interface
   print*, "Convert bc3d_interface.inp to bc3d_interface.inc ..."
   open(88,file="bc3d_interface.inp")
   read(88,*)
   read(88,*) NB
   allocate(Block(NB))
   do m=1,NB
    B => Block(m)
    read(88,*) B%nx,B%ny,B%nz
    read(88,*)
    read(88,*) B%subface
    allocate(B%bc_msg(B%subface))
    do ksub=1, B%subface
      Bc => B%bc_msg(ksub)
      Bc%f_no=ksub
      read(88,*)  kb(1),ke(1),kb(2),ke(2),kb(3),ke(3),Bc%bc
      if(is_interface_bc(Bc%bc)) then
        read(88,*) kb1(1),ke1(1),kb1(2),ke1(2),kb1(3),ke1(3),Bc%nb1
      else
        kb1(:)=0; ke1(:)=0; Bc%nb1=0
      endif
      call Convert_bc(Bc,kb,ke,kb1,ke1)
    enddo
   enddo
   close(88)
!  Search for linked block f_no1
   do m=1,NB
     B => Block(m)
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(is_interface_bc(Bc%bc) .and. Bc%nb1 .gt. 0) then
         Bc%f_no1=0
         B1=>Block(Bc%nb1)
         do ksub1=1,B1%subface
           Bc1=>B1%bc_msg(ksub1)
           if(Bc%ib1==Bc1%ib .and. Bc%ie1==Bc1%ie .and. Bc%jb1==Bc1%jb .and. Bc%je1==Bc1%je &
              .and. Bc%kb1==Bc1%kb .and. Bc%ke1==Bc1%ke .and. Bc1%nb1==m) then
             Bc%f_no1=ksub1
             exit
           endif
         enddo
         if(Bc%f_no1 ==0) then
           print*, " Error in find linked block number !!!"
           print*, "Block, subface=",m, ksub
           stop
         endif
       endif
     enddo
   enddo
   open(99,file="bc3d_interface.inc")
    write(99,*) " Inp-liked file of OpenCFD-EC (interface)"
    write(99,*) NB
    do m=1,NB
     B => Block(m)
     write(99,*) B%nx, B%ny, B%nz
     write(99,*) "Block ", m
     write(99,*) B%subface
      do ksub=1,B%subface
        Bc => B%bc_msg(ksub)
       write(99,"(9I6)") Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
       write(99,"(12I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
      enddo
    enddo
   close(99)
   print*, "Convert bc3d_interface.inp to bc3d_interface.inc OK"
  end
!----------------------------------------------------------------------
! Auto-prepare bc3d.inp from bc3d_interface.inp + material.in
!
! bc3d_interface.inp is treated as the single source of truth for the block
! connections (and normally for the physical faces as well).  The explicit
! cross-class interface codes (11..19, see const_var) are derived from the
! block types given in material.in, so the interface codes of bc3d.inp no
! longer have to be kept in sync by hand.
!
! Iflag_bc_check (control.ec):
!   0 : disabled (legacy behaviour, bc3d.inp must be supplied)
!   1 : if bc3d.inp is missing -> generate it from bc3d_interface.inp;
!       if bc3d.inp is present -> validate it and, when inconsistent or
!       unreadable, back it up (bc3d.inp.bak) and regenerate it
!   2 : strict - validate bc3d.inp and stop when inconsistent
!
! On regeneration the physical (non-interface) face codes of the existing
! bc3d.inp are preserved (read back from bc3d.inp.bak); only the interface
! codes are re-derived.  Explicit codes 11..19 carry no inline continuation
! line; negative codes keep the inline neighbour line, matching convert_inp_inc.
!----------------------------------------------------------------------
  subroutine prepare_bc3d_input
   use Global_var
   use const_var
   implicit none
   logical:: ex_if, ex_bc, ex_mat, have_types, ok
   integer:: nt, m, ios
   integer,allocatable:: btl(:)

   if(Iflag_bc_check .le. 0) return

   inquire(file='bc3d_interface.inp',exist=ex_if)
   if(.not. ex_if) return
   inquire(file='bc3d.inp',exist=ex_bc)
   inquire(file='material.in',exist=ex_mat)

   have_types = .false.
   nt = 0
   if(ex_mat) then
     open(97,file='material.in')
     read(97,*,iostat=ios) nt
     if(ios==0 .and. nt>0) then
       allocate(btl(nt))
       read(97,*,iostat=ios) (btl(m),m=1,nt)
       if(ios==0) have_types=.true.
     endif
     close(97)
   endif
   if(.not. have_types) then
     print*, 'prepare_bc3d_input: material.in unavailable -> interface faces ', &
             'written as BC_Inner (-1); coupling still dispatches by block type'
     if(.not. allocated(btl)) then
       nt = 1
       allocate(btl(1))
       btl(1) = BLOCK_FLUID
     endif
   endif

   if(.not. ex_bc) then
     call gen_bc3d(btl, nt, have_types, '')
     print*, 'prepare_bc3d_input: bc3d.inp generated from bc3d_interface.inp'
   else
     call check_bc3d(btl, nt, have_types, ok)
     if(ok) then
       print*, 'prepare_bc3d_input: bc3d.inp consistent with bc3d_interface.inp, keep it'
     else if(Iflag_bc_check .eq. 2) then
       print*, 'prepare_bc3d_input: ERROR - bc3d.inp is inconsistent with ', &
               'bc3d_interface.inp/material.in (Iflag_bc_check=2)'
       stop
     else
       print*, 'prepare_bc3d_input: bc3d.inp inconsistent -> save as bc3d.inp.bak, ', &
               'regenerate interface codes (physical codes kept from the backup)'
       call backup_bc3d
       call gen_bc3d(btl, nt, have_types, 'bc3d.inp.bak')
     endif
   endif

   if(allocated(btl)) deallocate(btl)

  contains

!   Unordered pair of block types -> canonical explicit interface code.
!   Returns 0 when no explicit code is defined (same-class pair).
   integer function canonical_interface_code(ta,tb)
     integer,intent(in):: ta, tb
     canonical_interface_code = 0
     if((ta==BLOCK_FLUID .and. tb==BLOCK_SOLID) .or. &
        (ta==BLOCK_SOLID .and. tb==BLOCK_FLUID)) then
       canonical_interface_code = BC_Interface_FluidSolid
     else if((ta==BLOCK_FLUID .and. tb==BLOCK_LOWSPEED) .or. &
             (ta==BLOCK_LOWSPEED .and. tb==BLOCK_FLUID)) then
       canonical_interface_code = BC_Interface_FluidLow
     else if((ta==BLOCK_LOWSPEED .and. tb==BLOCK_SOLID) .or. &
             (ta==BLOCK_SOLID .and. tb==BLOCK_LOWSPEED)) then
       canonical_interface_code = BC_Interface_LowSolid
     else if(ta==BLOCK_SOLID .and. tb==BLOCK_SOLID) then
       canonical_interface_code = BC_Interface_SolidSolid
     else if((ta==BLOCK_LOWSPEED .and. tb==BLOCK_POROUS) .or. &
             (ta==BLOCK_POROUS .and. tb==BLOCK_LOWSPEED)) then
       canonical_interface_code = BC_Interface_LowPorous
     else if((ta==BLOCK_SOLID .and. tb==BLOCK_POROUS) .or. &
             (ta==BLOCK_POROUS .and. tb==BLOCK_SOLID)) then
       canonical_interface_code = BC_Interface_SolidPorous
     else if(ta==BLOCK_POROUS .and. tb==BLOCK_POROUS) then
       canonical_interface_code = BC_Interface_PorousPorous
     else if((ta==BLOCK_FLUID .and. tb==BLOCK_POROUS) .or. &
             (ta==BLOCK_POROUS .and. tb==BLOCK_FLUID)) then
       canonical_interface_code = BC_Interface_FluidPorous
     endif
   end function canonical_interface_code

   subroutine backup_bc3d
     integer:: ios
     character(len=256):: line
     open(96,file='bc3d.inp',status='old')
     open(95,file='bc3d.inp.bak',status='replace')
     do
       read(96,'(a)',iostat=ios) line
       if(ios/=0) exit
       write(95,'(a)') trim(line)
     enddo
     close(96); close(95)
   end subroutine backup_bc3d

!   Read bc3d_interface.inp (continuation line after every interface code)
!   and write bc3d.inp with the canonical interface codes.  Same-class faces
!   become BC_Inner (-1).  Physical face codes come from the interface file,
!   unless ovname is non-empty: then they are copied from that file where the
!   (normalised) face ranges match (used to preserve user physical BCs).
!   That overlay file is read entirely up front; if it cannot be parsed it is
!   simply ignored so a malformed bc3d.inp never aborts the run.
   subroutine gen_bc3d(tlist, nblk_t, has_types, ovname)
     integer,intent(in):: nblk_t
     integer,intent(in):: tlist(*)
     logical,intent(in):: has_types
     character(len=*),intent(in):: ovname
     integer:: NB, NBu, m, k, ksub, nsub, nx,ny,nz, oc, ios, tot, t
     integer:: ib,ie,jb,je,kb,ke,bc
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1
     integer:: ru(6), kk(6), cc, nxr,nyr,nzr, nsr
     character(len=256):: line
     integer,allocatable:: oblk(:), ofac(:,:)
     logical:: ov, found

     ov = (len_trim(ovname) > 0)
     NBu = 0
     tot = 0
     if(ov) then
       open(86,file=trim(ovname),status='old',iostat=ios)
       if(ios/=0) then
         ov = .false.
       else
         read(86,'(a)',iostat=ios) line
         if(ios==0) read(86,*,iostat=ios) NBu
         if(ios/=0) then
           ov = .false.; NBu = 0
         endif
       endif
       if(ov) then
         ios = 0
         do m=1,NBu
           read(86,*,iostat=ios) nxr,nyr,nzr
           if(ios==0) read(86,'(a)',iostat=ios) line
           if(ios==0) read(86,*,iostat=ios) nsr
           if(ios/=0) exit
           do k=1,nsr
             read(86,*,iostat=ios) kk(1),kk(2),kk(3),kk(4),kk(5),kk(6),cc
             if(ios/=0) exit
             if(cc .lt. 0) read(86,'(a)',iostat=ios) line
             if(ios/=0) exit
             tot = tot + 1
           enddo
           if(ios/=0) exit
         enddo
         if(ios/=0) ov = .false.
       endif
       if(ov) then
         allocate(oblk(tot), ofac(7,tot))
         rewind(86)
         read(86,'(a)') line
         read(86,*) NBu
         t = 0
         do m=1,NBu
           read(86,*) nxr,nyr,nzr
           read(86,'(a)') line
           read(86,*) nsr
           do k=1,nsr
             read(86,*) kk(1),kk(2),kk(3),kk(4),kk(5),kk(6),cc
             t = t + 1
             oblk(t) = m
             ofac(1,t)=min(abs(kk(1)),abs(kk(2))); ofac(2,t)=max(abs(kk(1)),abs(kk(2)))
             ofac(3,t)=min(abs(kk(3)),abs(kk(4))); ofac(4,t)=max(abs(kk(3)),abs(kk(4)))
             ofac(5,t)=min(abs(kk(5)),abs(kk(6))); ofac(6,t)=max(abs(kk(5)),abs(kk(6)))
             ofac(7,t)=cc
             if(cc .lt. 0) read(86,'(a)') line
           enddo
         enddo
       endif
       close(86)
     endif

     open(88,file='bc3d_interface.inp')
     open(99,file='bc3d.inp',status='replace')
     read(88,'(a)') line
     write(99,'(a)') 'OpenCFD-EC bc3d.inp generated from bc3d_interface.inp'
     read(88,*) NB
     write(99,*) NB
     do m=1,NB
       read(88,*) nx,ny,nz
       write(99,'(3I10)') nx,ny,nz
       read(88,'(a)') line
       write(99,'(a)') trim(line)
       read(88,*) nsub
       write(99,*) nsub
       do ksub=1,nsub
         read(88,*) ib,ie,jb,je,kb,ke,bc
         if(is_interface_bc(bc)) then
           read(88,*) ib1,ie1,jb1,je1,kb1,ke1,nb1
           if(has_types .and. m<=nblk_t .and. nb1>=1 .and. nb1<=nblk_t) then
             if(tlist(m) .ne. tlist(nb1)) then
               oc = canonical_interface_code(tlist(m), tlist(nb1))
               if(oc .le. 0) oc = BC_Inner
             else
               oc = BC_Inner
             endif
           else
             oc = BC_Inner
           endif
         else
           oc = bc
           if(ov .and. tot>0) then
             ru(1)=min(abs(ib),abs(ie)); ru(2)=max(abs(ib),abs(ie))
             ru(3)=min(abs(jb),abs(je)); ru(4)=max(abs(jb),abs(je))
             ru(5)=min(abs(kb),abs(ke)); ru(6)=max(abs(kb),abs(ke))
             do t=1,tot
               if(oblk(t)==m .and. ofac(1,t)==ru(1) .and. ofac(2,t)==ru(2) .and. &
                  ofac(3,t)==ru(3) .and. ofac(4,t)==ru(4) .and. &
                  ofac(5,t)==ru(5) .and. ofac(6,t)==ru(6)) then
                 oc = ofac(7,t); exit
               endif
             enddo
           endif
         endif
         write(99,'(6I10,I10)') ib,ie,jb,je,kb,ke,oc
         if(oc .lt. 0) write(99,'(7I10)') ib1,ie1,jb1,je1,kb1,ke1,nb1
       enddo
     enddo
     close(88); close(99)
     if(allocated(oblk)) deallocate(oblk,ofac)
   end subroutine gen_bc3d

!   Validate a user-supplied bc3d.inp against bc3d_interface.inp + block types.
!   Face ranges are normalised (abs + min/max) exactly like Convert_bc, so the
!   '-1/-2' and '1/2' encodings compare equal.  Continuation lines: interface
!   file -> after every is_interface_bc code; user file -> after a negative
!   code only (see convert_inp_inc).  Any read error (e.g. an inline line after
!   a positive code, which convert_inp_inc cannot consume) marks the file as
!   inconsistent so it gets regenerated instead of aborting.
   subroutine check_bc3d(tlist, nblk_t, has_types, ok)
     integer,intent(in):: nblk_t
     integer,intent(in):: tlist(*)
     logical,intent(in):: has_types
     logical,intent(out):: ok
     integer:: NBi, NBu, m, k, j, nsubi, nsubu, nx,ny,nz, ios
     integer:: ib,ie,jb,je,kb,ke,bc
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1
     integer:: cu, expected
     character(len=256):: line
     integer,allocatable:: aki(:,:), aci(:), anbi(:)
     integer,allocatable:: aku(:,:), acu(:)
     logical:: found

     ok = .true.
     open(88,file='bc3d_interface.inp',iostat=ios)
     if(ios/=0) then
       ok = .false.; return
     endif
     open(87,file='bc3d.inp',iostat=ios)
     if(ios/=0) then
       ok = .false.; close(88); return
     endif
     read(88,'(a)',iostat=ios) line
     if(ios==0) read(88,*,iostat=ios) NBi
     if(ios==0) read(87,'(a)',iostat=ios) line
     if(ios==0) read(87,*,iostat=ios) NBu
     if(ios/=0) then
       print*, 'check_bc3d: cannot read header of bc3d.inp/bc3d_interface.inp'
       ok=.false.; close(87); close(88); return
     endif
     if(NBi .ne. NBu) then
       print*, 'check_bc3d: block count differs: bc3d_interface=',NBi,' bc3d.inp=',NBu
       ok = .false.
     endif
     do m=1,min(NBi,NBu)
       read(88,*,iostat=ios) nx,ny,nz
       if(ios==0) read(88,'(a)',iostat=ios) line
       if(ios==0) read(88,*,iostat=ios) nsubi
       if(ios/=0) then
         print*, 'check_bc3d: cannot parse bc3d_interface.inp block',m
         ok=.false.; close(87); close(88); return
       endif
       allocate(aki(6,nsubi), aci(nsubi), anbi(nsubi))
       do k=1,nsubi
         read(88,*,iostat=ios) ib,ie,jb,je,kb,ke,bc
         if(ios/=0) exit
         aki(1,k)=min(abs(ib),abs(ie)); aki(2,k)=max(abs(ib),abs(ie))
         aki(3,k)=min(abs(jb),abs(je)); aki(4,k)=max(abs(jb),abs(je))
         aki(5,k)=min(abs(kb),abs(ke)); aki(6,k)=max(abs(kb),abs(ke))
         aci(k)=bc; anbi(k)=0
         if(is_interface_bc(aci(k))) then
           read(88,*,iostat=ios) ib1,ie1,jb1,je1,kb1,ke1,nb1
           if(ios/=0) exit
           anbi(k) = nb1
         endif
       enddo
       if(ios/=0) then
         print*, 'check_bc3d: cannot parse bc3d_interface.inp block',m
         ok=.false.; close(87); close(88); return
       endif
       read(87,*,iostat=ios) nx,ny,nz
       if(ios==0) read(87,'(a)',iostat=ios) line
       if(ios==0) read(87,*,iostat=ios) nsubu
       if(ios/=0) then
         print*, 'check_bc3d: cannot parse bc3d.inp block',m
         ok=.false.; close(87); close(88); return
       endif
       allocate(aku(6,nsubu), acu(nsubu))
       do k=1,nsubu
         read(87,*,iostat=ios) ib,ie,jb,je,kb,ke,bc
         if(ios/=0) exit
         aku(1,k)=min(abs(ib),abs(ie)); aku(2,k)=max(abs(ib),abs(ie))
         aku(3,k)=min(abs(jb),abs(je)); aku(4,k)=max(abs(jb),abs(je))
         aku(5,k)=min(abs(kb),abs(ke)); aku(6,k)=max(abs(kb),abs(ke))
         acu(k)=bc
         if(bc .lt. 0) read(87,'(a)',iostat=ios) line
         if(ios/=0) exit
       enddo
       if(ios/=0) then
         print*, 'check_bc3d: cannot parse bc3d.inp block',m, &
                 '(inline continuation after a positive code?)'
         ok=.false.; close(87); close(88); return
       endif
       do k=1,nsubi
         if(.not. is_interface_bc(aci(k))) cycle
         found = .false.
         j = 0
         do while(.not. found .and. j < nsubu)
           j = j + 1
           if(aku(1,j)==aki(1,k) .and. aku(2,j)==aki(2,k) .and. &
              aku(3,j)==aki(3,k) .and. aku(4,j)==aki(4,k) .and. &
              aku(5,j)==aki(5,k) .and. aku(6,j)==aki(6,k)) found = .true.
         enddo
         if(.not. found) then
           print*, 'check_bc3d: block',m,' interface face missing in bc3d.inp:', &
                   aki(1,k),aki(2,k),aki(3,k),aki(4,k),aki(5,k),aki(6,k)
           ok = .false.
           cycle
         endif
         cu = acu(j)
         if(has_types .and. m<=nblk_t .and. anbi(k)>=1 .and. anbi(k)<=nblk_t) then
           if(tlist(m) .ne. tlist(anbi(k))) then
             expected = canonical_interface_code(tlist(m), tlist(anbi(k)))
             if(cu .ne. expected .and. cu .ne. BC_Wall) then
               print*, 'check_bc3d: block',m,' face',aki(1,k),aki(2,k),aki(3,k), &
                       aki(4,k),aki(5,k),aki(6,k),' code',cu,' expected',expected
               ok = .false.
             endif
           else
             if(cu .ge. 0 .and. cu .ne. BC_Wall) then
               if((tlist(m)==BLOCK_FLUID .or. tlist(m)==BLOCK_LOWSPEED) .or. &
                  .not. is_interface_bc(cu)) then
                 print*, 'check_bc3d: block',m,' same-class face',aki(1,k), &
                         aki(2,k),' code',cu,' should be an internal (negative) code'
                 ok = .false.
               endif
             endif
           endif
         else
           if(cu .ge. 0 .and. cu .ne. BC_Wall .and. .not. is_interface_bc(cu)) then
             print*, 'check_bc3d: block',m,' interface face',aki(1,k),aki(2,k), &
                     'coded as physical boundary',cu
             ok = .false.
           endif
         endif
       enddo
       deallocate(aki,aci,anbi,aku,acu)
     enddo
     close(88); close(87)
   end subroutine check_bc3d

  end subroutine prepare_bc3d_input

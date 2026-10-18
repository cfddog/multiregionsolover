
  module Flow_Var
   use precision_EC
   real(PRE_EC), save,pointer,dimension(:,:,:)::  d,uu,v,w,T,p,cc ! ????x-????y-????z-?????????????
   real(PRE_EC), save,pointer,dimension(:,:,:,:):: Flux                 ! i- ,j-??k-????????
   real(PRE_EC), save,pointer,dimension(:,:,:):: Lvi,Lvj,Lvk,Lci,Lcj,Lck   ! ??????????Jocabian????? 

  end module Flow_Var

!------------------------------------------------------------------------------------------
! ?????? ?????? ??????
!-----------------------------------------------------------------------------------------
  program main
   use Global_Var
   implicit none
   integer:: ierr,mBlock,kt_in,m
   real(PRE_EC):: Sfac,Sfac1,max_res
   logical:: use_staggered
   logical:: tight_continuation
   call Init_mpi
  
   if(my_id .eq. 0) then
    print*,  "--------- MultiPhysMixSolver ver 1.0 (MPI, based on OpenCFD-EC 1.16a) --------"
    print*,  "        Multi-physics multi-block solver: compressible / low-speed AC /       "
    print*,  "        solid / porous regions with cross-region coupling                    "
    print*,  "----------------------------------------------------------------------------- " 
   endif

   call read_parameter                     ! ????????????????????
!$ call omp_set_num_threads(NUM_THREADS)   ! ????OpenMP??????????? ????????????? ??????openmp?????????????!

!$ if(my_id ==0) then         ! ??????????��???? ??openmp???????????????
!$OMP Parallel
!$  print*, "omp run ..."
!$OMP END parallel
!$ endif 

   allocate( Mesh(Num_Mesh) )                                 ! ????????? ?????? ?????????????�B??
   
   if(my_id .eq. 0)  call prepare_bc3d_input
   if(my_id .eq. 0)  call check_mesh_multigrid               ! ??????????????????????????,???څ?????????????
   
   call Init                               ! ?????????????????��???????
   call set_control_para                   ! ?څ????????????????????????????????????????????????????????
   call check_mesh_quality                 ! ???????????,?????????????????????????
   call Init_flow                          ! ????????? ???????
!  Halo-exchange self test (IF_Debug = 3): fills U(1:NVAR) and the low-speed /
!  porous pressure B%p with a unique per-(block,i,j,k) sentinel, runs the
!  same-class buffer exchange and verifies EVERY ghost cell against the mapped
!  neighbour value (same-process and cross-process subfaces alike).  The state
!  is restored before the run continues, so this is non-destructive.
!  Set IF_Debug=3 in control.ec and read the "[bufchk]" lines at start-up.
   if(IF_Debug == 3) then
     if(my_id == 0) print*, " Halo-exchange self test (IF_Debug=3) ......"
     call check_buffer_exchange(1)
   endif
!  Cross-process safety net: abort automatically when this mesh mixes block
!  types whose interface has no MPI version (12/14/17/18/19) and the two
!  sides are on different processes.  np = 1 returns immediately.
   call check_remote_coupling_support(1)
!  Parallel AC march: whole-job pressure-BC flag + AC_MPI_Sync safety checks
!  (start-up only; np=1 and jobs without AC blocks are unaffected)
   call init_ac_parallel_flags(1)

   if(my_id .eq. 0) print*, " Start ......"

   if(IF_Debug == 1) call dbg_dump_interfaces(1)

!------------------------------------------------------------------------
! Staggered segmented coupling: dispatches by detected cross-region pair
! (11 FLUID<->SOLID CHT, 12 FLUID<->LOWSPEED matching block-GS, 13
! LOWSPEED<->SOLID CHT, 19 FLUID<->POROUS thermal).  Otherwise the
! standard per-step simultaneous multi-block time loop runs.
   if(Iflag_Couple_Scheme == 1) then
     use_staggered = .true.
     tight_continuation = .false.
!    Restart-time coupling-state decision (2026-09-25, control.ec $couple_ec
!    Iflag_Couple_Restart):
!      0 = auto (default): if the previous run's staggered coupling had already
!          satisfied Tol_Couple_Tw, do NOT replay the staggered schedule --
!          continue with the per-step (tight) coupling in the normal time loop
!          (t_end then applies again).
!      1 = force the tight (per-step) coupling after a restart
!      2 = force the staggered driver after a restart (pre-2026-09-25 behaviour)
!     -1 = ignore the coupling state stored in the restart file completely
     if(Iflag_Couple_Restart >= 0 .and. restart_found == 1) then
       if(Iflag_Couple_Restart == 1) then
         use_staggered = .false.
       else if(Iflag_Couple_Restart == 0 .and. Couple_State_Mode >= 0 .and. &
               Couple_State_Conv == 1) then
         use_staggered = .false.
       endif
     endif
     if(.not. use_staggered) then
!      ---- continue with the tight (per-step) coupling -----------------------
!      The interface T_w/q_w are recomputed from the restored U/Ts on every
!      step, so changing modes introduces no jump.  The solid GS budget is
!      capped: with the $solid_ec defaults (20000 sweeps / tol 1e-9) a single
!      per-step solid solve was measured at ~15-20 s for cases/fluid_solid,
!      which makes a long continuation impractical.  Only the sweep limit and
!      the tolerance are touched (relaxation factor and physics unchanged) and
!      the effective values are printed / echoed in output_para.out.
       Iflag_Couple_Scheme = 0
       if(Solid_Max_Iter > 20) Solid_Max_Iter = 20
       if(Solid_Tol < 1.d-6)   Solid_Tol     = 1.d-6
!      report what was read from the restart file BEFORE set_couple_tight_state
!      below overwrites the globals
       if(my_id .eq. 0) then
         if(Couple_State_Mode == 0) then
           print*, ' restart: previous run already used the per-step (tight)', &
                   ' coupling -> keep marching to t_end=', t_end, &
                   ' (Solid_Max_Iter=', Solid_Max_Iter, ' Solid_Tol=', Solid_Tol, ')'
         else
           print*, ' restart: staggered coupling already satisfied (last metric=', &
                   Couple_State_Twmax, ' <= tol=', Couple_State_Tol, ', pair=', &
                   Couple_State_Pair, ' after', Couple_State_Iter, ' outer iters)'
           print*, '   -> switch to per-step (tight) coupling:', &
                   ' Iflag_Couple_Scheme=0, Solid_Max_Iter=', Solid_Max_Iter, &
                   ' Solid_Tol=', Solid_Tol, ' ; marching to t_end=', t_end
         endif
       endif
       call set_couple_tight_state()   ! keep this state in later restart files
       tight_continuation = .true.
     else
       if(my_id .eq. 0) print*, ' Staggered segmented coupling mode (Iflag_Couple_Scheme=1) ...'
       call run_staggered_multiregion(1)
       call mpi_finalize(ierr)
       stop
     endif
   endif

!------------------------------------------------------------------------
! ??????????????????????????????????? ????1??Euler??3??RK
   do while(Mesh(1)%tt .lt. t_end )
    
     call show_Wall_time()
     if(Num_Mesh .eq. 1) then                          ! ???????????1?????
!==========================================================================================
!   ??????????????: ????????????
!   ??????: ???????? block ??? ?��????+???????(??????),
!   ?????????????????????��????????????????????�D????????
       select case(Time_Method)
       case(Time_RK3)                                  ! 3??RK: ???3?????(stage)
         call comput_Sfac(Sfac,Sfac1)                  ! =0
         do mBlock=1,Mesh(1)%Num_Block
           call Set_Un_oneblock(1,mBlock)              ! Un=U
         enddo
         do KRK=1,3
           do mBlock=1,Mesh(1)%Num_Block               ! ??????? block, ?????????
             call solver_one_block(1,mBlock,Sfac,Sfac1)
           enddo
           if(Total_proc > 1 .and. Solid_Sync_Chunk > 0) &
             call solid_solve_all_blocks(1)   ! P3: same solid sweep count on all ranks
           call couple_fluid_solid_interfaces(1)       ! ????-???��???????
           call couple_solid_solid_interfaces(1)       ! ????-???��??????????
            call couple_lowspeed_porous_interfaces(1)   ! low-speed-porous (16)

           if( IFLAG_LIMIT_FLOW == 1) call limit_flow(1)
           call Boundary_condition_onemesh(1)          ! ??????? (?څGhost Cell???)
           call update_buffer_onemesh(1)               ! ?????????????
           call update_p_buffer_onemesh(1)             ! low-speed AC: p 鬼点(每外层步一次)

           call update_Ts_buffer_onemesh(1)            ! ???????????????
         enddo

       case(Time_Dual_LU_SGS)                          ! ????????: ?????
         do kt_in=1, step_inner_Limit
           call comput_Sfac(Sfac,Sfac1)                ! ??????????? (????????????, ?????????)
           do mBlock=1,Mesh(1)%Num_Block               ! ??????? block, ?????????
             call solver_one_block(1,mBlock,Sfac,Sfac1)
           enddo
           if(Total_proc > 1 .and. Solid_Sync_Chunk > 0) &
             call solid_solve_all_blocks(1)   ! P3: same solid sweep count on all ranks
           call couple_fluid_solid_interfaces(1)       ! ????-???��???????
            call couple_lowspeed_porous_interfaces(1)   ! low-speed-porous (16)

           call couple_solid_solid_interfaces(1)       ! ????-???��??????????
           if( IFLAG_LIMIT_FLOW == 1) call limit_flow(1)
           call Boundary_condition_onemesh(1)          ! ??????? (?څGhost Cell???)
           call update_buffer_onemesh(1)               ! ?????????????
           call update_p_buffer_onemesh(1)             ! low-speed AC: p 鬼点(每外层步一次)

           call update_Ts_buffer_onemesh(1)            ! ???????????????
           call comput_max_Res_onemesh(1)              ! ???????/???????��?
           max_res=Mesh(1)%Res_rms(1)
           do m=1,Mesh(1)%NVAR
             max_res=max(max_res,Mesh(1)%Res_rms(m))
           enddo
           if( max_res .le. Res_Inner_Limit) exit      ! ???��??????????????
         enddo
         if(my_id .eq. 0) then
           print*, "Inner step ... ", kt_in
           print*, "rms residual eq =", Mesh(1)%Res_rms(1:Mesh(1)%NVAR)
         endif
         do mBlock=1,Mesh(1)%Num_Block
           call Set_Un1_Un_oneblock(1,mBlock)          ! Un1=Un; Un=U
         enddo

       case default                                    ! LU_SGS + 1 step Euler
         call comput_Sfac(Sfac,Sfac1)                  ! =0
         call Set_Un(1)                                ! Un=U
         do mBlock=1,Mesh(1)%Num_Block
           call solver_one_block(1,mBlock,Sfac,Sfac1)
!          NB: the low-porous interface is coupled once per step (after the
!          block loop) so that the MPI call count is independent of how many
!          blocks this process owns; keeping this in-loop call would make the
!          1-proc and N-proc coupling cadence differ.
!          if(LS_Algorithm .eq. 3) call couple_lowspeed_porous_interfaces(1)
         enddo
         if(Total_proc > 1 .and. Solid_Sync_Chunk > 0) &
           call solid_solve_all_blocks(1)   ! P3: same solid sweep count on all ranks
         call couple_fluid_solid_interfaces(1)       ! ????-???��???????
         call couple_solid_solid_interfaces(1)       ! ????-???��??????????
         call couple_lowspeed_porous_interfaces(1)   ! low-speed-porous (16)
         if( IFLAG_LIMIT_FLOW == 1) call limit_flow(1)
         call Boundary_condition_onemesh(1)            ! ??????? (?څGhost Cell???)
         call update_buffer_onemesh(1)                 ! ?????????????
         call update_p_buffer_onemesh(1)               ! low-speed AC: p 鬼点(每外层步一次，各进程计数恒定)
         call update_Ts_buffer_onemesh(1)            ! ???????????????
       end select

       call force_vt_kw(1)                             ! ????vt,Kt,Wt???
       Mesh(1)%tt=Mesh(1)%tt+dt_global                 ! ???
       Mesh(1)%Kstep=Mesh(1)%Kstep+1                   ! ??????
	 else  if(Num_Mesh .eq. 2)  then                    ! 2?????????1?????
  	  call NS_2stge_multigrid
     else                                               ! 3?????????1?????
  	  call NS_3stge_multigrid
     endif
 
 !  ??? ,????????????. ??Kstep_Filter=0????????   
	if(Kstep_smooth .gt. 0) then
 	  if(mod(Mesh(1)%Kstep, Kstep_smooth).eq.0)   call Filtering_oneMesh(1)                      ! ???            
    endif


!  ??????????????????????��??????????????: force.log, Residual.dat??
     if(mod(Mesh(1)%Kstep, Kstep_show).eq.0) then
      call comput_force
      call output_Res(1)
     endif
! ????????????????????(flow3d.dat, PLOT3D ???)
      if(mod(Mesh(1)%Kstep, Kstep_Save).eq.0) then
	     call output_flow 
!         if(If_debug == 1 ) call output_vt                ! ?????????????????debug???   ! Bug 2017-5-11
          if(If_debug == 1 .and.  If_viscous==1 .and.  Iflag_turbulence_model .ne. 0) call output_vt                ! ?????????????????debug???
          call output_Ts                                    ! ???????????
          call output_vtk                                   ! ???VTK???(Paraview??????)
	  endif
!  restart file field_restart.dat + node-centred SI flow field flow3d_node.dat
!  (self-scheduled by Kstep_restart / Kstep_save; must run on ALL ranks)
      call restart_step_output(1)
 
 ! ??????????
    if(Kstep_average > 0) then     
      if(mod(Mesh(1)%Kstep, Kstep_average) .eq.0) then
         call Time_average           ! ??????
	  endif
      if(mod(Mesh(1)%Kstep, Kstep_Save).eq.0) then
	      call output_flow_average   ! ????????,PLOT3D???
	  endif    
    endif 
   
   enddo

!  交错→逐步强耦合的续算路径：退出前补写一次重启文件（否则只有 Kstep 命中
!  Kstep_restart/Kstep_save 时才写，t_end 落在两次写出之间就没有文件可续）
   if(tight_continuation .and. Iflag_restart >= 0) call write_restart

   call mpi_finalize(ierr)   ! 结束 MPI，避免 mpirun 报异常退出
  end
!----------------------------------------------------------------------------------------------
  subroutine Init_mpi
   use Global_var
   implicit none
      integer,parameter:: IBuffer_Size=10000000    
!      real*8,allocatable,dimension(:)::  Buffer_mpi    ! Buffer for MPI  message transfer (used by MPI_Bsend)
       real*8:: Buffer_mpi(IBuffer_Size)
      integer   ierr, status(MPI_status_size)

!------------------------------------------------
       call mpi_init(ierr)                                     ! ?????MPI
       call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)           ! ???????????
       call mpi_comm_size(MPI_COMM_WORLD,Total_proc,ierr)      
!       allocate(Buffer_mpi(IBuffer_Size))
	   call MPI_BUFFER_ATTACH(Buffer_mpi,8*IBuffer_Size,ierr)   ! ????????????????????MPI_Bsend()???
   end subroutine Init_mpi


!  ???????????????????MPI????��??	  
    subroutine show_Wall_time()
      use Global_var
	  real*8:: wtime
	  real*8,save:: wtime0,wtime1    ! ????????????????
	  integer,save:: KP=0  ! ????
      if(my_id .eq. 0) then
	    wtime=MPI_Wtime()
        if(KP .eq. 0) then   
		  wtime0=wtime   ! ???CPU???
		else
          if(mod(Mesh(1)%Kstep, Kstep_Show).eq.0) then
		  print*, "CPU wall time in this step:", wtime-wtime1 
          print*, "Averaged CPU wall time is:", (wtime-wtime0)/KP 
          endif
        endif
		 wtime1=wtime  
         KP=KP+1    ! ??????
      endif

    end subroutine show_wall_time
       

!----------------------------------------------------------------------------
!  include "sub_interfaces.f90"
!  include "sub_init.f90"
!  include "sub_Residual.f90"
!  include "sub_flux_split.f90"
!  include "sub_boundary.f90"
!  include "sub_turbulence_BL.f90"
!  include "sub_turbulence_SA.f90"
!  include "sub_turbulence_SST.f90" 
!  include "sub_time_advance.f90"
!  include "sub_updata_buffer_mpi.f90"
!  include "sub_Post.f90"
!  include "sub_time_acceleraction.f90"
!  include "sub_LU_SGS.f90"
!  include "sub_IO.f90"
!  include "sub_geometry.f90"
!  include "sub_partation_mpi.f90"
!  include "sub_convert_inp.f90"
!  include "sub_read_parameter.f90"

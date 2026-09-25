!==============================================================================
! Restart file (field_restart.dat) for OpenCFD-EC 1.16a
!
!  File layout (unformatted, double precision):
!    record 1 : iver, Total_block, Num_Mesh, NVAR, Kstep, LAP          (integers)
!    record 2 : tt, Ma, Re, gamma, T_inf, Lscale, dt_global             (reals)
!    then, for every block m = 1..Total_block (GLOBAL block order, i.e. the
!    same order as flow3d.dat / Mesh3d.x):
!      record : nx, ny, nz
!      record : U (1:NVAR, 1-LAP:nx+LAP-1, 1-LAP:ny+LAP-1, 1-LAP:nz+LAP-1)
!      record : Ts  (1-LAP:.., ..)      solid/porous frame temperature
!      record : Tsn (1-LAP:.., ..)      previous time step (= Ts if unused)
!      record : p   (1-LAP:.., ..)      low-speed / porous pressure
!
!    iver = 2 additionally appends ONE trailer record describing the
!    cross-region coupling state (see set_couple_state):
!      record : Couple_State_Mode, Conv, Pair, Iter, nf, nk
!      record : Couple_State_Twmax, Couple_State_Tol
!      record : Couple_Tw_save, Couple_qw_save, Couple_u_save, Couple_pw_save
!               (only when nf,nk > 0; interface face-cell arrays)
!    The trailer lets a restarted run decide whether the staggered coupling
!    is already satisfied (-> continue with the per-step tight coupling) and
!    lets a continued staggered run restore the interface T_w/q_w instead of
!    re-seeding them from control.ec (zero jump across the restart).
!    Files written before this change (iver = 1) are still read fine; the
!    coupling state is then simply unknown (old behaviour).
!
!  The WHOLE LAP-deep buffer (ghost-cell) box is stored, so the buffer
!  information of the previous run is recovered exactly; inter-block buffers are
!  additionally re-exchanged once after reading.
!
!  Controls (control.ec, namelist $flow_ec):
!    Iflag_restart  = 0 auto (read when present) | 1 force | -1 never
!    Kstep_restart  = write interval; <=0 -> use Kstep_save
!==============================================================================

!------------------------------------------------------------------------------
! Write the restart file.  Called by ALL ranks: rank 0 owns the file and
! receives the blocks owned by the other ranks (same pattern as output_flow).
!------------------------------------------------------------------------------
  subroutine write_restart
   use Global_var
   implicit none
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: NVAR1,m,m1,i,j,k,nx,ny,nz,mt,nxe,nye,nze,NU,NT,tag,ierr
   integer:: Send_to_ID, status(MPI_status_size)
   real(PRE_EC),allocatable,dimension(:,:,:,:):: Ubuf
   real(PRE_EC),allocatable,dimension(:,:,:):: Tsbuf,Tnbuf,Pbuf

   MP=>Mesh(1)
   NVAR1=MP%NVAR

   if(my_id .eq. 0) then
     open(99,file=RESTART_FILE,form="unformatted",status="replace")
     write(99) 2, Total_block, Num_Mesh, NVAR1, Mesh(1)%Kstep, LAP
     write(99) Mesh(1)%tt, Ma, Re, gamma, T_inf, Lscale, dt_global
   endif

   do m=1,Total_block
     nx=bNi(m); ny=bNj(m); nz=bNk(m)
     nxe=nx+2*LAP-1; nye=ny+2*LAP-1; nze=nz+2*LAP-1
     NU=NVAR1*nxe*nye*nze; NT=nxe*nye*nze
     if(B_proc(m) .eq. 0) then
!      this block is owned by rank 0 -> write it directly
       if(my_id .eq. 0) then
         mt=B_n(m); B=>MP%Block(mt)
         write(99) nx,ny,nz
         write(99) ((((B%U(m1,i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1),m1=1,NVAR1)
         write(99) (((B%Ts(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
         write(99) (((B%Tsn(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
         write(99) (((B%p(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
       endif
     else
!      the block lives on another rank
       if(my_id .eq. 0) then
         allocate(Ubuf(NVAR1,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Tsbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Tnbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Pbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         tag=B_n(m)+200
         call MPI_Recv(Ubuf,NU,OCFD_DATA_TYPE,B_proc(m),tag,MPI_COMM_WORLD,status,ierr)
         call MPI_Recv(Tsbuf,NT,OCFD_DATA_TYPE,B_proc(m),tag+1,MPI_COMM_WORLD,status,ierr)
         call MPI_Recv(Tnbuf,NT,OCFD_DATA_TYPE,B_proc(m),tag+2,MPI_COMM_WORLD,status,ierr)
         call MPI_Recv(Pbuf,NT,OCFD_DATA_TYPE,B_proc(m),tag+3,MPI_COMM_WORLD,status,ierr)
         write(99) nx,ny,nz
         write(99) ((((Ubuf(m1,i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1),m1=1,NVAR1)
         write(99) (((Tsbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
         write(99) (((Tnbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
         write(99) (((Pbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                      k=1-LAP,nz+LAP-1)
         deallocate(Ubuf,Tsbuf,Tnbuf,Pbuf)
       endif
     endif
   enddo

!  coupling-state trailer (iver = 2): mode/conv/pair/iter/nf/nk + scalars + the
!  interface T_w/q_w/u/p_w arrays.  Written only when a staggered driver (or the
!  tight-coupling switch) registered a state; rank 0 still owns unit 99 here.
   if(my_id .eq. 0 .and. Couple_State_Mode .ge. 0 .and. &
      Couple_State_nf .gt. 0 .and. Couple_State_nk .gt. 0 .and. &
      allocated(Couple_Tw_save)) then
     write(99) Couple_State_Mode, Couple_State_Conv, Couple_State_Pair, &
               Couple_State_Iter, Couple_State_nf, Couple_State_nk
     write(99) Couple_State_Twmax, Couple_State_Tol
     write(99) Couple_Tw_save, Couple_qw_save, Couple_u_save, Couple_pw_save
   else if(my_id .eq. 0 .and. Couple_State_Mode .ge. 0) then
     write(99) Couple_State_Mode, Couple_State_Conv, Couple_State_Pair, &
               Couple_State_Iter, 0, 0
     write(99) Couple_State_Twmax, Couple_State_Tol
   endif

!  send phase for the ranks that do not own the file
   if(my_id .ne. 0) then
     do m=1,MP%Num_Block
       B=>MP%Block(m)
       nx=B%nx; ny=B%ny; nz=B%nz
       nxe=nx+2*LAP-1; nye=ny+2*LAP-1; nze=nz+2*LAP-1
       NU=NVAR1*nxe*nye*nze; NT=nxe*nye*nze
       tag=m+200
       call MPI_Send(B%U,  NU,OCFD_DATA_TYPE,0,tag,  MPI_COMM_WORLD,ierr)
       call MPI_Send(B%Ts, NT,OCFD_DATA_TYPE,0,tag+1,MPI_COMM_WORLD,ierr)
       call MPI_Send(B%Tsn,NT,OCFD_DATA_TYPE,0,tag+2,MPI_COMM_WORLD,ierr)
       call MPI_Send(B%p,  NT,OCFD_DATA_TYPE,0,tag+3,MPI_COMM_WORLD,ierr)
     enddo
   else
     close(99)
   endif

   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0) print*, " write ", trim(RESTART_FILE), " OK  (Kstep=", &
       Mesh(1)%Kstep, ", tt=", Mesh(1)%tt, ")"
  end subroutine write_restart

!------------------------------------------------------------------------------
! Read the restart file (if present and consistent with the current case).
!   ok = .true.  : state (U, Ts, Tsn, p, Kstep, tt, buffer layers) restored
!   ok = .false. : no usable restart file -> the caller keeps its normal
!                  initialisation (uniform flow / Iflag_init path)
! With Iflag_restart = 1 a missing/inconsistent file is a FATAL error.
!------------------------------------------------------------------------------
  subroutine read_restart(ok)
   use Global_var
   implicit none
   logical,intent(out):: ok
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   logical:: ex
   integer:: iver,nb_file,nmesh_file,nvar_file,kst_file,lap_file
   integer:: nx,ny,nz,nx1,ny1,nz1,nxe,nye,nze,NU,NT,tag,ierr,ios
   integer:: m,m1,i,j,k,mt,status(MPI_status_size),iok,kst_read
   real(PRE_EC):: r4(7),tt_read
   real(PRE_EC),allocatable,dimension(:,:,:,:):: Ubuf
   real(PRE_EC),allocatable,dimension(:,:,:):: Tsbuf,Tnbuf,Pbuf

   MP=>Mesh(1)
   ok=.false.; iok=0; kst_read=0; tt_read=0.d0

   if(my_id .eq. 0) then
     inquire(file=RESTART_FILE,exist=ex)
     if(.not. ex) then
       print*, " read_restart: '", trim(RESTART_FILE), "' not found"
       if(Iflag_restart .eq. 1) iok=-1
     else
       open(99,file=RESTART_FILE,form="unformatted",status="old")
       read(99,iostat=ios) iver,nb_file,nmesh_file,nvar_file,kst_file,lap_file
       if(ios .ne. 0) then
         print*, " read_restart: cannot read the file header (iostat=",ios,")"
         if(Iflag_restart .eq. 1) iok=-1
       else if(iver .ne. 1 .and. iver .ne. 2) then
         print*, " read_restart: file version", iver, " is not supported"
         if(Iflag_restart .eq. 1) iok=-1
       else if(nb_file .ne. Total_block .or. nvar_file .ne. MP%NVAR .or. &
               lap_file .ne. LAP) then
         print*, " read_restart: file/case mismatch  nblock,nvar,LAP =", &
                 nb_file, nvar_file, lap_file, "  (now", Total_block, MP%NVAR, LAP,")"
         if(Iflag_restart .eq. 1) iok=-1
       else
         read(99,iostat=ios) r4(1:7)
         read(99,iostat=ios) nx1,ny1,nz1
         if(ios .ne. 0) then
           print*, " read_restart: cannot read the block-1 header"
           if(Iflag_restart .eq. 1) iok=-1
         else if(nx1 .ne. bNi(1) .or. ny1 .ne. bNj(1) .or. nz1 .ne. bNk(1)) then
           print*, " read_restart: file was written for another mesh  block1 =", &
                   nx1,ny1,nz1, "  (now", bNi(1),bNj(1),bNk(1),")"
           if(Iflag_restart .eq. 1) iok=-1
         else
           iok=1
           kst_read=kst_file; tt_read=r4(1)
           if(abs(r4(2)-Ma) .gt. 1.d-8*max(1.d0,abs(Ma)) .or. &
              abs(r4(6)-Lscale) .gt. 1.d-8*max(1.d0,abs(Lscale))) then
             print*, " read_restart: WARNING Ma/Lscale in the file =", r4(2), r4(6), &
                     " but control.ec has", Ma, Lscale
           endif
           rewind(99)
           read(99) iver,nb_file,nmesh_file,nvar_file,kst_file,lap_file
           read(99) r4(1:7)
         endif
       endif
       if(iok .ne. 1) close(99)
     endif
   endif

   call MPI_Bcast(iok,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   if(iok .lt. 0) then
     if(my_id .eq. 0) print*, " Iflag_restart=1 requires a usable restart file -> stop"
     stop 1
   endif
   if(iok .eq. 0) then
     if(my_id .eq. 0) print*, " -> start from the uniform / Iflag_init state"
     return
   endif

!---- transfer the per-block state (rank 0 owns the file) ---------------------
!  Master: walk the GLOBAL block order; read the file and push remote blocks to
!  their owners.  Workers: receive their OWN blocks in local order (the block
!  partition assigns a contiguous global range per rank, so the orders match --
!  the same convention as read_flow_data / output_flow).
   if(my_id .eq. 0) then
     do m=1,Total_block
       nx=bNi(m); ny=bNj(m); nz=bNk(m)
       nxe=nx+2*LAP-1; nye=ny+2*LAP-1; nze=nz+2*LAP-1
       NU=MP%NVAR*nxe*nye*nze; NT=nxe*nye*nze
       tag=B_n(m)+200
       if(B_proc(m) .eq. 0) then
         mt=B_n(m); B=>MP%Block(mt)
         read(99) nx1,ny1,nz1
         read(99) ((((B%U(m1,i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1),m1=1,MP%NVAR)
         read(99) (((B%Ts(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
         read(99) (((B%Tsn(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
         read(99) (((B%p(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
       else
         allocate(Ubuf(MP%NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Tsbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Tnbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         allocate(Pbuf(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
         read(99) nx1,ny1,nz1
         read(99) ((((Ubuf(m1,i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1),m1=1,MP%NVAR)
         read(99) (((Tsbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
         read(99) (((Tnbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
         read(99) (((Pbuf(i,j,k),i=1-LAP,nx+LAP-1),j=1-LAP,ny+LAP-1), &
                    k=1-LAP,nz+LAP-1)
         call MPI_Send(Ubuf, NU,OCFD_DATA_TYPE,B_proc(m),tag,  MPI_COMM_WORLD,ierr)
         call MPI_Send(Tsbuf,NT,OCFD_DATA_TYPE,B_proc(m),tag+1,MPI_COMM_WORLD,ierr)
         call MPI_Send(Tnbuf,NT,OCFD_DATA_TYPE,B_proc(m),tag+2,MPI_COMM_WORLD,ierr)
         call MPI_Send(Pbuf, NT,OCFD_DATA_TYPE,B_proc(m),tag+3,MPI_COMM_WORLD,ierr)
         deallocate(Ubuf,Tsbuf,Tnbuf,Pbuf)
       endif
     enddo
!    coupling-state trailer (iver = 2 only).  A missing/short trailer (e.g. a
!    file written by a killed run) simply leaves the state unknown -> the old
!    behaviour is kept.
     if(iver .ge. 2) then
       read(99,iostat=ios) Couple_State_Mode, Couple_State_Conv, &
                           Couple_State_Pair, Couple_State_Iter, &
                           Couple_State_nf, Couple_State_nk
       if(ios .eq. 0) then
         read(99,iostat=ios) Couple_State_Twmax, Couple_State_Tol
         if(ios .eq. 0 .and. Couple_State_nf .gt. 0 .and. &
            Couple_State_nk .gt. 0) then
           allocate(Couple_Tw_save(Couple_State_nf,Couple_State_nk))
           allocate(Couple_qw_save(Couple_State_nf,Couple_State_nk))
           allocate(Couple_u_save(Couple_State_nf,Couple_State_nk))
           allocate(Couple_pw_save(Couple_State_nf,Couple_State_nk))
           read(99,iostat=ios) Couple_Tw_save, Couple_qw_save, &
                               Couple_u_save, Couple_pw_save
           if(ios .ne. 0) then
             deallocate(Couple_Tw_save,Couple_qw_save,Couple_u_save,Couple_pw_save)
             Couple_State_nf=0; Couple_State_nk=0
           endif
         endif
         Couple_State_Found=1
         print*, " read_restart: coupling state  mode=", Couple_State_Mode, &
                 " (0=tight,1=staggered)  satisfied=", Couple_State_Conv, &
                 "  pair=", Couple_State_Pair, "  last max|dT_w|=", Couple_State_Twmax, &
                 "  interface cells=", Couple_State_nf, "x", Couple_State_nk
       else
         Couple_State_Found=0; Couple_State_Mode=-1
         if(Iflag_restart .eq. 1) iok=-1
         print*, " read_restart: no coupling-state trailer in the file"
       endif
     endif
     close(99)
   else
     do m=1,MP%Num_Block
       B=>MP%Block(m)
       nx=B%nx; ny=B%ny; nz=B%nz
       nxe=nx+2*LAP-1; nye=ny+2*LAP-1; nze=nz+2*LAP-1
       NU=MP%NVAR*nxe*nye*nze; NT=nxe*nye*nze
       tag=m+200
       call MPI_Recv(B%U,  NU,OCFD_DATA_TYPE,0,tag,  MPI_COMM_WORLD,status,ierr)
       call MPI_Recv(B%Ts, NT,OCFD_DATA_TYPE,0,tag+1,MPI_COMM_WORLD,status,ierr)
       call MPI_Recv(B%Tsn,NT,OCFD_DATA_TYPE,0,tag+2,MPI_COMM_WORLD,status,ierr)
       call MPI_Recv(B%p,  NT,OCFD_DATA_TYPE,0,tag+3,MPI_COMM_WORLD,status,ierr)
     enddo
   endif

!---- coupling state: rank 0 read it from the trailer, share it with all ranks -
   if(my_id .eq. 0 .and. Couple_State_Found .ne. 1) then
     Couple_State_Mode=-1; Couple_State_Conv=0; Couple_State_Pair=0
     Couple_State_Iter=0; Couple_State_nf=0; Couple_State_nk=0
     Couple_State_Twmax=0.d0; Couple_State_Tol=0.d0
   endif
   call MPI_Bcast(Couple_State_Found,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Mode, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Conv, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Pair, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Iter, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_nf,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_nk,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Twmax,1,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Couple_State_Tol,  1,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   if(Couple_State_Found .eq. 1 .and. Couple_State_nf .gt. 0 .and. &
      Couple_State_nk .gt. 0) then
     if(my_id .ne. 0) then
       allocate(Couple_Tw_save(Couple_State_nf,Couple_State_nk))
       allocate(Couple_qw_save(Couple_State_nf,Couple_State_nk))
       allocate(Couple_u_save(Couple_State_nf,Couple_State_nk))
       allocate(Couple_pw_save(Couple_State_nf,Couple_State_nk))
     endif
     call MPI_Bcast(Couple_Tw_save,Couple_State_nf*Couple_State_nk, &
                    OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Couple_qw_save,Couple_State_nf*Couple_State_nk, &
                    OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Couple_u_save, Couple_State_nf*Couple_State_nk, &
                    OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Couple_pw_save,Couple_State_nf*Couple_State_nk, &
                    OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   endif

!---- step/time bookkeeping + refresh the inter-block buffers -----------------
   Mesh(1)%Kstep = kst_read
   Mesh(1)%tt    = tt_read
   call MPI_Bcast(Mesh(1)%Kstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(Mesh(1)%tt,1,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   restart_found       = 1
   restart_Kstep_saved = kst_read
   restart_tt_saved    = tt_read
   call update_buffer_onemesh(1)
   call update_Ts_buffer_onemesh(1)
   ok=.true.
   if(my_id .eq. 0) then
     print*, " restart from '", trim(RESTART_FILE), "': Kstep=", Mesh(1)%Kstep, &
             "  tt=", Mesh(1)%tt, "  (flow/Ts/p + LAP buffer layers restored)"
     if(Num_Mesh .gt. 1) print*, " restart NOTE: only the finest mesh is restored"
   endif
  end subroutine read_restart

!------------------------------------------------------------------------------
! Per-save-step dispatcher: writes the restart file and (optionally) the
! node-centred SI flow field.  Called by ALL ranks right after the flow3d.dat
! output; it does its own Kstep-based due check so the call sites stay simple.
!------------------------------------------------------------------------------
  subroutine restart_step_output(nMesh)
   use Global_var
   implicit none
   integer,intent(in):: nMesh
   integer:: kst,kr
   kst = Mesh(nMesh)%Kstep
!  Iflag_restart = -1 disables the restart feature completely (neither read nor
!  written); 0 = auto (read when present, write periodically); 1 = force read.
   if(Iflag_restart .ge. 0) then
     kr  = Kstep_restart
     if(kr .le. 0) kr = Kstep_save
     if(kr .gt. 0) then
       if(mod(kst,kr) .eq. 0) call write_restart
     endif
   endif
   if(Iflag_flow_node .eq. 1 .and. Kstep_save .gt. 0) then
     if(mod(kst,Kstep_save) .eq. 0) call output_flow_node
   endif
  end subroutine restart_step_output

!------------------------------------------------------------------------------
! Register the coupling state that write_restart puts into the iver=2 trailer.
!   mode     : 0 = tight (per-step) coupling, 1 = staggered
!   conv     : 1 = the coupling already satisfied "tol" (Tol_Couple_Tw)
!   paircode : 11/12/13/19 (informational; validated when restoring)
!   iter     : outer iteration reached
!   twmax    : last convergence metric (max|dT_w| for 11/13; dp/du-type for 12;
!              max|dT_w| for 19)
!   tol      : tolerance "conv" was measured against
! The interface arrays are taken from the GLOBAL fp_Tw/fp_qw/fp_u/fp_pw (that is
! where the staggered drivers keep them); only arrays whose shape matches
! fp_Tw are copied, the rest stay zero, so the trailer layout is fixed.
! Call it at the end of every outer staggered iteration, right after the
! convergence metric is known, so that any periodic restart write carries a
! current state.
!------------------------------------------------------------------------------
  subroutine set_couple_state(mode, conv, paircode, iter, twmax, tol)
   use Global_var
   implicit none
   integer,intent(in):: mode, conv, paircode, iter
   real(PRE_EC),intent(in):: twmax, tol
   integer:: nf, nk

   Couple_State_Mode  = mode
   Couple_State_Conv  = conv
   Couple_State_Pair  = paircode
   Couple_State_Iter  = iter
   Couple_State_Twmax = twmax
   Couple_State_Tol   = tol
   Couple_State_Found = 1

   nf = 0; nk = 0
   if(allocated(fp_Tw)) then
     nf = size(fp_Tw,1); nk = size(fp_Tw,2)
   endif
   if(nf <= 0 .or. nk <= 0) then        ! no interface arrays in this driver (yet)
     Couple_State_nf = 0; Couple_State_nk = 0
     if(allocated(Couple_Tw_save)) deallocate(Couple_Tw_save,Couple_qw_save, &
                                              Couple_u_save,Couple_pw_save)
     return
   endif
   if(allocated(Couple_Tw_save)) then
     if(size(Couple_Tw_save,1) /= nf .or. size(Couple_Tw_save,2) /= nk) then
       deallocate(Couple_Tw_save,Couple_qw_save,Couple_u_save,Couple_pw_save)
     endif
   endif
   if(.not. allocated(Couple_Tw_save)) then
     allocate(Couple_Tw_save(nf,nk)); allocate(Couple_qw_save(nf,nk))
     allocate(Couple_u_save(nf,nk));  allocate(Couple_pw_save(nf,nk))
   endif
   Couple_Tw_save = 0.d0; Couple_qw_save = 0.d0
   Couple_u_save  = 0.d0; Couple_pw_save = 0.d0
   if(allocated(fp_Tw) .and. size(fp_Tw,1) == nf .and. size(fp_Tw,2) == nk) &
       Couple_Tw_save = fp_Tw
   if(allocated(fp_qw) .and. size(fp_qw,1) == nf .and. size(fp_qw,2) == nk) &
       Couple_qw_save = fp_qw
   if(allocated(fp_u)  .and. size(fp_u,1)  == nf .and. size(fp_u,2)  == nk) &
       Couple_u_save = fp_u
   if(allocated(fp_pw) .and. size(fp_pw,1) == nf .and. size(fp_pw,2) == nk) &
       Couple_pw_save = fp_pw
   Couple_State_nf = nf
   Couple_State_nk = nk
  end subroutine set_couple_state

!------------------------------------------------------------------------------
! Register "the staggered coupling is settled, we continue with the per-step
! (tight) coupling from now on".  Written into the trailer so that every later
! restart also stays tight.  No interface arrays are needed: in tight mode the
! interface quantities are recomputed from the restored U/Ts on every step, so
! the continued run has no jump by construction.
!------------------------------------------------------------------------------
  subroutine set_couple_tight_state()
   use Global_var
   implicit none
   Couple_State_Mode  = 0
   Couple_State_Conv  = 1
   Couple_State_Pair  = 0
   Couple_State_Iter  = 0
   Couple_State_Twmax = 0.d0
   Couple_State_Tol   = Tol_Couple_Tw
   Couple_State_Found = 1
   Couple_State_nf    = 0
   Couple_State_nk    = 0
   if(allocated(Couple_Tw_save)) deallocate(Couple_Tw_save,Couple_qw_save, &
                                            Couple_u_save,Couple_pw_save)
  end subroutine set_couple_tight_state


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
     write(99) 1, Total_block, Num_Mesh, NVAR1, Mesh(1)%Kstep, LAP
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
       else if(iver .ne. 1) then
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


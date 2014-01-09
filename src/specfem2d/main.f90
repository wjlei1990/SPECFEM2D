program wave2d

  use set_src_and_rec
  use wave2d_variables
  use wave2d_solver
  use mesh_io
  use wave2d_io_subs
  use simulation_type
  !use wave2d_write_seismo_subs

  use mpi

  implicit none
  !include "constants.h"

  !src and receiver
  integer :: ns, nrec
  integer, dimension(:), allocatable :: sglob, rglob
  double precision, dimension(:), allocatable :: x_src, z_src, x_rec, z_rec
  
  !time series
  double precision, dimension(:,:,:), allocatable :: samp, syn

  double precision :: f0(3)
  double precision :: source_time_function
  !character(len=200) :: snap_name, last_frame_name
  integer :: i
  double precision :: stf
  !double precision, dimension(:,:,:,:,:), allocatable :: absorb_field

  !Par_file var
  !integer :: NEX, NEZ, NELE, NPROC_XI, NPROC_ZI, NSTEP
  character(len=300) :: LOCAL_PATH, OUTPUT_PATH
  !double precision :: DENSITY, IMCOMPRESSIBILITY, RIGIDITY
  !logical :: DEBUG

  !MPI env var
  integer :: rank, nproc, comm
  integer :: ierr
  integer :: itime

  !********* PROGRAM STARTS HERE *********************
  call MPI_init(ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD,comm,ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nproc, ierr)

  call read_parfile(LOCAL_PATH, OUTPUT_PATH, rank, nproc, comm)
  stop

  !=================
  !setup
  !which type of simulation
  if(rank.eq.0)then
    print *, 'The type of simulation:',SIMUL_TYPE
    print *, "1)dynamic relax  2)wave propagation"
    call set_simulation_flag(SIMUL_TYPE,comm)
    ! set up model and mesher
    print *, "ALPHA=",S_ALPHA
    print *, "BETA=",S_BETA
  endif

  !===================
  !read mesh
  call read_mesh_init(rank,NEX,NEZ,nglob,nspec,nspecb,ninterface)
  NELE=max(NEX, NEZ)
  allocate(x(nglob))
  allocate(z(nglob))
  allocate(ibool(NGLLX,NGLLZ,nspec))
  allocate(ibelm(4,NELE))
  allocate(my_neighbours(8))
  allocate(nibool_interfaces(8))
  allocate(ibool_interfaces(8,NELE))
  
  call read_array(rank,NEX,NEZ,nglob,nspec,nspecb,ninterface,&
        x,z,ibool,ibelm,my_neighbours,nibool_interfaces,ibool_interfaces)
  !call set_model_property()
  !stop

  call read_src_and_rec(ns, x_src, z_src, nrec, x_rec, z_rec, &
        rank, nproc, comm)
  call locate_src_and_rec(x, z, ibool, ns, sglob, x_src, z_src, &
        nrec, rglob, x_rec, z_rec, rank, nproc, comm)

  !===================
  !locate source and receiver
  allocate(samp(NSTEP,3,ns))
  samp(:,:,:) = 0.0
  allocate(syn(NSTEP,3,nrec))
  syn(:,:,:) = 0.0
  ! setup source time function
  print *, 'Setup source time functions ...'
  hdur = 100 * DT 
  print *, 'hdur = ',sngl(hdur),' s'
  !f0(1) = 0.d10; f0(2) = 1.0d10; f0(3) = 0.d10 ! for SH source
  f0(1) = 1.d10; f0(2) = 0.0d10; f0(3) = 1.d10  ! for SV source
  samp = 0.0
  do itime = 1, NSTEP
    stf = source_time_function(dble(itime-1)*DT-hdur,hdur)
    do i = 1, ns
      samp(itime, :, i) = stf * f0(:)  
    enddo
  enddo
  
  !=========================
  ! solver for forward wavefield
  print *, 'Start simulation ...'
  call solver(ns, sglob, samp, nrec, rglob, syn, &
    rank, nproc, comm)

  print *, 'Write out seismograms at the receiver ...'
  call write_seismogram(syn, rglob, nrec, DT, NSTEP, OUTPUT_PATH)

  ! deallocate variables
  deallocate(sglob,samp,x_src,z_src,rglob,syn,x_rec,z_rec)

  call MPI_Barrier(comm,ierr)
  call MPI_Finalize(ierr)

end program wave2d
 

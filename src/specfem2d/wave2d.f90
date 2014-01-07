program wave2d

  use wave2d_variables
  use wave2d_solver
  use wave2d_sub
  use simulation_type

  implicit none

  integer :: ns, nrec
  double precision, dimension(:,:,:), allocatable :: samp, data, syn, adj_syn
  integer, dimension(:), allocatable :: sglob, rglob
  double precision, dimension(:), allocatable :: x_source, z_source, x_rec, z_rec

  double precision :: t0, f0(3), ft(3), pxpf(3), f1(3), f2(3), dxdf(3)
  double precision :: p_1, p_2, df, phi
  double precision :: source_time_function
  character(len=200) :: snap_name, last_frame_name
  integer :: i, index, l, irec, itime, icomp, iter, isave, isolver
  double precision :: xold(3), fold, gold(3), pold(3), xnew(3), fnew, gnew(3), pnew(3)
  double precision :: lam, beta, tstart, tend, stf
  double precision, dimension(:,:,:,:,:), allocatable :: absorb_field
  double precision, dimension(:), allocatable :: rho_kernel, mu_kernel, kappa_kernel

  !********* PROGRAM STARTS HERE *********************
  call MPI_init(ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD,comm,ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, nproc, ierr)

  !=================
  !setup
  out_dir = "OUTPUT_FILES/"
  !which type of simulation
  SIMUL_TYPE=2
  IM_TRUE=.false.
  print *, 'The type of simulation:',SIMUL_TYPE
  print *, "1)dynamic relax  2)wave propagation"
  call set_simulation_flag()
  ! set up model and mesher
  print *, "ALPHA=",S_ALPHA
  print *, "BETA=",S_BETA

  !===================
  !read mesh
  call read_mesh()
  !call set_model_property()
  !stop

  call read_src_and_rec(ns, x_src, z_src,&
        nrec, x_rec, z_rec, rank, nproc, comm)
  call locate_src_and_rec(ns, sglob, x_src, z_src, &
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
  call write_seismogram(syn, nrec, trim(out_dir)//'seismograms')

  ! deallocate variables
  deallocate(sglob,samp,x_source,z_source,rglob,syn,x_rec,z_rec)

  call MPI_Barrier(comm,ierr)
  call MPI_Finalize(ierr)

end program wave2d
 

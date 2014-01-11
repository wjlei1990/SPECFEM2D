program wave2d_kernel

  use wave2d_variables
  use wave2d_solver
  use wave2d_sub

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
  double precision :: lam, beta, eps, tstart, tend, stf
  double precision, dimension(:,:,:,:,:), allocatable :: absorb_field
  double precision, dimension(:), allocatable :: rho_kernel, mu_kernel, kappa_kernel

    
  !********* PROGRAM STARTS HERE *********************
  
  out_dir = "OUTPUT_FILES/"

  ! set up model and mesher
  call mesher()
  call set_model_property()

  ! allocate memory for source and receiver related variables
  ns = 1
  allocate(sglob(ns))
  allocate(x_source(ns),z_source(ns))
  allocate(samp(NSTEP,3,ns))

  nrec = 1
  allocate(rglob(nrec))
  allocate(x_rec(nrec),z_rec(nrec))
  allocate(syn(NSTEP,3,nrec), data(NSTEP,3,nrec), adj_syn(NSTEP,3,nrec))

  ! initialization
  sglob = 0;  rglob = 0
  samp(:,:,:) = 0.0;  syn(:,:,:) = 0.0

  ! set up source and receiver geometry
  x_source(1) = LENGTH/4; z_source(1) = HEIGHT/2
  x_rec(1) = 3 * LENGTH/4; z_rec(1) = HEIGHT/2
  call set_glob(nrec, x_rec, z_rec, rglob)
  call set_glob(ns,x_source,z_source,sglob)

  open(12,file=trim(out_dir)//'sr.txt',status='unknown',iostat=ios)
  if (ios /= 0) stop 'Error openning out_dir/sr.txt'
  do i = 1, ns
    write(12,*) 'S ', sngl(x_source(i)/LENGTH), sngl(z_source(i)/LENGTH)
  enddo
  do i = 1, nrec
    write(12,*) 'R ', sngl(x_rec(i)/LENGTH), sngl(z_rec(i)/LENGTH)
  enddo
  close(12)

  hdur = 100 * DT 
  print *, 'hdur = ',sngl(hdur),' s'

  ! solve for forward wavefield
  print *, 'Calculate forward wavefield'
  f0(1) = 1.d10; f0(2) = 0.0d10; f0(3) = 0

  samp = 0.0
  do itime = 1, NSTEP
    stf = source_time_function(dble(itime-1)*DT-hdur,hdur)
    do i = 1, ns
      samp(itime, :, i) = stf * f0(:)  
    enddo
  enddo

! ************** following are different for different types of simulation *******************
  isolver = 1; last_frame_name = trim(out_dir)//'last_frame.txt'
  allocate(absorb_field(NSTEP, 3, NGLL, NELE, 3))
  call solver(isolver, ns, sglob, samp, nrec, rglob, syn, trim(last_frame_name), absorb_field)

  print *, 'Write out seismograms at the receiver'
  call write_seismogram(syn, nrec, trim(out_dir)//'forward')

! choose P
!  tstart = 15.0; tend = 22.0; irec = 1
! choose PP
!  tstart = 22.0; tend = 27.0; irec = 1
! choose PS+SP
!  tstart = 28.0; tend = 34.0; irec = 1
! choose S
!  tstart = 30.0; tend = 35.0; irec = 2
! choose SS
  tstart = 37.0; tend = 47.0; irec = 1
  call traveltime_adjoint_source(nrec, irec, syn, tstart, tend, adj_syn)
   
  ! kernels
  isolver = 3
  allocate(rho_kernel(NGLOB), mu_kernel(NGLOB), kappa_kernel(NGLOB))
  call solver(isolver, ns, sglob, samp, nrec, rglob, adj_syn, trim(last_frame_name), absorb_field, &
       rho_kernel, mu_kernel, kappa_kernel)

  ! deallocate variables
  deallocate(sglob,samp,x_source,z_source,rglob,syn,data,adj_syn,x_rec,z_rec)
  deallocate(absorb_field, rho_kernel, mu_kernel,  kappa_kernel)

end program wave2d_kernel
 !adjoint field
  !  call solver(2, ns,sglob,samp, nrec, rglob, syn)
  !  print *, 'write out adjoint seismogram at the source '
  !  call write_seismogram(samp, ns, 'adjoint')
 

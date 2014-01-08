program main
!mesh main program

  !use wave2d_variables
  use wave2d_mesh

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

  !output dir of mesh file
  out_dir = "OUTPUT_FILES/"

  !==========
  !start mesh
  read_parameters(NEX,NEZ)

  nex_one_proc=NEX/NPROC_X
  nez_one_proc=NEZ/NPROC_Z
  nspec=nex_one_proc*nez_one_proc
  call acllocate_array()
  call get_ibool()
  call calcu_xigll_zgll(xigll,zigll)

!do each processor

  do i=1, NRPOC_Z
    do j=1, NPROC_X
       iproc=(i-1)*NPROC_X
       mesh_anchor(1,1)=MODEL_X1
       mesh_anchor(1,2)=MODEL_Z1
       mesh_anchor(2,1)=mesh_anchor+LENGTH
       mesh_anchor(2,2)=mesh_anchor+HEIGHT
       call init_variables(nspec,x,z,max_nex_nez,ibelm,niboll_interfaces,ibool_interfaces,rho,&
              kappa,mu)
       call mesh_one_proc(SIMUL_TYPE,mesh_anchor,nex_one_proc,nez_one_proc,nspec,LENGTH,HEIGHT,&
              x,z,nglob,ibool,xigll,zigll,ibelm,MODEL_X1,MODEL_X2,&
              MODEL_Z1,MODEL_Z2,ID,DEBUG)

       call build_mesh_interface(nex_one_proc,nez_one_proc,anchor, ibool, max_interface_size, ninterface,my_neighbour, &
              nibool_interfaces,ibool_interfaces,anchor,MODEL_X1,MODEL_X2,MODEL_Z1,&
              MODEL_Z2,NPROC_X)
       call set_model_property(NSPEC,nglob,ibool,x,z,rho,kappa,mu)

       call save_mesh(iproc,nex_one_proc,nez_one_proc,nglob,ibool,x,z,nspec,nspecb,ibelm, &
              ninterface,my_neighbour,nibool_interfaces,ibool_interfaces,ID,&
              rho,kappa,mu)
    end do
  enddo


  !which type of simulation
  !SIMUL_TYPE=2
  !IM_TRUE=.false.
  !print *, 'The type of simulation:',SIMUL_TYPE
  !print *, "1)dynamic relax  2)wave propagation"
  !call set_simulation_flag()
  ! set up model and mesher
  !print *, "ALPHA=",S_ALPHA
  !print *, "BETA=",S_BETA
  
  !print *, 'Allocate source and receivers ...'
  ! sources
  !ns = 1
  !allocate(sglob(ns))
  !allocate(x_source(ns),z_source(ns))
  !allocate(samp(NSTEP,3,ns))
  !sglob = 0; samp(:,:,:) = 0.0
  !x_source(1) = LENGTH/4.; z_source(1) = 3./4.*HEIGHT
  !call set_glob(ns,x_source,z_source,sglob)
  !print *,"x_source,z_source,"

  ! receivers
  !nrec = 3
  !allocate(rglob(nrec))
  !allocate(x_rec(nrec),z_rec(nrec))
  !allocate(syn(NSTEP,3,nrec))
  !rglob(:) = 0; syn(:,:,:) = 0.0
  !x_rec(1) = 1 * LENGTH/4; z_rec(1) = HEIGHT
  !x_rec(2) = 2 * LENGTH/4; z_rec(2) = HEIGHT
  !x_rec(3) = 3 * LENGTH/4; z_rec(3) = HEIGHT
  !call set_glob(nrec, x_rec, z_rec, rglob)

  ! output source and receiver location xy files
  !open(12,file=trim(out_dir)//'sr.txt',status='unknown',iostat=ios)
  !if (ios /= 0) stop 'Error openning out_dir/sr.txt'
  !do i = 1, ns
  !  write(12,*) 'S ', sngl(x_source(i)/LENGTH), sngl(z_source(i)/HEIGHT)
  !enddo
  !do i = 1, nrec
  !  write(12,*) 'R ', sngl(x_rec(i)/LENGTH), sngl(z_rec(i)/HEIGHT)
  !enddo
  !close(12)

  ! setup source time function
  !print *, 'Setup source time functions ...'
  !hdur = 100 * DT 
  !print *, 'hdur = ',sngl(hdur),' s'

   !f0(1) = 0.d10; f0(2) = 1.0d10; f0(3) = 0.d10 ! for SH source
   !f0(1) = 1.d10; f0(2) = 0.0d10; f0(3) = 1.d10  ! for SV source

  !samp = 0.0
  !do itime = 1, NSTEP
  !  stf = source_time_function(dble(itime-1)*DT-hdur,hdur)
  !  do i = 1, ns
  !    samp(itime, :, i) = stf * f0(:)  
  !  enddo
  !enddo

  
  ! solve for forward wavefield
  !print *, 'Start simulation ...'
  !call solver(ns, sglob, samp, nrec, rglob, syn)

  !print *, 'Write out seismograms at the receiver ...'
  !call write_seismogram(syn, nrec, trim(out_dir)//'seismograms')

  ! deallocate variables
  !deallocate(sglob,samp,x_source,z_source,rglob,syn,x_rec,z_rec)

end program main 
 

subroutine read_parameters(nex,nez,model_x1,model_x2,model_z1,model_z2,nproc_x,nproc_z,&
        simul_type,debug,density0,incompressibility0,rigidity0)

   implicit none

   integer ::nex,nez
   double precision ::model_x1,model_x2,model_z1,model_z2
   integer ::nproc_x,nproc_z
   integer ::simul_type
   logical ::debug
   double precision::density0,incompressibility0,rigidity0

   !local variables or not used in mesh

  integer :: idummy=5, i, ierr
  character(len=30) :: dstring
  integer :: IIN=110

  character(len=10) :: model_type
  character(len=200)::local_path,output_path
  double precision  ::record_length_in_seconds,dt

!   nex=8
!   nez=8
!   nproc_x=2
!   nproc_z=2
!   model_x1=0
!   model_x2=8000
!   model_z1=0
!   model_z2=8000
!   simul_type=2
!   debug=.false.

  open(unit=IIN,file="DATA/Par_file")

  do i=1,idummy
    read(IIN,*)
  enddo
  read(IIN,5) dstring, simul_type
  read(IIN,*)
  read(IIN,*)
  read(IIN,5) dstring, nex
  read(IIN,5) dstring, nez
  read(IIN,*)
  read(IIN,*)
  read(IIN,5) dstring, nproc_x
  read(IIN,5) dstring, nproc_z
  read(IIN,*)
  read(IIN,*)
  read(IIN,6) dstring, local_path
  read(IIN,*)
  read(IIN,*)
  read(IIN,6) dstring, output_path
  read(IIN,*)
  read(IIN,*)
  read(IIN,7) dstring, record_length_in_seconds
  read(IIN,*)
  read(IIN,*)
  read(IIN,7) dstring, dt
  read(IIN,*)
  read(IIN,*)
  read(IIN,7) dstring, model_x1
  read(IIN,7) dstring, model_x2
  read(IIN,7) dstring, model_z1
  read(IIN,7) dstring, model_z2
  read(IIN,*)
  read(IIN,*)
  read(IIN,6) dstring, model_type
  read(IIN,7) dstring, density0
  read(IIN,7) dstring, incompressibility0
  read(IIN,7) dstring, rigidity0
  read(IIN,*)
  read(IIN,*)
  read(IIN,8) dstring, debug
  close(IIN)


  print  *,'nex',nex,nez,nproc_x,nproc_z,model_x1,model_x2,model_z1,model_z2,simul_type
!   nex=8
!   nez=8
!   nproc_x=2
!   nproc_z=2
!   model_x1=0
!   model_x2=8000
!   model_z1=0
!   model_z2=8000
!   simul_type=2

5 format(A,I5)
6 format(A,A)
7 format(A,E15.5)
8 format(A,L5)

end subroutine read_parameters

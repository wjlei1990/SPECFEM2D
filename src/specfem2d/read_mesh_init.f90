subroutine read_mesh_init(iproc,local_path)
  
  use wave2d_variables
  use constants
  implicit none
 
  integer ::iproc
  character(len=*) ::local_path
!  integer ::NEX,NEZ
!  integer ::nglob
!  integer, dimension(NGLLX,NGLLZ,nspec)::ibool
!  double precision,dimension(nglob) ::x,z
!  integer ::nspec
!  integer ::max_interface_size
!  integer, dimension(4) ::nspecb
!  integer,dimension(4,max_interface_size) ::ibelm
!  integer ::ninterface
!  integer, dimension(ninterface) ::my_neighbour,nibool_interfaces
!  integer ::max_interface_size_node
!  integer,dimension(max_interface_size_node,8) ::ibool_interfaces
!  integer,dimension(3,nglob) ::id
!  double precision,dimension(NGLLX,NGLLZ,nspec) ::rho,kappa,mu

  !local variable
  character(len=200) ::outfile
  integer ::ios
  

  write(*,*) 'save the database'
  write(outfile,"(('DATABASES_MPI/MESHINFO_',i5.5))") iproc
  open(iproc,FILE=trim(outfile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  read(iproc,*) NEX,NEZ
  read(iproc,*) nglob
  read(iproc,*) nspec
  read(iproc,*) nspecb(:)
  read(iproc,*) NELE
  read(iproc,*) ninterface
  read(iproc,*) max_ibool_interface_size
  close(iproc)
end subroutine read_mesh_init

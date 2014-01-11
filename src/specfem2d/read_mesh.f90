!mesh io file
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

  !write(*,*) 'save the database'
  write(outfile,"(('MESHINFO_',i5.5))") iproc
  !print *, "outfile:", trim(outfile)
  outfile=trim(local_path)//'/'//trim(outfile)
  !print *, "outfile:", trim(outfile)
  open(iproc,FILE=trim(outfile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  read(iproc,*) NEX,NEZ
  read(iproc,*) nglob
  read(iproc,*) nspec
  read(iproc,*) nspecb(:)
  read(iproc,*) NELE
  read(iproc,*) ninterface
  read(iproc,*) max_ibool_interfaces_size
  close(iproc)
end subroutine read_mesh_init


subroutine read_mesh(iproc,local_path)
  
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
  integer ::iglob,ispec
  integer ::i,j
  

  
  write(*,*) 'read the database:',iproc
  write(outfile,"(('MESHINFO_',i5.5))")iproc
  outfile=trim(local_path)//'/'//trim(outfile)
  open(iproc,FILE=trim(outfile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  read(iproc,*) NEX,NEZ
  read(iproc,*) nglob
  read(iproc,*) nspec
  read(iproc,*) nspecb(:)
  read(iproc,*) NELE
  read(iproc,*) ninterface
  read(iproc,*) max_ibool_interfaces_size
  N_EQ=3*nglob

  do iglob=1,nglob
    read(iproc,*) x(iglob),z(iglob)
  end do

  do ispec=1,nspec
    do j=1,NGLLZ
      read(iproc,*) ibool(:,j,ispec)
    end do
  end do

  do iglob=1,nglob
    read(iproc,*)id(:,iglob)
  end do

  do ispec=1,nspec
    do i=1,NGLLZ
      do j=1,NGLLX
        read(iproc,*)rho(i,j,ispec),kappa(i,j,ispec),mu(i,j,ispec)
      end do
    end do
  end do

  do i=1,4
    read(iproc,*) ibelm(:,i)
  end do
  do i=1,ninterface
    read(iproc,*)my_neighbours(i),nibool_interfaces(i)
    read(iproc,*)ibool_interfaces(:,i)
  end do
  close(iproc)
end subroutine read_mesh

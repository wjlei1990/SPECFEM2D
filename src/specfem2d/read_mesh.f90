module mesh_io

 use constants
 implicit none

contains

subroutine read_mesh_init(iproc,NEX,NEZ,nglob,nspec,nspecb, &
    ninterface)

  integer ::iproc
  integer ::NEX,NEZ
  integer ::nglob
  integer ::nspec
  integer, dimension(4) ::nspecb
  integer ::ninterface

  !local variable
  character(len=200) ::infile
  integer ::ios
  

  write(*,*) 'read the database'
  write(infile,"(('MESHINFO',i5.5))") iproc
  open(iproc,FILE=trim(infile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  read(iproc,*) NEX,NEZ
  read(iproc,*) nglob
  read(iproc,*) nspec
  write(iproc,*) nspecb(:)
  write(iproc,*)ninterface
end subroutine read_mesh_init
 
subroutine read_array(iproc,NEX,NEZ,nglob,nspec,nspecb, &
    ninterface,x,z,ibool,ibelm,my_neighbour,nibool_interfaces,&
    ibool_interfaces)

  integer ::iproc
  integer ::NEX,NEZ
  integer ::nglob
  integer, dimension(NGLLX,NGLLZ,nspec)::ibool
  double precision,dimension(nglob) ::x,z
  integer ::nspec
  integer, dimension(4) ::nspecb
  integer,dimension(:,:) ::ibelm
  integer ::ninterface
  integer, dimension(ninterface) ::my_neighbour,nibool_interfaces
  integer,dimension(:,:) ::ibool_interfaces

  !local variable
  integer ::iglob,ispec
  integer ::i,j

  do iglob=1,nglob
    read(iproc,*) x(iglob),z(iglob)
  end do
  do ispec=1,nspec
    do j=1,NGLLZ
      read(iproc,*) ibool(:,j,ispec)
    end do
  end do
  do i=1,4
    read(iproc,*) ibelm(i,:)
  end do
  do i=1,ninterface
    read(iproc,*)my_neighbour(i),nibool_interfaces(i)
    read(iproc,*)ibool_interfaces(i,1:nibool_interfaces(i))
  end do
  close(iproc)
end subroutine read_array

end module mesh_io

subroutine save_mesh(iproc,NEX,NEZ,nglob,ibool,x,z,nspec,nspecb,max_interface_size,ibelm, &
    ninterface,my_neighbour,nibool_interfaces,max_interface_size_node,ibool_interfaces,ID,&
    rho,kappa,mu)
  
  use constants
  implicit none
 
  integer ::iproc
  integer ::NEX,NEZ
  integer ::nglob
  integer, dimension(NGLLX,NGLLZ,nspec)::ibool
  double precision,dimension(nglob) ::x,z
  integer ::nspec
  integer ::max_interface_size
  integer, dimension(4) ::nspecb
  integer,dimension(4,max_interface_size) ::ibelm
  integer ::ninterface
  integer, dimension(ninterface) ::my_neighbour,nibool_interfaces
  integer ::max_interface_size_node
  integer,dimension(max_interface_size_node,8) ::ibool_interfaces
  integer,dimension(3,nglob) ::ID
  double precision,dimension(NGLLX,NGLLZ,nspec) ::rho,kappa,mu

  !local variable
  character(len=200) ::outfile
  integer ::ios
  integer ::iglob,ispec
  integer ::i,j
  

  write(*,*) 'save the database'
  write(outfile,"(('OUTPUT_FILES/MESHINFO_',i5.5))") iproc
  open(iproc,FILE=trim(outfile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  write(iproc,*) NEX,NEZ
  write(iproc,*) nglob
  write(iproc,*) nspec
  write(iproc,*) nspecb(:)
  write(iproc,*) max_interface_size
  write(iproc,*) ninterface
  write(iproc,*) max_interface_size_node


  do iglob=1,nglob
    write(iproc,*) x(iglob),z(iglob)
  end do

  do ispec=1,nspec
    do j=1,NGLLZ
      write(iproc,*) ibool(:,j,ispec)
    end do
  end do

  do iglob=1,nglob
    write(iproc,*)id(:,iglob)
  end do

  do ispec=1,nspec
    do i=1,NGLLZ
      do j=1,NGLLX
        write(iproc,*)rho(i,j,ispec),kappa(i,j,ispec),mu(i,j,ispec)
      end do
    end do
  end do

  do i=1,4
    write(iproc,*) ibelm(i,:)
  end do
  do i=1,ninterface
    write(iproc,*)my_neighbour(i),nibool_interfaces(i),i
    write(iproc,*)ibool_interfaces(:,i)
  end do
  close(iproc)
end subroutine save_mesh

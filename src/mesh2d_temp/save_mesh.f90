subroutine save_mesh(iproc,NEX,NEZ,nglob,ibool,x,z,nspec,nspecb,ibelm, &
    ninterface,my_neighbour,nibool_interfaces,ibool_interfaces)

  integer ::iproc
  integer ::NEX,NEZ
  integer ::nglob
  integer, dimension(NGLLX,NGLLZ,nspec)::ibool
  double precision,dimension(nglob) ::x,z
  integer ::nspec
  integer, dimension(4) ::nsepcb
  integer,dimension(:,:), intent(in) ::ibelm
  integer ::ninterface
  integer, dimension(ninterface) ::my_neighbour,nibool_interfaces
  integer,dimension(:,:),intent(in) ::ibool_interfaces

  !local variable
  character(len=200) ::outfile
  integer ::ios
  integer ::iglob,ispec
  integer ::i,j
  

  write(*,*) 'save the database'
  write(outfile,"(('MESHINFO',i5.5))") iproc
  open(iproc,FILE=trim(outfile),status='unknown',iostat=ios)
  if(ios /= 0) stop 'Error in openning the database file'
  write(iproc,*) NEX,NEZ
  write(iproc,*) nglob
  write(iproc,*) nspec
  write(iproc,*) nspecb(:)
  write(iproc,*)ninterface


  do iglob=1,nglob
    write(iproc,*) x(iglob),z(iglob)
  end do
  do ispec=1,nspec
    do j=1,NGLLZ
      write(iproc,*) ibool(:,j,ispec)
    end do
  end do
  do i=1,4
    write(iproc,*) ibelm(i,:)
  end do
  do i=1,ninterface
    write(iproc,*)my_neighbour(i),nibool_interfaces(i)
    write(iproc,*)ibool_interfaces(i,1:nibool_interfaces(i))
  end do
  close(iproc)
end subroutine save_mesh

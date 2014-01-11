module build_mesh_interface_subs

  use constants

contains

  subroutine build_mesh_interface(iproc,NEX,NEZ,anchor, ibool, max_interface_size, my_neighbour, &
      nibool_interfaces, ibool_interfaces,model_x1,model_x2,model_z1,&
      model_z2,nproc_x)
    implicit none

    integer ::iproc
    double precision,intent(in) ::model_x1,model_x2,model_z1,model_z2
    double precision,dimension(:,:), intent(in) :: anchor
    integer, intent(in) ::NEX,NEZ
    integer, intent(in) ::nproc_x
    integer, dimension(:,:,:), intent(in) ::ibool

    integer, dimension(MAX_INTERFACE) :: my_neighbour
    integer, dimension(MAX_INTERFACE) :: nibool_interfaces
    integer ::max_interface_size
    integer, dimension(max_interface_size, MAX_INTERFACE),intent(out) :: ibool_interfaces

    logical :: boundary_bool(4)
    integer :: i,j,inode

    boundary_bool(:)=.true. 
    nibool_interfaces=0

    !calculate the number of interfaces 
    if(abs(anchor(1,1)-model_x1).gt.TOL)then
      boundary_bool(1)=.false.
      !n_model_boundary=n_model_boundary+1
    endif
    if(abs(anchor(1,2)-model_x2).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(3)=.false.
    endif
    if(abs(anchor(2,1)-model_z1).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(2)=.false.
    endif
    if(abs(anchor(2,2)-model_z2).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(4)=.false.
    endif
    print *,'anchor',iproc,'anchor',anchor(:,1),anchor(:,2)
    my_neighbour(1)=iproc-1
    my_neighbour(2)=iproc-nproc_x-1
    my_neighbour(3)=iproc-nproc_x
    my_neighbour(4)=iproc-nproc_x+1
    my_neighbour(5)=iproc+1
    my_neighbour(6)=iproc+nproc_x+1
    my_neighbour(7)=iproc+nproc_x
    my_neighbour(8)=iproc+nproc_x-1

    if(.not.boundary_bool(1)) then
      nibool_interfaces(1)=NEZ*NGLLZ-(NEZ-1)
      inode=0
      do i=1,NEZ
        do j=1,NGLLZ-1
          inode=inode+1
          ibool_interfaces(inode,1)=ibool(1,j,NEX*(i-1)+1)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,1)=ibool(1,NGLLZ,NEX*(NEZ-1)+1)
      if(.not.boundary_bool(2)) then
        nibool_interfaces(2)=1
        ibool_interfaces(1,2)=ibool(1,1,1)
      end if
    end if

    
    if(.not.boundary_bool(2)) then
      nibool_interfaces(3)=NEX*NGLLX-(NEX-1)
      inode=0
      do i=1,NEX
        do j=1,NGLLX-1
          inode=inode+1
          ibool_interfaces(inode,3)=ibool(j,1,i)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,3)=ibool(NGLLX,1,NEX)
      if(.not.boundary_bool(3)) then
        nibool_interfaces(4)=1
        ibool_interfaces(1,4)=ibool(NGLLX,1,NEX)
      end if
    end if

      
    if(.not.boundary_bool(3)) then
      nibool_interfaces(5)=NEZ*NGLLZ-(NEZ-1)
      inode=0
      do i=1,NEZ
        do j=1,NGLLZ-1
          inode=inode+1
          ibool_interfaces(inode,5)=ibool(NGLLX,j,NEX*i)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,5)=ibool(NGLLX,NGLLZ,NEX*NEZ)
      if(.not.boundary_bool(4)) then
        nibool_interfaces(6)=1
        ibool_interfaces(1,6)=ibool(NGLLX,NGLLZ,NEX*NEZ)
      end if
    end if
        

    if(.not.boundary_bool(4)) then
      nibool_interfaces(7)=NEX*NGLLX-(NEX-1)
      inode=0
      do i=1,NEX
        do j=1,NGLLX-1
          inode=inode+1
          ibool_interfaces(inode,7)=ibool(j,NGLLZ,NEZ*(NEX-1)+i)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,7)=ibool(NGLLX,NGLLZ,NEZ*NEX)
      if(.not.boundary_bool(1)) then
        nibool_interfaces(8)=1
        ibool_interfaces(1,8)=ibool(1,NGLLZ,NEX*(NEZ-1)+1)
      end if
    end if
      
    do i=1,8
      if(my_neighbour(i)<0) nibool_interfaces(i)=0
    end do
    


    !if(n_model_boundary.eq.1)then
    !  ninterface=5
    !else if(n_model_boundary.eq.2)then
    !  ninterface=3
    !else if(n_model_boundary.eq.0)then
    !  ninterface=8
    !else
    !  print *, "Error in calculating ninterface"
    !  stop
    !endif

  end subroutine build_mesh_interface




end module build_mesh_interface_subs

module build_mesh_interface_subs

  use constants

contains

  subroutine build_mesh_interface(NEX,NEZ,anchor, ibool, max_interface_size, ninterface, my_neighbour, &
      nibool_interfaces, ibool_interfaces)

    double precision, intent(in) :: anchor(:)
    integer, intent(in) ::NEX,NEZ
    integer, dimension(:,:,:), intent(in) ::ibool

    integer :: ninterface
    integer, dimension(MAX_INTERFACE) :: my_neighbour
    integer, dimension(MAX_INTERFACE) :: nibool_interfaces
    integer, dimension(max_interface_size, MAX_INTERFACE) :: ibool_interfaces

    logical :: boundary_bool(4)
    integer :: i,j 

    boundary_bool(:)=.true. 
    nibool_interfaces=0
    ninterface=8

    !calculate the number of interfaces 
    if(abs(anchor(1,1)-MODEL_X1).gt.TOL)then
      boundary_bool(1)=.false.
      !n_model_boundary=n_model_boundary+1
    endif
    if(abs(anchor(1,2)-MODEL_X2).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(3)=.false.
    endif
    if(abs(anchor(2,1)-MODEL_Z1).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(2)=.false.
    endif
    if(abs(anchor(2,4)-MODEL_Z2).gt.TOL)then
      !n_model_boundary=n_model_boundary+1
      boundary_bool(4)=.false.
    endif
    my_neighbour(1)=iproc-1
    my_neighbour(2)=iproc-NPROC_X-1
    my_neighbour(3)=iproc-NPROC_X
    my_neighbour(4)=iproc-NPROC_X+1
    my_neighbour(5)=iproc+1
    my_neighbour(6)=iproc+NPROC_X+1
    my_neighbour(7)=iproc+NPROC_X
    my_neighbour(8)=iproc+NPROC_X-1

    if(.not.boundary_bool(1)) then
      nbool_interfaces(1)=NEZ*NGLLZ-(NEZ-1)
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
        nbool_interfaces(2)=1
        ibool_interfaces(1,2)=ibool(1,1,1)
      end if
    end if

    
    if(.not.boundary_bool(2)) then
      nbool_interfaces(3)=NEX*NGLLX-(NEX-1)
      inode=0
      do i=1,NEX
        do j=1,NGLLX-1
          inode=inode+1
          ibool_interfaces(inode,1)=ibool(1,j,NEZ*(i-1)+1)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,1)=ibool(1,NGLLX,NEZ*(NEX-1)+1)
      if(.not.boundary_bool(3)) then
        nbool_interfaces(4)=1
        ibool_interfaces(1,4)=ibool(NGLLX,1,NEX)
      end if
    end if

      
    if(.not.boundary_bool(3)) then
      nbool_interfaces(5)=NEZ*NGLLZ-(NEZ-1)
      inode=0
      do i=1,NEZ
        do j=1,NGLLZ-1
          inode=inode+1
          ibool_interfaces(inode,1)=ibool(1,j,NEX*(i-1)+1)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,1)=ibool(1,NGLLZ,NEX*(NEZ-1)+1)
      if(.not.boundary_bool(4)) then
        nbool_interfaces(6)=1
        ibool_interfaces(1,6)=ibool(NGLLX,NGLLZ,NEX*NEZ)
      end if
    end if
        

    if(.not.boundary_bool(4)) then
      nbool_interfaces(7)=NEX*NGLLX-(NEX-1)
      inode=0
      do i=1,NEX
        do j=1,NGLLX-1
          inode=inode+1
          ibool_interfaces(inode,1)=ibool(1,j,NEZ*(i-1)+1)         
        end do
      end do
      inode=inode+1
      ibool_interfaces(inode,1)=ibool(1,NGLLX,NEZ*(NEX-1)+1)
      if(.not.boundary_bool(1)) then
        nbool_interfaces(8)=1
        ibool_interfaces(1,8)=ibool(1,NGLLZ,NEX*(NEZ-1)+1)
      end if
    end if
      
    


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

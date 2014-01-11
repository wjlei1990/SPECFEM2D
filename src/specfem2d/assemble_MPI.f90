
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
! These subroutines are for the most part not used in the sequential version.
!

module assemble_subs

!#ifdef USE_MPI

!-----------------------------------------------
! Assembling the mass matrix.
!-----------------------------------------------
contains

  subroutine assemble_MPI_mass(array_val,npoin_val, &
              ninterface,  max_ibool_interfaces_size, &
              ibool_interfaces, nibool_interfaces, my_neighbours,&
              rank, comm)

  use mpi
  use constants

  implicit none

  !include 'constants.h'

  integer :: npoin_val
  double precision, dimension(npoin_val), intent(inout) :: array_val

  integer, intent(in)  :: ninterface, max_ibool_interfaces_size
  integer, dimension(max_ibool_interfaces_size, ninterface), intent(in)  :: &
    ibool_interfaces
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces
  integer, dimension(ninterface), intent(in)  :: my_neighbours
  ! array to assemble
  ! acoustic

  integer  :: ipoin, num_interface
  integer  :: ier,ierr
  integer  :: i
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
  double precision, dimension(max_ibool_interfaces_size, ninterface) :: &
       buffer_send_faces_scalar, buffer_recv_faces_scalar
  integer, dimension(MPI_STATUS_SIZE) :: msg_status
  integer, dimension(ninterface)  :: msg_requests
  integer :: rank, comm

  buffer_send_faces_scalar(:,:) = 0.d0
  buffer_recv_faces_scalar(:,:) = 0.d0

  do num_interface = 1, ninterface
     ipoin = 0
     do i = 1, nibool_interfaces(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val(ibool_interfaces(i,num_interface))
     enddo

     ! non-blocking send
     if(nibool_interfaces(num_interface).ne.0)then
      call MPI_ISEND( buffer_send_faces_scalar(1,num_interface), &
          max_ibool_interfaces_size,&
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          comm, msg_requests(num_interface), ier)
      endif
  enddo

  do num_interface = 1, ninterface

     ! starts a blocking receive
     if(nibool_interfaces(num_interface).ne.0)then
      call MPI_RECV ( buffer_recv_faces_scalar(1,num_interface), &
          max_ibool_interfaces_size, &
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          comm, msg_status(1), ier)
      endif
      call MPI_Barrier(comm,ierr)

     ipoin = 0
     do i = 1, nibool_interfaces(num_interface)
        ipoin = ipoin + 1
        array_val(ibool_interfaces(i,num_interface)) = &
             array_val(ibool_interfaces(i,num_interface))  &
             + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

  enddo


!     do num_interface = 1, ninterface
!       ipoin = 0
!       do i = 1, nibool_interfaces(num_interface)
!         ipoin = ipoin + 1
!         print  &
!           *,'mass',array_val(ibool_interfaces(i,num_interface)),num_interface,i,rank
!        if(rank.eq.) 
!       end do
!     end do


  ! synchronizes MPI processes
  call MPI_BARRIER(comm,ier)

  end subroutine assemble_MPI_mass


!-----------------------------------------------
! Assembling accel_elastic for elastic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimisations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_vector(array_val,npoin,ID, &
                ninterface, max_ibool_interfaces_size,&
                ibool_interfaces, nibool_interfaces, &
  !              buffer_send_faces_vector, &
  !              buffer_recv_faces_vector, &
                my_neighbours, rank, comm)

  use mpi
  use constants

  implicit none

  !include 'constants.h'

  integer, intent(in)  :: npoin
  double precision, dimension(1,npoin), intent(inout) :: array_val

  integer, dimension(:,:), intent(in) ::ID
  integer, intent(in)  :: ninterface
  integer, intent(in)  :: max_ibool_interfaces_size
  integer, dimension(max_ibool_interfaces_size,ninterface), intent(in)  :: ibool_interfaces
  integer, dimension(ninterface)   :: nibool_interfaces
  !integer, dimension(ninterface*2), intent(inout)  :: tab_requests_send_recv
  integer, dimension(ninterface*2)  :: tab_requests_send_recv
  !double precision, dimension(max_ibool_interfaces_size,ninterface), intent(inout) :: &
  !     buffer_send_faces_vector
  !double precision, dimension(max_ibool_interfaces_size,ninterface), intent(inout) :: &
  !     buffer_recv_faces_vector
  ! array to assemble
  double precision, dimension(3*max_ibool_interfaces_size,ninterface) :: &
       buffer_send_faces_vector
  double precision, dimension(3*max_ibool_interfaces_size,ninterface) :: &
       buffer_recv_faces_vector
  integer, dimension(ninterface), intent(in) :: my_neighbours

  integer  :: ipoin, num_interface, iinterface, ier, i, j,ierr
  integer, dimension(MPI_STATUS_SIZE)  :: stat
  integer :: rank, comm

  buffer_send_faces_vector(:,:)=0
  buffer_recv_faces_vector(:,:)=0


  do iinterface = 1, ninterface
    ipoin = 0
    do i = 1, nibool_interfaces(iinterface)
      do j=1,3
!        print *,'id',ID(:,:),rank
!        print *,'ibool i ii',ibool_interfaces(i,iinterface),i,iinterface,rank
!        print    *,'all',ipoin,iinterface,ID(j,ibool_interfaces(i,iinterface)),rank
        buffer_send_faces_vector(ipoin+1,iinterface) = &
             array_val(1,ID(j,ibool_interfaces(i,iinterface)))
        ipoin = ipoin + 1
      enddo
    enddo
    if(nibool_interfaces(iinterface).gt.0)then
!      print *, "max send:", maxval(buffer_send_faces_vector),rank, iinterface
    endif
  enddo
!  print *,"finish asse",rank
  call MPI_Barrier(comm,ierr)

  do iinterface = 1, ninterface

!    print *,'begin exe send',iinterface,nibool_interfaces(iinterface),rank
    if(nibool_interfaces(iinterface).gt.0)then
!      print *,'send exe',iinterface,rank,my_neighbours(iinterface)
      call MPI_ISEND( buffer_send_faces_vector(1,iinterface), &
        3*nibool_interfaces(iinterface), MPI_DOUBLE_PRECISION, &
        my_neighbours(iinterface), 12, comm, &
        tab_requests_send_recv(iinterface), ier)
        !print *, "bsfv:", buffer_send_faces_vector(1:3*nibool_interfaces(iinterface),iinterface), rank, my_neighbours(iinterface)
!        print *,'send exe1',iinterface,rank,my_neighbours(iinterface)
    endif
  end do
    
!  print *,"finish send 1",rank
  call MPI_Barrier(comm,ierr)
  do iinterface = 1, ninterface

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi(rank,comm,'MPI_ISEND unsuccessful in assemble_MPI_vector')
    endif
    
!    print *,"finish send",rank
    call MPI_Barrier(comm,ierr)


    if(nibool_interfaces(iinterface).gt.0)then
      call MPI_IRECV ( buffer_recv_faces_vector(1,iinterface), &
             3*nibool_interfaces(iinterface), MPI_DOUBLE_PRECISION, &
             my_neighbours(iinterface), 12, comm, &
             tab_requests_send_recv(ninterface+iinterface), ier)
       ! print *, "brfv:", buffer_recv_faces_vector, rank, my_neighbours(iinterface)
!      print *, "max rec:", maxval(buffer_recv_faces_vector(:,iinterface)),rank, iinterface
    endif

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi(rank,comm,'MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif

  enddo
!  print *,"finish rs",rank
  call MPI_Barrier(comm,ierr)

  do iinterface = 1, ninterface
    if(nibool_interfaces(iinterface).gt.0) then
!    print *,'wait sucess',iinterface,rank
       call MPI_Wait (tab_requests_send_recv(iinterface), stat, ier)
       call MPI_Wait (tab_requests_send_recv(iinterface+ninterface), stat, ier)
!    print *,'wait sucess1',iinterface,rank
     end if

  enddo
!  print *,"finish wait",rank
  call MPI_Barrier(comm,ierr)


  do iinterface = 1, ninterface
    ipoin = 0
    do i = 1, nibool_interfaces(iinterface)
      do j=1,3
        array_val(1,ID(j, ibool_interfaces(i,iinterface))) = &
            array_val(1,ID(j, ibool_interfaces(i,iinterface)))  &
            + buffer_recv_faces_vector(ipoin+1,iinterface)
        ipoin = ipoin + 1
      enddo
     enddo
  enddo
!  print *,"finish **",rank
  call MPI_Barrier(comm,ierr)

  end subroutine assemble_MPI_vector



!#endif
end module assemble_subs

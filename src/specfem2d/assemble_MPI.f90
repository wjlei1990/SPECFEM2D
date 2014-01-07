
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


#ifdef USE_MPI

!-----------------------------------------------
! Assembling the mass matrix.
!-----------------------------------------------
  subroutine assemble_MPI_scalar(array_val,npoin_val, &
                              ninterface, max_interface_size, &
                              max_ibool_interfaces_size, &
                              ibool_interfaces_elastic, &
                              nibool_interfaces_elastic, my_neighbours)

  use :: mpi

  implicit none

  include 'constants.h'

  integer, intent(in)  :: ninterface
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_ac,max_ibool_interfaces_size_el, &
    max_ibool_interfaces_size_po
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: &
    ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic,nibool_interfaces_elastic, &
    nibool_interfaces_poroelastic
  integer, dimension(ninterface), intent(in)  :: my_neighbours
  ! array to assemble
  ! acoustic
  integer :: npoin_val1
  real(kind=CUSTOM_REAL), dimension(npoin_val1), intent(inout) :: array_val1
  ! elastic
  integer :: npoin_val2
  real(kind=CUSTOM_REAL), dimension(npoin_val2), intent(inout) :: array_val2,array_val5
  ! poroelastic
  integer :: npoin_val3
  real(kind=CUSTOM_REAL), dimension(npoin_val3), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface
  integer  :: ier
  integer  :: i
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
  double precision, dimension(max_ibool_interfaces_size_ac+2*max_ibool_interfaces_size_el+&
       2*max_ibool_interfaces_size_po, ninterface)  :: &
       buffer_send_faces_scalar, &
       buffer_recv_faces_scalar
  integer, dimension(MPI_STATUS_SIZE) :: msg_status
  integer, dimension(ninterface)  :: msg_requests

  buffer_send_faces_scalar(:,:) = 0.d0
  buffer_recv_faces_scalar(:,:) = 0.d0

  do num_interface = 1, ninterface

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val1(ibool_interfaces_acoustic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val2(ibool_interfaces_elastic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val5(ibool_interfaces_elastic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val3(ibool_interfaces_poroelastic(i,num_interface))
     enddo
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val4(ibool_interfaces_poroelastic(i,num_interface))
     enddo

     ! non-blocking send
     call MPI_ISEND( buffer_send_faces_scalar(1,num_interface), &
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
          nibool_interfaces_acoustic(num_interface)+2*nibool_interfaces_elastic(num_interface)+&
          2*nibool_interfaces_poroelastic(num_interface), &
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, msg_requests(num_interface), ier)

  enddo

  do num_interface = 1, ninterface

     ! starts a blocking receive
     call MPI_RECV ( buffer_recv_faces_scalar(1,num_interface), &
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
          nibool_interfaces_acoustic(num_interface)+2*nibool_interfaces_elastic(num_interface)+&
          2*nibool_interfaces_poroelastic(num_interface), &
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, msg_status(1), ier)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        array_val1(ibool_interfaces_acoustic(i,num_interface)) = &
            array_val1(ibool_interfaces_acoustic(i,num_interface))  &
             + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        array_val2(ibool_interfaces_elastic(i,num_interface)) = &
            array_val2(ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
        array_val5(ibool_interfaces_elastic(i,num_interface)) = &
            array_val5(ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        array_val3(ibool_interfaces_poroelastic(i,num_interface)) = &
            array_val3(ibool_interfaces_poroelastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        array_val4(ibool_interfaces_poroelastic(i,num_interface)) = &
            array_val4(ibool_interfaces_poroelastic(i,num_interface)) &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

  enddo

  ! synchronizes MPI processes
  call MPI_BARRIER(mpi_comm_world,ier)

  end subroutine assemble_MPI_scalar


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
  subroutine assemble_MPI_vector_el(array_val2,npoin, &
                                   ninterface, &
                                   inum_interfaces_elastic, &
                                   max_interface_size, max_ibool_interfaces_size_el,&
                                   ibool_interfaces_elastic, nibool_interfaces_elastic, &
                                   tab_requests_send_recv_elastic, &
                                   buffer_send_faces_vector_el, &
                                   buffer_recv_faces_vector_el, &
                                   my_neighbours)

  use :: mpi

  implicit none

  include 'constants.h'
  include 'precision.h'

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_el
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  real(CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(inout)  :: &
       buffer_send_faces_vector_el
  real(CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(inout)  :: &
       buffer_recv_faces_vector_el
  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(3,npoin), intent(inout) :: array_val2
  integer, dimension(ninterface), intent(in) :: my_neighbours

  integer  :: ipoin, num_interface, iinterface, ier, i
  integer, dimension(MPI_STATUS_SIZE)  :: status_elastic


  do iinterface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        buffer_send_faces_vector_el(ipoin+1:ipoin+3,iinterface) = &
             array_val2(:,ibool_interfaces_elastic(i,num_interface))
        ipoin = ipoin + 3
     enddo

  enddo

  do iinterface = 1, ninterface_elastic

    num_interface = inum_interfaces_elastic(iinterface)

    call MPI_ISEND( buffer_send_faces_vector_el(1,iinterface), &
             3*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_ISEND unsuccessful in assemble_MPI_vector_el')
    endif

    call MPI_IRECV ( buffer_recv_faces_vector_el(1,iinterface), &
             3*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(ninterface_elastic+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_IRECV unsuccessful in assemble_MPI_vector_el')
    endif

  enddo

  do iinterface = 1, ninterface_elastic*2

    call MPI_Wait (tab_requests_send_recv_elastic(iinterface), status_elastic, ier)

  enddo

  do iinterface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        array_val2(:,ibool_interfaces_elastic(i,num_interface)) = &
            array_val2(:,ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_vector_el(ipoin+1:ipoin+3,iinterface)
        ipoin = ipoin + 3
     enddo

  enddo

  end subroutine assemble_MPI_vector_el



#endif


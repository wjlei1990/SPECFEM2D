module set_src_and_rec

  use mpi
  use constants

  implicit none
  !include "constants.h"

contains

  subroutine read_src_and_rec(ns, x_src, z_src,&
        nrec, x_rec, z_rec, rank, nproc, comm)

    integer :: ns, nrec
    double precision, allocatable :: x_src(:), z_src(:)
    double precision, allocatable :: x_rec(:), z_rec(:)
    integer :: rank, nproc, comm

    integer :: IIN=110
    integer :: i, j, ierr

    !=========
    !read src in master node
    if(rank.eq.0)then
      open(unit=IIN,file="DATA/source.txt")
      read(IIN,*) ns
    endif

    call MPI_Bcast(ns, 1, MPI_INTEGER, 0, comm, ierr)
    allocate(x_src(ns), z_src(ns))

    if(rank.eq.0)then
      do i=1,ns
        read(IIN, *) x_src(i), z_src(i)
      enddo
    endif

    call MPI_Bcast(x_src, ns, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(z_src, ns, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    !===========
    !read rec in master node
    if(rank.eq.0)then
      open(unit=IIN,file="DATA/receiver.txt")
      read(IIN,*) nrec
    endif

    call MPI_Bcast(nrec, 1, MPI_INTEGER, 0, comm, ierr)
    allocate(x_rec(nrec), z_rec(nrec))

    if(rank.eq.0)then
      do i=1,nrec
        read(IIN, *) x_rec(i), z_rec(i)
      enddo
    endif

    call MPI_Bcast(x_rec, nrec, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(z_rec, nrec, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    !Bcast

  end subroutine read_src_and_rec


  subroutine locate_src_and_rec(x, z, ibool, ns, sglob, x_src, z_src, &
        nrec, rglob, x_rec, z_rec, rank, nproc, comm)

    double precision :: x(:), z(:)
    integer :: ibool(:,:,:)
    integer :: ns, nrec
    integer :: sglob(:), rglob(:)
    double precision :: x_src(:), z_src(:), x_rec(:), z_rec(:)
    integer :: rank, nproc, comm
    integer :: min_proc
    
    integer :: i, ierr
    double precision :: short_dist
    character(len=100) :: string_proc, fn

    do i=1, ns
      call find_close_grid_point(x,z,ibool, x_src(i),z_src(i), &
            short_dist,sglob(i),min_proc,rank,nproc,comm)
    enddo 
!    print *,"rank,sglob:",rank, sglob(1)

    write(string_proc,'(i5.5)')rank
    fn="OUTPUT_FILES/source_proc"//trim(string_proc)
!    print *,"fn", fn, rank
    !call MPI_Barrier(comm, ierr)
    !stop
    open(11,file=fn)
    do i=1,ns
      if(sglob(i).gt.0)then
        write(11,*) x_src(i), z_src(i), sglob(i), rank, x(sglob(i)), z(sglob(i))
      endif
    enddo
    close(11)

    do i=1, nrec
      call find_close_grid_point(x,z,ibool,x_rec(i),z_rec(i), &
            short_dist,rglob(i),min_proc,rank,nproc,comm)
    enddo
    fn="OUTPUT_FILES/receiver_proc"//trim(string_proc)
    open(11,file=fn)
    do i=1,nrec
      if(rglob(i).gt.0)then
        write(11,*) x_rec(i), z_rec(i), rglob(i), rank, x(rglob(i)), z(rglob(i))
      endif
    enddo
    close(11)
!    print *,"rank,rglob:",rank, rglob(1)

  end subroutine locate_src_and_rec

  subroutine find_close_grid_point(x,z,ibool,x_target,z_target,short_dist,&
                target_glob, min_proc, rank, nproc, comm)

    use wave2d_variables, only:NEX, NEZ

    double precision,intent(in) :: x(:), z(:)
    integer,intent(in) :: ibool(:,:,:)
    double precision, intent(in) :: x_target, z_target
    double precision :: short_dist
    integer :: target_glob
    integer :: rank, nproc, comm

    integer :: i, j,ii,jj, ix, iz, ispec, iglob
    integer :: min_proc
    double precision :: dist_temp
    double precision :: dist_temp_array(nproc)

    integer :: ierr

    !init
    short_dist=10000.0
    ispec=0
    !loop all elements to find the closest point
    do j=1, NEZ
      do i=1, NEX
        ispec=ispec+1
        do jj=1, NGLLZ
          do ii=1, NGLLX
            iglob=ibool(ii,jj,ispec)
            dist_temp=distance(x(iglob), z(iglob), x_target, z_target)
            if(dist_temp.lt.short_dist)then
              short_dist=dist_temp
              target_glob=ibool(ii,jj,ispec)
            endif
          enddo
        enddo
      enddo
    enddo

    !print *, "rank, short_dist:", rank, short_dist
    call MPI_Barrier(comm,ierr)

    call MPI_Gather(short_dist,1, MPI_DOUBLE_PRECISION, dist_temp_array, 1, &
      MPI_DOUBLE_PRECISION, 0, comm, ierr)

    if(rank.eq.0) then
      !print *,"gather dist:",dist_temp_array(:)
      dist_temp=dist_temp_array(1)
      min_proc=0
      do i=1, nproc
        if(dist_temp_array(i).lt.dist_temp)then
          dist_temp=dist_temp_array(i)
          min_proc=i-1
        endif
      enddo
    endif
    
    call MPI_Bcast(min_proc, 1, MPI_INTEGER, 0, comm, ierr);

    if(rank.ne.min_proc)then
      target_glob=0
    endif

  end subroutine find_close_grid_point

  double precision function distance(x1,z1,x2,z2)

    double precision :: x1, z1, x2, z2

    distance=sqrt((x2-x1)**2+(z2-z1)**2)

  end function distance

end module set_src_and_rec

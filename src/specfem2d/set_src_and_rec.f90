module set_src_and_rec


contains

  subroutine read_src_and_rec(ns, x_src, z_src,&
        nrec, x_rec, z_rec, rank, nproc, comm)

    integer :: ns, nrec
    double precision :: x_src(:), z_src(:)
    double precision :: x_rec(:), z_rec(:)

    integer :: IIN=110
    integer :: i,j

    !=========
    !read src in master node
    if(rank.eq.0)then
      open(unit=IIN,file="source.txt")
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
      open(unit=IIN,file="receiver.txt")
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


  subroutine locate_src_and_rec(ns, sglob, x_src, z_src, &
        nrec, rglob, x_rec, z_rec, rank, nproc, comm)


  end subroutine locate_src_and_rec



end module set_src_and_rec

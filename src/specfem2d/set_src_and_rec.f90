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


  subroutine locate_src_and_rec(x, z, ibool, ns, sglob, x_src, z_src, &
        nrec, rglob, x_rec, z_rec, rank, nproc, comm)

    double precision :: x(:), z(:)
    integer :: ibool(:,:,:)
    integer :: ns, nrec
    integer :: sglob(:), rglob(:)
    double precision :: x_src(:), z_src(:), x_rec(:), z_rec(:)
    integer :: rank, nproc, comm
    
    integer :: i
    double precision :: short_dist

    do i=1, ns
      call find_close_grid_point(x,z,ibool, x_src(i),z_src(i), &
            short_dist,sglob(i),min_proc,rank,nproc,comm)
    enddo 

    do i=1, nrec
      call find_close_grid_point(x,z,ibool,x_rec(i),z_rec(i), &
            short_dist,rglob(i),min_proc,rank,nproc,comm)
    enddo

  end subroutine locate_src_and_rec

  subroutine find_close_grid_point(x,z,ibool,x_target,z_target,short_dist,&
                target_glob, min_proc, rank, nproc, comm)

    double precision,intent(in) :: x(:), z(:)
    integer,intent(in) :: ibool(:,:,:)
    double precision, intent(in) :: x_target, z_target, short_dist
    integer :: target_glob
    integer :: rank, nproc, comm

    integer :: i, j, ix, iz, ispec, iglob
    integer :: min_proc
    double precision :: dist_temp
    double precision :: dist_temp_array(nproc)

    !init
    short_dist=10000.0
    ispec=0
    !loop all elements to find the closest point
    do j=1, NEZ
      do i=1, NEX
        ispec=ispec+1
        do j=1, NGLLZ
          do i=1, NGLLX
            iglob=ibool(i,j,ispec)
            dist_temp=distance(x(iglob), z(iglob), x_rec, z_rec)
            if(dist_temp.lt.short_dist)then
              short_dist=dits_temp
              target_glob=ibool(i,j,ispec)
            endif
          enddo
        enddo
      enddo
    enddo

    call MPI_Gather(short_dist,1, MPI_DOUBLE_PRECISION, dist_temp_array, &
      MPI_DOUBLE_PRECISION, 0, comm)

    if(rank.eq.0) then
      dist_temp=dist_temp_array(1)
      min_proc=1
      do i=2, nproc
        if(dist_temp_array(i).lt.dist_temp)then
          dist_temp=dist_temp_array(i)
          min_proc=i
        endif
      enddo
    endif
    
    call MPI_Bcast(min_proc, 1, MPI_INTEGER, 0, comm);

    if(rank.ne.min_proc)then
      target_glob=0
    endif

  end subroutine find_close_grid_point

  double precision function distance(x1,z1,x2,z2)

    double precision :: x1, z1, x2, z2

    distance=sqrt((x2-x1)**2+(z2-z1)**2)

  end function distance

end module set_src_and_rec

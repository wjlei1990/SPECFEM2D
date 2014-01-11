module wave2d_io_subs

  !use wave2d_constants
  !use wave2d_variables

  implicit none

contains
  
  subroutine read_parfile(LOCAL_PATH, OUTPUT_PATH, &
         rank, nproc, comm)

    use wave2d_variables
    use mpi

    character(len=*) :: LOCAL_PATH, OUTPUT_PATH

    integer :: rank, nproc, comm

    integer :: idummy=5, i, ierr
    character(len=30) :: dstring
    integer :: IIN=110

    INTEGER :: NEX_XI, NEX_ZI, NPROC_XI, NPROC_ZI
    double precision :: MODEL_X1, MODEL_X2,MODEL_Z1,MODEL_Z2
    double precision ::DENSITY, INCOMPRESSIBILITY, RIGIDITY
    character(len=10) :: MODEL_TYPE


    open(unit=IIN,file="DATA/Par_file")

    if(rank.eq.0)then
      do i=1,idummy
        read(IIN,*)
      enddo
      read(IIN,5) dstring, SIMUL_TYPE 
      read(IIN,*)
      read(IIN,*)
      read(IIN,5) dstring, NEX_XI
      read(IIN,5) dstring, NEX_ZI
      read(IIN,*)
      read(IIN,*)
      read(IIN,5) dstring, NPROC_XI
      read(IIN,5) dstring, NPROC_ZI
      read(IIN,*)
      read(IIN,*)
      read(IIN,6) dstring, LOCAL_PATH
      read(IIN,*)
      read(IIN,*)
      read(IIN,6) dstring, OUTPUT_PATH
      read(IIN,*)
      read(IIN,*)
      read(IIN,7) dstring, RECORD_LENGTH_IN_SECONDS
      read(IIN,*)
      read(IIN,*)
      read(IIN,7) dstring, DT
      read(IIN,*)
      read(IIN,*)
      read(IIN,7) dstring, MODEL_X1
      read(IIN,7) dstring, MODEL_X2
      read(IIN,7) dstring, MODEL_Z1
      read(IIN,7) dstring, MODEL_Z2
      read(IIN,*)
      read(IIN,*)
      read(IIN,6) dstring, MODEL_TYPE
      read(IIN,7) dstring, DENSITY
      read(IIN,7) dstring, INCOMPRESSIBILITY
      read(IIN,7) dstring, RIGIDITY
      read(IIN,*)
      read(IIN,*)
      read(IIN,8) dstring, DEBUG
      !check
      if(nproc.ne.(NPROC_XI*NPROC_ZI))then
        call exit_mpi(rank, comm, "error in number of proc")
      endif
    endif
    call MPI_Bcast(SIMUL_TYPE,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(NEX_XI,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(NEX_ZI,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(NPROC_XI,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(NPROC_ZI,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(LOCAL_PATH,300,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(OUTPUT_PATH,300,MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(RECORD_LENGTH_IN_SECONDS,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(DT,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(DEBUG,1,MPI_LOGICAL,0,comm,ierr)

    NEX=NEX_XI/NPROC_XI
    NEZ=NEX_ZI/NPROC_ZI
    !NELE=max(NEX, NEZ)
    !NSPEC=NEX*NEZ 
    NGLLSQUARE=NGLLX*NGLLZ
    NGLL=max(NGLLX, NGLLZ)
    !NGLOB=((NGLLX-1)*NEX +1)*((NGLLZ-1)*NEZ+1)

    NSTEP=RECORD_LENGTH_IN_SECONDS/DT

5 format(A,I5)
6 format(A,A)
7 format(A,E15.5)
8 format(A,L5)

  end subroutine read_parfile

  !subroutine write_snapshot(disp, filename)

   ! double precision,intent(in) :: disp(3,NGLOB)
   ! character(len=200),intent(in) :: filename
!
!    integer :: icomp, iglob, ios
!
!    open(unit = 11, file = trim(filename), status = 'unknown',iostat=ios)
!    if (ios /= 0) stop 'Error writing snapshot to disk'
!    do iglob = 1, NGLOB
!       write(11,'(5e12.3)') x(iglob)/LENGTH, z(iglob)/LENGTH, &
!                  sngl(disp(1,iglob)),sngl(disp(2,iglob)), sngl(disp(3,iglob))
!    enddo
!    close(11)
!
!  end subroutine write_snapshot

  !-----------------------------------------------------   

  subroutine write_seismogram(seis, rglob, nrec, DT, NSTEP, OUTPUT_PATH)

    integer, intent(in) :: nrec, NSTEP
    integer, intent(in) :: rglob(:)
    double precision, intent(in) :: DT
    double precision, intent(in) ::  seis(NSTEP,3,nrec)
    character(len=*),intent(in) :: OUTPUT_PATH 
    integer :: ios

    character(len=200) :: file_name
    integer :: irec, icomp, itime

    call system('mkdir -p '//trim(OUTPUT_PATH)//'')

    do irec = 1, nrec
      if(rglob(irec).gt.0)then
        do icomp = 1, 3
          write(file_name,'(a,a,i2.2,a,i1.1)') trim(OUTPUT_PATH), &
            "/seismogram_", irec, '_', icomp
          open(unit = 12, file = file_name, status = 'unknown', iostat=ios)
          if (ios /= 0) stop 'Error opening seismogram to write'
          do itime = 1, NSTEP
            write(12, *) DT * itime, seis(itime, icomp, irec)
          enddo
          close(12)
        enddo
      endif
    enddo

  end subroutine write_seismogram


  subroutine write_stf(samp,nstep,ns,OUTPUT_PATH)
    integer ::nstep,ns
    double precision,dimension(nstep,3,ns)::samp

    !local 
    integer ::istep,is,ios
    character(len=*),intent(in) :: OUTPUT_PATH
    character(len=200) ::outfile

    do is=1,ns
      write(outfile,'(a,a,i2.2,a,i1.1)') trim(OUTPUT_PATH), &
                  "/source_", is
 
      open(unit=12,file = outfile, status= 'unknown',iostat=ios)
      do istep=1,nstep
        write(12,*) samp(istep,1,is),samp(istep,2,is),samp(istep,3,is) 
      end do
      close(12)

    end do
  end subroutine write_stf
    

  !------------------------------------------------

!end module wave2d_sub
end module wave2d_io_subs

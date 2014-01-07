module wave2d_sub

  use wave2d_constants
  use wave2d_variables

  implicit none

contains

  ! -------------------------------------------

  subroutine set_glob(nrec, x_rec, z_rec, rglob)

    integer, intent(in) :: nrec
    double precision, intent(in) :: x_rec(nrec), z_rec(nrec)
    integer, intent(out) :: rglob(nrec)
    
    double precision :: d_min_rec, d
    integer :: irec, ispec, i, j, iglob

    do irec = 1, nrec

       d_min_rec = sqrt(LENGTH**2+HEIGHT**2)

       do ispec = 1,NSPEC
          do j = 1,NGLLZ
             do i = 1,NGLLX
                iglob = ibool(i,j,ispec)
                d = sqrt((x_rec(irec)-x(iglob))**2+(z_rec(irec)-z(iglob))**2)
                if(d < d_min_rec) then
                  d_min_rec  = d
                  rglob(irec) = ibool(i,j,ispec)
                endif
             enddo
          enddo
       enddo

    enddo


  end subroutine set_glob


  !---------------------------------------

  subroutine write_snapshot(disp, filename)

    double precision,intent(in) :: disp(3,NGLOB)
    character(len=200),intent(in) :: filename

    integer :: icomp, iglob, ios

    open(unit = 11, file = trim(filename), status = 'unknown',iostat=ios)
    if (ios /= 0) stop 'Error writing snapshot to disk'
    do iglob = 1, NGLOB
       write(11,'(5e12.3)') x(iglob)/LENGTH, z(iglob)/LENGTH, &
                  sngl(disp(1,iglob)),sngl(disp(2,iglob)), sngl(disp(3,iglob))
    enddo
    close(11)

  end subroutine write_snapshot

  !-----------------------------------------------------   

  subroutine write_seismogram(seis, nrec, seis_name)

    integer, intent(in) :: nrec
    double precision, intent(in) ::  seis(NSTEP,3,nrec)
    character(len=*),intent(in) :: seis_name

    character(len=200) :: file_Name
    integer :: irec, icomp, itime

    do irec = 1, nrec
       do icomp = 1, 3
          write(file_name,'(a,a,i2.2,a,i1.1)') trim(seis_name), '_', irec, '_', icomp
          open(unit = 12, file = file_name, status = 'unknown', iostat=ios)
          if (ios /= 0) stop 'Error opening seismogram to write'
          do itime = 1, NSTEP
             write(12, *) DT * itime, seis(itime, icomp, irec)
          enddo
          close(12)
       enddo
    enddo

  end subroutine write_seismogram

  !------------------------------------------------

end module wave2d_sub

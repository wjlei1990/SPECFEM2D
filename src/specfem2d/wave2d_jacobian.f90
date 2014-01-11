module wave2d_jacobian

  use wave2d_variables
  use wave2d_define_der_matrices

  implicit none

contains

  subroutine cal_jacobian()

    integer ispec,ib,i,j,k,l,iglob,iglob1,itime,ix,iz,ii,jj
    integer N_EQ,n_fix
    integer,dimension(5) :: topo_array
    integer :: ver_loc
    double precision :: ver_fac


    !topo_array=(/3,5,7,10/)

    ! set up grid and spectral elements
    call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

    ispec = 0
    iglob = 0
    !nspecb(:) = 0
    !N_EQ=0
    n_fix=0

    dxdxi(:,:,:)=0.
    dzdxi(:,:,:)=0.
    dxdgamma(:,:,:)=0.
    dzdgamma(:,:,:)=0.
    jacobian(:,:,:)=0.

    !------------------------------------------------------------------------
    ! loop over all elements
    do iz = 1,NEZ
      do ix = 1,NEX
        ispec = ispec+1
      !------------------------------------------------------------------------

        ! loop over GLL points to calculate jacobian, and set up numbering
        !
        ! jacobian = | dx/dxi dx/dgamma | = (z2-z1)*(x2-x1)/4  as dx/dgamma=dz/dxi = 0
        !            | dz/dxi dz/dgamma | 
        !
        do j=1,NGLLZ
          do i=1, NGLLX
            ! jacobian, integration weight
            do k=1,NGLLX
              dxdxi(i,j,ispec) = dxdxi(i,j,ispec)+hprime_xx(k,i)*x(ibool(k,j,ispec))
              dzdxi(i,j,ispec) = dzdxi(i,j,ispec)+hprime_xx(k,i)*z(ibool(k,j,ispec))
            enddo
            !print *, "dxdxi,dzdxi",dxdxi(i,j,ispec),dzdxi(i,j,ispec)
            do k=1,NGLLZ
               dxdgamma(i,j,ispec) = dxdgamma(i,j,ispec)+hprime_zz(k,j)*x(ibool(i,k,ispec))
               dzdgamma(i,j,ispec) = dzdgamma(i,j,ispec)+hprime_zz(k,j)*z(ibool(i,k,ispec))
            enddo
            !print *,"dxdgamma,dzdgamma",dxdgamma(i,j,ispec),dzdgamma(i,j,ispec)
            !jacobian(i,j,ispec) = (z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4. 
            jacobian(i,j,ispec)=dxdxi(i,j,ispec)*dzdgamma(i,j,ispec)-dzdxi(i,j,ispec)*dxdgamma(i,j,ispec)
            !inverse the jacobian for each GLL point in the element
            dxidx(i,j,ispec)=dzdgamma(i,j,ispec)/jacobian(i,j,ispec)
            dgammadx(i,j,ispec)=-dzdxi(i,j,ispec)/jacobian(i,j,ispec)
            dxidz(i,j,ispec)=-dxdgamma(i,j,ispec)/jacobian(i,j,ispec)
            dgammadz(i,j,ispec)=dxdxi(i,j,ispec)/jacobian(i,j,ispec)
            !print *, "11,12,21,22"
            !print *, dxidx(i,j,ispec),dgammadx(i,j,ispec)
            !print *, dxidz(i,j,ispec),dgammadz(i,j,ispec)
            !end loop over GLL points
            !if(jacobian(i,j,ispec).ne.((z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4.) ) then
            ! print *, "Yes, it is right!"
            !endif
          end do
        end do
      enddo
    enddo

    jacobianb(:,:,:)=0.0
    do i=1,4
      do j=1,nspecb(i)
        if(i.eq.1.or.i.eq.2) then
          do k=1,NGLLZ
            jacobianb(i,k,j)=abs(z(ibool(1,1,ibelm(j,i)))-z(ibool(1,NGLLZ,ibelm(j,i))))/2.0
          end do
        else if(i.eq.3) then
          do k=1,NGLLX
            jacobianb(3,k,j)=abs(x(ibool(1,1,ibelm(j,i)))-x(ibool(NGLLX,1,ibelm(j,i))))/2.0
          end do
        else if(i.eq.4) then
          do k=1,NGLLX
            do l=1,NGLLX
               jacobianb(4,k,j) = jacobianb(4,k,j) &
                                  + hprime_xx(l,k)*x(ibool(l,NGLLZ,ibelm(j,i)))
            end do
          end do
        end if
      end do
    end do

    ! estimate the time step
    !dh = HEIGHT/dble((NGLLZ-1)*NEZ)
    !if(dh > LENGTH/dble((NGLLX-1)*NEX)) dh = LENGTH/dble((NGLLX-1)*NEX)
    !c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
    !time_step = 0.2*dh/c
    !print *
    !print *,'time step estimate from courant = 0.2: ',sngl(time_step),' seconds'
    print *,'  actual time step: ',sngl(DT),' seconds'
    !print *,'Jacobi',maxval(jacobian),maxval(jacobianb),rank
   
  end subroutine cal_jacobian

end module wave2d_jacobian

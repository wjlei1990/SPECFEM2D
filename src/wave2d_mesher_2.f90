module wave2d_mesher

  use wave2d_constants
  use wave2d_variables
  use wave2d_define_der_matrices

  implicit none

contains

  subroutine mesher

    integer ispec,ib,i,j,k,iglob,iglob1,itime,ix,iz,ii,jj

    ! set up grid and spectral elements
    call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

    ispec = 0
    iglob = 0
    nspecb(:) = 0

    ! loop over all elements
    do iz = 1,NEZ
       do ix = 1,NEX
          ispec = ispec+1

          ! evenly spaced achors between 0 and 1
          x1(ispec) = LENGTH*dble(ix-1)/dble(NEX)
          x2(ispec) = LENGTH*dble(ix)/dble(NEX)
          z1(ispec) = HEIGHT*dble(iz-1)/dble(NEZ)
          z2(ispec) = HEIGHT*dble(iz)/dble(NEZ)

          ! loop over GLL points to calculate jacobian, and set up numbering
          !
          ! jacobian = | dx/dxi dx/dgamma | = (z2-z1)*(x2-x1)/4  as dx/dgamma=dz/dxi = 0
          !            | dz/dxi dz/dgamma | 
          !
          do j = 1,NGLLZ
             do i = 1,NGLLX
                ! set up local to global numbering
                if ( (i.eq.1).and.(ix.gt.1) ) then
                   ibool(i,j,ispec) = ibool(NGLLX,j,ispec-1)
                else if ( (j.eq.1).and.(iz.gt.1) ) then
                   ibool(i,j,ispec) = ibool(i,NGLLZ,ispec-NEX)
                else
                   iglob = iglob + 1
                   ibool(i,j,ispec) = iglob         
                endif
                ! get the global grid points
                iglob1 = ibool(i,j,ispec)
                x(iglob1) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
                z(iglob1) = 0.5*(1.-zigll(j))*z1(ispec)+0.5*(1.+zigll(j))*z2(ispec)
             enddo
          enddo

          do j=1,NGLLZ
            do i=1, NGLLX
               ! jacobian, integration weight
               do k=1,NGLLX
                 dxdxi(i,j,ispec) = dxdxi(i,j,ispec)+hprime_xx(k,j)*x(ibool(k,j,ispec))
                 dzdxi(i,j,ispec) = dzdxi(i,j,ispec)+hprime_xx(k,j)*z(ibool(k,j,ispec))
               enddo
               do k=1,NGLLZ
                 dxdgamma(i,j,ispec) = dxdgamma(i,j,ispec)+hprime_zz(i,k)*x(ibool(i,k,ispec))
                 dzdgamma(i,j,ispec) = dzdgamma(i,j,ispec)+hprime_zz(i,k)*z(ibool(i,k,ispec))
               enddo
                 !jacobian(i,j,ispec) = (z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4. 
                 jacobian(i,j,ispec) = dxdxi(i,j,ispec)*dzdgamma(i,j)-dzdxi(i,j,ispec)*dxdgamma(i,j,ispec)
                 ! end loop over GLL points
                 if(jacobian(i,j,ispec)==(z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4.)
                  print *, "Yes, it is right!"
                endif
             end do
          end do

          !print *,ibool(1,1,1), ibool(1,2,1),ibool(2,1,1)

          ! if boundary element
          if (ix.eq.1) then         
             nspecb(1) = nspecb(1) + 1
             ibelm(1,nspecb(1)) = ispec
             do j = 1,NGLLZ
                jacobianb(1,j,nspecb(1))= (z2(ispec)-z1(ispec))/2.
             end do
          endif
          if (ix.eq.NEX) then
             nspecb(2) = nspecb(2) + 1
             ibelm(2,nspecb(2)) = ispec
             do j = 1,NGLLZ
                jacobianb(2,j,nspecb(2))= (z2(ispec)-z1(ispec))/2.
             end do
          endif
          if (iz.eq.1) then
             nspecb(3) = nspecb(3) + 1
             ibelm(3,nspecb(3)) = ispec
             do i = 1,NGLLX
                jacobianb(3,i,nspecb(3))= (x2(ispec)-x1(ispec))/2.
             end do
          endif
          if (iz.eq.NEZ) then
             nspecb(4) = nspecb(4) + 1
             ibelm(4,nspecb(4)) = ispec
             do i = 1,NGLLX
                jacobianb(4,i,nspecb(4))= (x2(ispec)-x1(ispec))/2.
             end do
          endif
          ! end loop over elements
       end do
    end do


    ! estimate the time step
    dh = HEIGHT/dble((NGLLZ-1)*NEZ)
    if(dh > LENGTH/dble((NGLLX-1)*NEX)) dh = LENGTH/dble((NGLLX-1)*NEX)
    c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
    time_step = 0.2*dh/c
    print *
    print *,'time step estimate from courant = 0.2: ',sngl(time_step),' seconds'
    print *,'  actual time step: ',sngl(DT),' seconds'
   
  end subroutine mesher

  !-------------------------------------------------------


  subroutine set_model_property()

    integer :: ispec, i, j, iglob

    ! properties
    do ispec = 1,NSPEC
       !  get center of element
       iglob = ibool(NGLLX/2,NGLLZ/2,ispec)
       do j = 1,NGLLZ
          do i = 1,NGLLX
             if(z(iglob) >= HEIGHT*0.5) then
                ! crust
                rho(i,j,ispec) = DENSITY
                kappa(i,j,ispec) = INCOMPRESSIBILITY
                mu(i,j,ispec) = RIGIDITY
             else
                ! mantle
                rho(i,j,ispec) = 3380.
                kappa(i,j,ispec) = 1.3d+11
                mu(i,j,ispec) = 6.8d+10
             endif
          end do
       end do
    end do

  end subroutine set_model_property

!---------------------------------------------

end module wave2d_mesher

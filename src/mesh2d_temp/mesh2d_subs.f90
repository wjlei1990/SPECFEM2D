module wave2d_mesh

  use wave2d_constants
  !use wave2d_variables
  use wave2d_define_der_matrices

  implicit none

contains

  subroutine mesh_one_proc(mesh_anchor, NEX, NEZ)

    double precision :: mesh_anchor(2,4)
    integer :: NEX, NEZ

    double precision :: LENGTH, HEIGHT

    integer ispec,ib,i,j,k,iglob,iglob1,itime,ix,iz,ii,jj
    integer N_EQ,n_fix

    ! set up grid and spectral elements
    call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

    ispec = 0
    iglob = 0
    nspecb(:) = 0
    N_EQ=0
    n_fix=0

    dxdxi(:,:,:)=0.
    dzdxi(:,:,:)=0.
    dxdgamma(:,:,:)=0.
    dzdgamma(:,:,:)=0.

    !-----------------------------------------------------------
    LENGTH=abs(mesh_anchor(1,2)-mesh_anchor(1,1))
    HEIGHT=abs(mesh_anchor(2,4)-mesh_anchor(2,1))
    ! loop over all elements
    do iz = 1,NEZ
       do ix = 1,NEX
          ispec = ispec+1

          ! evenly spaced achors between 0 and 1
          x1(ispec) = LENGTH*dble(ix-1)/dble(NEX)+mesh_anchor(1,1)
          x2(ispec) = LENGTH*dble(ix)/dble(NEX)+mesh_anchor(1,1)
          z1(ispec) = HEIGHT*dble(iz-1)/dble(NEZ)+mesh_anchor(2,1)
          z2(ispec) = HEIGHT*dble(iz)/dble(NEZ)+mesh_anchor(2,1)

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

          !print *,ibool(1,1,1), ibool(1,2,1),ibool(2,1,1)

          ! if boundary element
          if (ix.eq.1) then         
            if(abs(x1(ispec)-MODEL_X1).lt.TOL)then
              nspecb(1) = nspecb(1) + 1
              ibelm(1,nspecb(1)) = ispec
              do j = 1,NGLLZ
                jacobianb(1,j,nspecb(1))= (z2(ispec)-z1(ispec))/2.
              end do
            endif
          endif
          if (ix.eq.NEX) then
            if(abs(x2(ispec)-MODEL_X2).lt.TOL)then
              nspecb(2) = nspecb(2) + 1
              ibelm(2,nspecb(2)) = ispec
              do j = 1,NGLLZ
                jacobianb(2,j,nspecb(2))= (z2(ispec)-z1(ispec))/2.
              enddo
            endif
          endif
          if (iz.eq.1) then
            if(abs(z1(ispec)-MODEL_Z1).lt.TOL)then
              nspecb(3) = nspecb(3) + 1
              ibelm(3,nspecb(3)) = ispec
              do i = 1,NGLLX
                jacobianb(3,i,nspecb(3))= (x2(ispec)-x1(ispec))/2.
              end do
            endif
          endif
          if (iz.eq.NEZ) then
            if(abs(z2(ispec)-MODEL_Z2).lt.TOL)then
              nspecb(4) = nspecb(4) + 1
              ibelm(4,nspecb(4)) = ispec
              do i = 1,NGLLX
                do k=1,NGLLX
                  jacobianb(4,i,nspecb(4)) = jacobianb(4,i,nspecb(4))+hprime_xx(k,i)*x(ibool(k,NGLLZ,ispec))
                enddo
               ! jacobianb(4,i,nspecb(4))= (x2(ispec)-x1(ispec))/2.
              end do
             !print *,"jacobianb(4,:,:)"
             !print *, jacobianb(4,:,:)
           endif
          endif
          ! end loop over elements
       end do
    end do
    
    call set_ID_array()
    !call set_IEN_array()
    !call set_LM_array()

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
  subroutine set_ID_array()

    integer ix,iz,i,j,ispec
    integer n_fix

    N_EQ=0
    n_fix=0
    !get the ID value
    if(SIMUL_TYPE==1)then
      !loop over element
      do iz = 1,NEZ
        do ix = 1,NEX
          ispec=(iz-1)*NEX+ix
          !loop over node
          do j=1,NGLLZ
            do i=1,NGLLX
              !do not do overlap things
              if((i==1).and.(ix.gt.1)) then
              elseif((j==1).and.(iz.gt.1))then
              else
                N_EQ=N_EQ+1
                ID(1,ibool(i,j,ispec))=N_EQ
              endif
            enddo
          enddo
        enddo
      enddo
      !loop over element
      do iz = 1,NEZ
        do ix = 1,NEX
          ispec=(iz-1)*NEX+ix
          !loop over node
          do j=1,NGLLZ
            do i=1,NGLLX
              !do not do overlap things
              if((i==1).and.(ix.gt.1)) then
              elseif((j==1).and.(iz.gt.1))then
              else
                if(iz==1.and.j==1)then
                  n_fix=n_fix+1
                  ID(2,ibool(i,j,ispec))=0
                else
                  N_EQ=N_EQ+1
                  ID(2,ibool(i,j,ispec))=N_EQ
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    !----------------------
    !wave simulation
    if(SIMUL_TYPE==2)then
      !loop over element: x
      do iz = 1,NEZ
        do ix = 1,NEX
          ispec=(iz-1)*NEX+ix
          !loop over node
          do j=1,NGLLZ
            do i=1,NGLLX
              !do not do overlap things
              if((i==1).and.(ix.gt.1)) then
              elseif((j==1).and.(iz.gt.1))then
              else
                N_EQ=N_EQ+1
                ID(1,ibool(i,j,ispec))=N_EQ
              endif
            enddo
          enddo
        enddo
      enddo
      !loop over element: y
      do iz = 1,NEZ
        do ix = 1,NEX
          ispec=(iz-1)*NEX+ix
          !loop over node
          do j=1,NGLLZ
            do i=1,NGLLX
              !do not do overlap things
              if((i==1).and.(ix.gt.1)) then
              elseif((j==1).and.(iz.gt.1))then
              else
                N_EQ=N_EQ+1
                ID(2,ibool(i,j,ispec))=N_EQ
              endif
            enddo
          enddo
        enddo
      enddo
      !loop over element: z
      do iz = 1,NEZ
        do ix = 1,NEX
          ispec=(iz-1)*NEX+ix
          !loop over node
          do j=1,NGLLZ
            do i=1,NGLLX
              !do not do overlap things
              if((i==1).and.(ix.gt.1)) then
              elseif((j==1).and.(iz.gt.1))then
              else
                N_EQ=N_EQ+1
                ID(3,ibool(i,j,ispec))=N_EQ
              endif
            enddo
          enddo
        enddo
      enddo
    endif


    !-----------------------------
    if(DEBUG)then
      print *,"-------------------"
      print *,"ID(1,:)"
      print *,ID(1,:)
      print *,"ID(2,:)"
      print *,ID(2,:)
      print *,"ID(3,:)"
      print *,ID(3,:)
      print *,"--------------------"
    endif

  end subroutine set_ID_array
!--------------------------------------------------------
  subroutine set_IEN_array()
    integer :: i,j,ix,iz,ispec
    
    ispec=0
    do iz = 1,NEZ
      do ix = 1,NEX
        ispec=(iz-1)*NEX+ix
        do j=1,NGLLZ
          do i=1,NGLLX
            IEN(i,j,ispec)=ibool(i,j,ispec)
          enddo
        enddo
      enddo
    enddo


    if(DEBUG)then
      do ispec=1,NSPEC
        print *, "Element",ispec
        print *,"IEN"
        print *,"--------------------"
        do j=1,NGLLZ
          print *,IEN(:,j,ispec)
        enddo
      enddo
    endif
  end subroutine set_IEN_array

!--------------------------------------------------------
  subroutine set_LM_array()
    integer :: i,j,ix,iz,ispec

    ispec=0
    do iz = 1,NEZ
      do ix = 1,NEX
        ispec=(iz-1)*NEX+ix
        do j=1,NGLLZ
          do i=1,NGLLX
            LMX(i,j,ispec)=ID(1,IEN(i,j,ispec))
            LMZ(i,j,ispec)=ID(2,IEN(i,j,ispec))
          enddo
        enddo
      enddo
    enddo

    if(DEBUG)then
      do ispec=1,NSPEC
        print *, "ELEment",ispec
        print *,"LMX"
        print *,"--------------------"
        do j=1,NGLLZ
          print *,LMX(:,j,ispec)
        enddo
        print *,"LMZ"
        print *,"--------------------"
        do j=1,NGLLZ
          print *,LMZ(:,j,ispec)
        enddo
      enddo
    endif

    !stop
  end subroutine set_LM_array

!--------------------------------------------------------
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

module wave2d_mesh

  use constants

  implicit none

contains

  subroutine mesh_one_proc(iproc,SIMUL_TYPE,mesh_anchor, NEX, NEZ,nspec,LENGTH,HEIGHT,&
              x,z,nglob,ibool,xigll,zigll,nspecb,ibelm,MODEL_X1,MODEL_X2,&
              MODEL_Z1,MODEL_Z2,id,DEBUG)

    integer ::iproc
    integer ::SIMUL_TYPE
    double precision,intent(in) :: mesh_anchor(2,2)
    integer,intent(in) :: NEX, NEZ
    integer,intent(in) ::nspec
    integer,dimension(NGLLX,NGLLZ,nspec) ::ibool
    integer ,intent(in) ::nglob
    integer,dimension(3,nglob),intent(out) ::id
    double precision, dimension(nglob) ::x,z
    double precision,dimension(NGLLX)::xigll
    double precision,dimension(NGLLZ)::zigll
    integer, dimension(4) ::nspecb
    integer,dimension(:,:)::ibelm
    double precision::MODEL_X1,MODEL_X2,MODEL_Z1,MODEL_Z2

    double precision :: LENGTH, HEIGHT
    logical ::DEBUG

    !local variable
    double precision ::x1,x2,z1,z2
    integer ispec,ib,i,j,k,iglob,iglob1,itime,ix,iz,ii,jj
    integer N_EQ,n_fix

    ispec = 0
    iglob = 0
    nspecb(:) = 0
    N_EQ=0
    n_fix=0


    !-----------------------------------------------------------
    !LENGTH=abs(mesh_anchor(1,2)-mesh_anchor(1,1))
    !HEIGHT=abs(mesh_anchor(2,2)-mesh_anchor(2,1))
    ! loop over all elements
    do iz = 1,NEZ
       do ix = 1,NEX
          ispec = ispec+1

          ! evenly spaced achors between 0 and 1
          x1 = LENGTH*dble(ix-1)/dble(NEX)+mesh_anchor(1,1)
          x2 = LENGTH*dble(ix)/dble(NEX)+mesh_anchor(1,1)
          z1 = HEIGHT*dble(iz-1)/dble(NEZ)+mesh_anchor(2,1)
          z2 = HEIGHT*dble(iz)/dble(NEZ)+mesh_anchor(2,1)

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
                x(iglob1) = 0.5*(1.-xigll(i))*x1+0.5*(1.+xigll(i))*x2
                z(iglob1) = 0.5*(1.-zigll(j))*z1+0.5*(1.+zigll(j))*z2
             enddo
          enddo

          !print *,ibool(1,1,1), ibool(1,2,1),ibool(2,1,1)

          ! if boundary element
          if (ix.eq.1) then         
            if(abs(x1-MODEL_X1).lt.TOL)then
              nspecb(1) = nspecb(1) + 1
              ibelm(nspecb(1),1) = ispec
            endif
          endif
          if (ix.eq.NEX) then
            if(abs(x2-MODEL_X2).lt.TOL)then
              nspecb(2) = nspecb(2) + 1
              ibelm(nspecb(2),2) = ispec
            endif
          endif
          if (iz.eq.1) then
            if(abs(z1-MODEL_Z1).lt.TOL)then
              nspecb(3) = nspecb(3) + 1
              ibelm(nspecb(3),3) = ispec
            endif
          endif
          if (iz.eq.NEZ) then
            if(abs(z2-MODEL_Z2).lt.TOL)then
              if(iproc.eq.0) print *,'z2',z2,MODEL_Z2,mesh_anchor(2,2)
              nspecb(4) = nspecb(4) + 1
              ibelm(nspecb(4),4) = ispec
           endif
          endif
          ! end loop over elements
       end do
    end do
    
  if(iglob.ne.nglob) then
        print *,iglob,nglob
        stop 'Error in getting ibool'
  end if

  call set_ID_array(nspec,nglob,SIMUL_TYPE,N_EQ,n_fix,NEX,NEZ,&
                ibool,id,DEBUG)
    !call set_IEN_array()
    !call set_LM_array()

    ! estimate the time step
    !dh = HEIGHT/dble((NGLLZ-1)*NEZ)
    !if(dh > LENGTH/dble((NGLLX-1)*NEX)) dh = LENGTH/dble((NGLLX-1)*NEX)
    !c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
    !time_step = 0.2*dh/c
    !print *
    !print *,'time step estimate from courant = 0.2: ',sngl(time_step),' seconds'
    !print *,'  actual time step: ',sngl(DT),' seconds'
   
  end subroutine mesh_one_proc

  !-------------------------------------------------------
  subroutine set_ID_array(nspec,nglob,SIMUL_TYPE,N_EQ,n_fix,NEX,NEZ,&
                ibool,ID,DEBUG)

    integer ::nspec,nglob
    integer ::SIMUL_TYPE
    integer ::N_EQ
    integer ::n_fix
    integer ::NEX,NEZ
    integer,dimension(3,nglob),intent(out) ::ID
    integer,dimension(NGLLX,NGLLZ,nspec) ::ibool
    logical ::DEBUG

    !local variable
    integer ix,iz,i,j,ispec

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
  subroutine set_model_property(nspec,nglob,ibool,x,z,rho,kappa,mu)
    integer ::nspec,nglob
    integer,dimension(NGLLX,NGLLZ,nspec)::ibool
    double precision,dimension(nglob) ::x,z
    double precision, dimension(NGLLX,NGLLZ,nspec)::rho,kappa,mu

    !local variables
    integer :: ispec, i, j, iglob

    ! properties
    do ispec = 1,nspec
       !  get center of element
       iglob = ibool(NGLLX/2,NGLLZ/2,ispec)
       do j = 1,NGLLZ
          do i = 1,NGLLX
             if(z(iglob) <= 30000.0) then
                ! crust
                rho(i,j,ispec) = 2600
                kappa(i,j,ispec) = 5.2e+10
                mu(i,j,ispec) = 2.6e+10
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

end module wave2d_mesh

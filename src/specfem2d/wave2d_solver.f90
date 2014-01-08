module wave2d_solver

  implicit none
  include "constants.h"

contains
!---------------------------------------------

  subroutine solver(nsources, sglob, samp, nreceivers, rglob, ramp, &
        DT, NSTEP, rank, nproc, comm)

    use wave2d_variables
    use wave2d_define_der_matrices
    !use wave2d_mesher
    !use solve_eg

    integer, intent(in) :: nsources, nreceivers
    integer, intent(in) :: sglob(:), rglob(:)
    double precision :: samp(:,:,:), ramp(:,:,:)

    !mpi env var
    integer, intent(in) :: rank, nproc, comm

    !mass matrix
    double precision :: mass_global(NGLOB), mass_local

    !stiffness matrix
    double precision :: Stiff(NGLOB, NGLOB)

    !motion info
    double precision, dimension(3, NGLOB) :: displ,veloc,accel
    double precision, dimension(1, 3*NGLOB) :: rhs
   
    !time marching info
    double precision :: DT
    integer :: NSTEP
    double precision :: deltat,deltatover2,deltatsqover2,deltatsq

    integer ispec,ib,i,j,k,iglob, iglob1, iglob2,itime,ix,iz, itime1, itime2, itime0
    integer isource, irec, icomp, isave
    character(len=100) :: filename


		!integer :: ilocal,ilocal1,ilocal2
    !integer :: iglobal1,iglobal2

    !integer :: ii,jj,itemp,jtemp

    !double precision :: M_B(3,2*NGLLSQUARE)
    !double precision,dimension(2*NGLOB,2*NGLOB) :: MK,inv_MK,temp_MK

    !double precision, dimension(3,2*NGLLSQUARE) :: temp_DB
    !double precision, dimension(2*NGLLSQUARE,2*NGLLSQUARE) :: temp_BDB
    double precision, dimension(1,3*NGLOB) :: d_global,v_global,a_global
    !double precision :: fac
    
    character(len=10) :: rhs3_name,rhs1_name,rhs2_name

    ! grid points per wavelength estimation 
    c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
    print *, 'number of grid points per wavelength for P: ', floor(hdur*c/dh)
    c = sqrt(RIGIDITY/DENSITY)
    print *, 'number of grid points per wavelength for S:', floor(hdur*c/dh)

    !NINT = NSTEP/NSAVE
    !if (NINT * NSAVE > NSTEP) stop 'NSTEP should equal to NINT * NSAVE'

    ! calculate the global mass matrix once and for all
    mass_global(:) = 0.
    do ispec = 1,NSPEC
      do j = 1,NGLLZ
        do i = 1,NGLLX
          mass_local = wxgll(i)*wzgll(j)*rho(i,j,ispec)*jacobian(i,j,ispec)
          iglob = ibool(i,j,ispec)
          mass_global(iglob) = mass_global(iglob) + mass_local
        end do
      end do
    end do

    ! time marching parameters
    deltat = DT
    deltatsq = deltat*deltat
    deltatover2 = deltat/2.
    deltatsqover2 = deltat*deltat/2.

    ! initialize
    displ(:,:) = 0.
    veloc(:,:) = 0.
    accel(:,:) = 0.

    !Init stiff matrix
    Stiff(:,:)=0.0
    call init_stiff_matrix_PSV(Stiff)
    call init_stiff_matrix_SH(Stiff)

    !MK=mass_global+Stiff
    !if(IM_TRUE)then
    !  MK=S_BETA*deltat*deltat*Stiff
    !  do i=1,NGLOB
    !    MK(i,i)=MK(i,i)+mass_global(i)
    !    MK(i+NGLOB,i+NGLOB)=MK(i+NGLOB,i+NGLOB)+mass_global(i)
    !  enddo

    !  print *,"inverse the matrix MK now..."
    !  if(SIMUL_TYPE.eq.2) then
    !    call cholsl(N_EQ,MK,inv_MK) 
    !  endif

    !  if(DEBUG_SOLVER)then
    !    temp_MK=matmul(MK,inv_MK)
    !    open(30,file="MK*INV_MK")
    !    do i=1,N_EQ
    !      do j=1,N_EQ
    !        write(30,*) i,j,temp_MK(i,j)
    !      enddo
    !    enddo
    !  endif
    !endif
    !stop

    if(DEBUG_SOLVER)then
      open(21,file="stiff")
      do i=1,NGLOB
        write(21,*) Stiff(i,:)
      enddo
    endif

    !
    ! START TIME MARCHING
    !
    do itime = 1, NSTEP

      ! 'predictor' update displacement using finite-difference time scheme (Newmark)
      do i = 1, 3
        displ(i,:) = displ(i,:) + deltat*veloc(i,:) + deltatsqover2*(1-2*S_BETA)*accel(i,:)
        veloc(i,:) = veloc(i,:) + (1-S_ALPHA)*deltat*accel(i,:)
        accel(i,:) = 0.
      enddo

      !print *,"displ prediction"
      !print *,displ(3,:)
      !print *,"veloc prediction"
      !print *,veloc(3,:)
      !print *,"accel prediction"
      !print *,accel(3,:)

      !accel_total: accelration vector
      !ensemble the real problem: drag the real value out
      do i=1,NGLOB
        do j=1,3 !3 dimension, x-y-z (P-SV and SH)
          if(ID(j,i)/=0) then
            d_global(1,ID(j,i))=displ(j,i)
          endif
        enddo
      enddo

      !N_EQ=max(max(ID(:,:)))

      rhs(1:1,1:N_EQ)=transpose(-matmul(Stiff(1:N_EQ,1:N_EQ),transpose(d_global(1:1,1:N_EQ))))


      if(DEBUG_SOLVER)then
        open(22,file="rhs1")
        write(22,*)rhs(:,:)
      endif
      !-----
      !print *, "begin a:"
      !print *, "--------"
      !print *, "accel 1"
      !print *, accel(1,:)
      !print *, "--------"
      !print *, "accel 3"
      !print *, accel(3,:)
      !print *, "--------end"
      !-----

      !
      ! force boundary conditions
      !
      !if(FORCE_FLAG) then
      !  print *, "Add external force"
      !  !surface
      !  ibb=4
      !  i=NGLLZ
      !  do ib=1,nspecb(ibb)
      !    if(ib==1)then
      !      j1=2; j2=NGLLX
      !    elseif(ib==nspecb(ibb))then
      !      j1=1; j2=NGLLX-1
      !    else
      !      j1=1; j2=NGLLX
      !    endif
      !    ispec = ibelm(ibb,ib)
      !    do j=j1,j2
      !      if(LMX(j,i,ispec)/=0) then
      !        rhs(1,LMX(j,i,ispec))=rhs(1,LMX(j,i,ispec))+force_x*wxgll(j)*jacobianb(ibb,j,ib)
      !      endif
      !      if(LMZ(j,i,ispec)/=0) then
      !        rhs(1,LMZ(j,i,ispec))=rhs(1,LMZ(j,i,ispec))+force_z*wxgll(j)*jacobianb(ibb,j,ib)
      !      endif
      !    enddo
      !  enddo
      !endif

      ! absorbing boudary condition
      if(ABSORB_FLAG) then
        print *,"add absorb boundary condition"
        call add_absorb_bc(rhs, veloc)
      enddo

      if(DEBUG_SOLVER)then
        write(rhs3_name,'(i4)') itime
        rhs3_name=adjustl(rhs3_name)
        rhs3_name='rhs2_'//trim(rhs3_name)
        write(*,*) rhs3_name
        open(23,file=rhs3_name)
        do i=1,NGLOB
          write(23,*) i,"x",rhs(1,i)
          write(23,*) i,"z",rhs(1,i+NGLOB)
        enddo
      endif
      !end absorb boudary condition
      endif

      
      ! add source
      !print *, "samp 1"
      !print *, samp(itime,:,1);
      !print *, "end samp"
      !print *,accel(:,sglob(1))
      !print *, "end accel source"

      if(SOURCE_FLAG) then
        print *,"add source"
        do isource = 1, nsources
          iglob = sglob(isource)
          !print *,"source,P_SV,ID",ID(1,iglob),ID(2,iglob)
          if(iglob.ne.0)then
            rhs(1,ID(1,iglob)) = rhs(1,ID(1,iglob)) + samp(itime,1,isource) 
            rhs(1,ID(2,iglob)) = rhs(1,ID(2,iglob)) + samp(itime,2,isource)
            rhs(1,ID(3,iglob)) = rhs(1,ID(3,iglob)) + samp(itime,3,isource)
          endif
        enddo
      endif

      if(DEBUG_SOLVER)then
        write(rhs3_name,'(i4)') itime
        rhs3_name=adjustl(rhs3_name)
        rhs3_name='rhs3_'//trim(rhs3_name)
        write(*,*) rhs3_name
        open(24,file=rhs3_name)
        do i=1,NGLOB
          write(24,*) i,"x",rhs(1,i)
          write(24,*) i,"z",rhs(1,i+NGLOB)
        enddo
        close(24)
      endif

      !print *,accel(:,sglob(1))
      !print *, "end accel source"

      if(DAMP_FLAG)then
        print *,"add damping"
        !damping the RHS
        do i=1,NGLOB
          if(ID(1,i)/=0)then
            rhs(1,ID(1,i))=rhs(1,ID(1,i))-damp_factor*sign1(veloc(1,i))*abs(rhs(1,ID(1,i)))
          endif
          if(ID(2,i)/=0)then
            rhs(1,ID(2,i))=rhs(1,ID(2,i))-damp_factor*sign1(veloc(3,i))*abs(rhs(1,ID(2,i)))
          endif
        enddo
      endif


      call assemble_MPI_mass(mass_global,nglob, &
                    ninterface,  max_ibool_interfaces_size, &
                    ibool_interfaces, nibool_interfaces,
                    my_neighbours)
      call assemble_MPI_vector(rhs,nglob,ID, &
                      ninterface, max_ibool_interfaces_size,&
                      ibool_interfaces, nibool_interfaces, &
                      tab_requests_send_recv, &
                      buffer_send_faces_vector,&
                      buffer_recv_faces_vector,&
                      my_neighbour)
      ! above here accel(:) are actually the RHS!
      ! divide by the mass matrix
      if(.not.IM_TRUE) then
        do i=1,NGLOB
          do j=1,3 !3 dimension
          if(ID(j,i)/=0)then
            accel(j,i) = rhs(1,ID(j,i))/mass_global(i)
          endif
        end do
      endif

      !if(IM_TRUE) then
      !  !do d
      !  a_global=transpose(matmul(inv_MK,transpose(rhs)))
      !  accel(1,:)=a_global(1,1:NGLOB)
      !  accel(3,:)=a_global(1,(NGLOB+1):(2*NGLOB))
      !endif

      !if(SIMUL_TYPE==2) then
      !  do i=1,NGLOB
      !    accel(1,i)=rhs(1,i)/mass_global(i)
      !    accel(3,i)=rhs(1,i+NGLOB)/mass_global(i)
      !  enddo
      !endif


      !print *,"------------------------------"
      !print *,"after modification,accel"
      !do i=1,NGLOB
      !  print *,i,ID(2,i),accel(3,i)
      !enddo
      !print *, accel(3,:)

      !print *, " true accel"

      ! `corrector' update velocity
      do  i = 1, 3
        veloc(i,:) = veloc(i,:) + S_ALPHA*deltat*accel(i,:)
        displ(i,:) = displ(i,:) + S_BETA*deltatsq*accel(i,:)
      enddo

      do irec = 1, nreceivers
        if(rglob(irec).ne.0)then
          ramp(itime,:,irec) = displ(:,rglob(irec))
        endif
      enddo

      ! save snapshots
      !print *, "time_step=",itime
      !if (mod(itime, NSAVE) == 0) then

      !  write(*,*) 'itime = ', itime
      !  write(filename,'(a,i4.4)') trim(out_dir)//'snapshot_',itime
!
 !       open(unit = 11, file = trim(filename),status = 'unknown',iostat=ios)
 !       if (ios /= 0) stop 'Error writing snapshot to disk'
 !       do iglob = 1, NGLOB
 !         write(11,'(5e12.3)') sngl(x(iglob)/LENGTH), sngl(z(iglob)/LENGTH), (sngl(displ(j,iglob)),j=1,3) 
 !       enddo
 !       close(11)
!
!      endif

      !stop
    end do ! end time loop

  end subroutine solver

!---------------------------------------------
  integer function sign1(x)
    double precision :: x
    if(x>0) then
      sign1=1
    elseif(x<0) then
      sign1=-1
    else
      sign1=0
    endif
    return
  end

!-----------------------------------------------
  subroutine init_stiff_matrix_PSV(Stiff)

    double precision :: Stiff(:,:)

    integer ispec,i,j,ii,jj
    integer :: ilocal
    double precision :: M_B(3,2*NGLLSQUARE)
    double precision, dimension(3,2*NGLLSQUARE) :: temp_DB
    double precision, dimension(2*NGLLSQUARE,2*NGLLSQUARE) :: temp_BDB
    double precision :: fac
    !calculate the stiffness matrix. once and for all

    do ispec = 1,NSPEC

      ! first double loop over GLL 
      ! compute and store B
      Stiff_e(:,:)=0.

      do j = 1,NGLLZ
        do i = 1,NGLLX
          !iglob2 = ibool(i,j,ispec)
          M_B(:,:) = 0.
          !compute M_B at every node in one element
          !integal them together: M_B(i,j)*jacobian*wgll
          do jj = 1,NGLLZ
            do ii=1,NGLLX
            !calculate B matrix for every node: M_Bode
              ilocal = NGLLX*(jj-1)+ii
              if(jj.eq.j) then
                !M_B(1,ilocal)=dxidxl(i,j,ispec)*hprime_xx(ii,i)
                M_B(1,ilocal)=dxidx(i,j,ispec)*hprime_xx(ii,i)
                !+dgammadxl(i,j,ispec)*hprime_zz(jj,j)
              end if
              if(ii.eq.i) then
                !M_B(2,ilocal+NGLLS) = dxidzl(i,j,ispec)*hprime_xx(ii,i)+dgammadzl(i,j,ispec)*hprime_zz(jj,j)
                M_B(2,ilocal+NGLLS) = dgammadz(i,j,ispec)*hprime_zz(jj,j)
              end if
            enddo
          enddo
          M_B(3,1:NGLLS) = M_B(2,(NGLLS+1):2*NGLLS)
          M_B(3,(NGLLS+1):2*NGLLS)= M_B(1,1:NGLLS)

           !print *,"-------------------"
           !print *, "dxidx"
           !print *, dxidx(:,:,ispec)
           !print *,"hprime_xx"
           !print *,hprime_xx(:,:)

           !print *,"-------------------"
           !print *,"-------------------"
           !print *, NGLLS
          !print *,"ispec=",ispec
           !print *,"i=:",i,"j=:",j
           !print *,"M_B:"
           !do itemp=1,3
           !  print *,M_B(itemp,1:NGLLS)
           !  print *,"-------------------"
           !  print *,M_B(itemp,(NGLLS+1):2*NGLLS)
           !print *,"-------------------"
           !enddo
           !print *,"-------------------"

          kappal = kappa(i,j,ispec)
          mul = mu(i,j,ispec)
          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2.*mul

          D = reshape((/lambdalplus2mul,lambdal,dble(0.),lambdal,lambdalplus2mul,dble(0.),dble(0.),dble(0.),mul/),(/3,3/))
          !print *,"-------------------"
          !print *,"D"
          !print *,D
          !print *,"-------------------"
            
          temp_DB=matmul(D,M_B)
          temp_BDB=matmul(transpose(M_B),temp_DB)

          fac=jacobian(i,j,ispec)*wxgll(i)*wzgll(j)

          Stiff_e=Stiff_e+temp_BDB*fac
        enddo
      enddo

      !ensemble the K matrix
      call ensem(ispec, Stiff_e, Stiff, "PSV")

    enddo
  end subroutine init_stiff_matrix_PSV

  subroutine init_stiff_matrix_SH(Stiff)

    double precision :: Stiff(:,:)

    integer ispec,i,j,ii,jj
    integer :: ilocal
    double precision :: M_B(2,NGLLSQUARE)
    double precision, dimension(2,NGLLSQUARE) :: temp_DB
    double precision, dimension(NGLLSQUARE,NGLLSQUARE) :: temp_BDB
    double precision :: fac
    !calculate the stiffness matrix. once and for all

    do ispec = 1,NSPEC

      ! first double loop over GLL 
      ! compute and store B
      Stiff_e(:,:)=0.

      do j = 1,NGLLZ
        do i = 1,NGLLX
          !iglob2 = ibool(i,j,ispec)
          M_B(:,:) = 0.
          !compute M_B at every node in one element
          !integal them together: M_B(i,j)*jacobian*wgll
          do jj = 1,NGLLZ
            do ii=1,NGLLX
            !calculate B matrix for every node: M_Bode
              ilocal = NGLLX*(jj-1)+ii
              if(jj.eq.j) then
                !M_B(1,ilocal)=dxidxl(i,j,ispec)*hprime_xx(ii,i)
                M_B(1,ilocal)=dxidx(i,j,ispec)*hprime_xx(ii,i)
                !+dgammadxl(i,j,ispec)*hprime_zz(jj,j)
              end if
              if(ii.eq.i) then
                !M_B(2,ilocal+NGLLS) = dxidzl(i,j,ispec)*hprime_xx(ii,i)+dgammadzl(i,j,ispec)*hprime_zz(jj,j)
                M_B(2,ilocal) = dgammadz(i,j,ispec)*hprime_zz(jj,j)
              end if
            enddo
          enddo

           !print *,"-------------------"
           !print *, "dxidx"
           !print *, dxidx(:,:,ispec)
           !print *,"hprime_xx"
           !print *,hprime_xx(:,:)

           !print *,"-------------------"
           !print *,"-------------------"
           !print *, NGLLS
          !print *,"ispec=",ispec
           !print *,"i=:",i,"j=:",j
           !print *,"M_B:"
           !do itemp=1,3
           !  print *,M_B(itemp,1:NGLLS)
           !  print *,"-------------------"
           !  print *,M_B(itemp,(NGLLS+1):2*NGLLS)
           !print *,"-------------------"
           !enddo
           !print *,"-------------------"

          kappal = kappa(i,j,ispec)
          mul = mu(i,j,ispec)
          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2.*mul

          D = reshape((/mul,dble(0.),dble(0.),mul/),(/2,2/))
          !print *,"-------------------"
          !print *,"D"
          !print *,D
          !print *,"-------------------"
            
          temp_DB=matmul(D,M_B)
          temp_BDB=matmul(transpose(M_B),temp_DB)

          fac=jacobian(i,j,ispec)*wxgll(i)*wzgll(j)

          Stiff_e=Stiff_e+temp_BDB*fac
        enddo
      enddo

      !ensemble the K matrix
      call ensem(ispec, Stiff_e, Stiff, "SH")

    enddo
  end subroutine init_stiff_matrix_PSV

!---------------------------------------------
  subroutine ensem(ispec, Stiff_e, Stiff, wave_type)

    use wave2d_variables

    integer :: ispec
    double precision :: Stiff_e(:,:), Stiff(:,:)
    character(len=*) :: wave_type
    
    integer :: i,j
    integer :: ii,jj
    integer :: iglobal1_x,iglobal1_z,iglobal2_x,iglobal2_z
    integer :: ilocal1,ilocal2

    !over one column 
    if(trim(wave_type).eq."PSV")then
      do j=1,NGLLZ
        do i=1,NGLLX
          !over one row
          do jj=1,NGLLZ
            do ii=1,NGLLX
!            iglobal1=ibool(i,j,ispec)
!            ilocal1=(j-1)*NGLLX+i
!            iglobal2=ibool(ii,jj,ispec)
!            ilocal2=(jj-1)*NGLLX+ii

!            Stiff(iglobal1,iglobal2)=Stiff(iglobal1,iglobal2)+Stiff_e(ilocal1,ilocal2)
!            Stiff(iglobal1+NGLOB,iglobal2+NGLOB)=Stiff(iglobal1+NGLOB,iglobal2+NGLOB)+Stiff_e(ilocal1+NGLLS,ilocal2+NGLLS)
!            Stiff(iglobal1+NGLOB,iglobal2)=Stiff(iglobal1+NGLOB,iglobal2)+Stiff_e(ilocal1+NGLLS,ilocal2)
!            Stiff(iglobal1,iglobal2+NGLOB)=Stiff(iglobal1,iglobal2+NGLOB)+Stiff_e(ilocal1,ilocal2+NGLLS)
              ilocal2=(jj-1)*NGLLX+ii
              ilocal1=(j-1)*NGLLX+i
              iglobal1_x=ID(1,ibool(i,j,ispec))
              iglobal2_x=ID(1,ibool(ii,jj,ispec))
              iglobal1_z=ID(3,ibool(i,j,ispec))
              iglobal2_z=ID(3,ibool(ii,jj,ispec))
            
              if((iglobal1_x/=0).and.(iglobal2_x/=0))then
                Stiff(iglobal1_x,iglobal2_x)=Stiff(iglobal1_x,iglobal2_x)+Stiff_e(ilocal1,ilocal2)
              endif
              if((iglobal1_x/=0).and.(iglobal2_z/=0))then
                Stiff(iglobal1_x,iglobal2_z)=Stiff(iglobal1_x,iglobal2_z)+Stiff_e(ilocal1,ilocal2+NGLLS)
              endif
              if((iglobal1_z/=0).and.(iglobal2_x/=0))then
                Stiff(iglobal1_z,iglobal2_x)=Stiff(iglobal1_z,iglobal2_x)+Stiff_e(ilocal1+NGLLS,ilocal2)
              endif
              if((iglobal1_z/=0).and.(iglobal2_z/=0))then
                Stiff(iglobal1_z,iglobal2_z)=Stiff(iglobal1_z,iglobal2_z)+Stiff_e(ilocal1+NGLLS,ilocal2+NGLLS)
              endif
            enddo
          enddo
        enddo
      enddo
    else if(trim(wave_type).eq."SH")then
      do j=1,NGLLZ
        do i=1,NGLLX
          !over one row
          do jj=1,NGLLZ
            do ii=1,NGLLX
              ilocal2=(jj-1)*NGLLX+ii
              ilocal1=(j-1)*NGLLX+i
              iglobal1_y=ID(2,ibool(i,j,ispec))
              iglobal2_y=ID(2,ibool(ii,jj,ispec))
              if((iglobal1_y/=0).and.(iglobal2_y/=0))then
                Stiff(iglobal1_y,iglobal2_y)=Stiff(iglobal1_y,iglobal_y)+Stiff_e(ilocal1,ilocal2)
              endif
            enddo
          enddo
        enddo
      enddo
    else
      print *, "Ensemble Erro: type incorrect"
      stop
    endif
            
  end subroutine ensem

  subroutine add_absorb_bc(rhs, veloc)

      !absorbing boundary conditions
      integer :: j1,j2, ib, ibb
      double precision :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight

      do ibb = 1, 3
        if (ibb == 1) then
          i = 1; nx = -1.; nz = 0.
        else if (ibb == 2) then
          i = NGLLX;  nx = 1.;  nz = 0.
        else if (ibb == 3) then
          i = 1;  nx = 0.; nz = -1.
        end if
        do ib = 1,nspecb(ibb)
          if (ibb == 1 .or. ibb == 2) then
            j1 = 1; j2 = NGLLZ
          else if (ib == 1) then ! left corner element
            j1 = 2; j2 = NGLLX 
          else if (ib == nspecb(ibb)) then ! right corner element
            j1 = 1; j2 = NGLLX-1  
          else
            j1 = 1; j2 = NGLLX
          endif
          ispec = ibelm(ibb,ib)

          do j = j1, j2
            if (ibb == 1 .or. ibb == 2) then
              iglob = ibool(i,j,ispec)
              rho_vp = dsqrt(rho(i,j,ispec)*(kappa(i,j,ispec)+FOUR_THIRDS*mu(i,j,ispec)))
              rho_vs = dsqrt(rho(i,j,ispec)*mu(i,j,ispec))
            else
              iglob = ibool(j,i,ispec)
              rho_vp = dsqrt(rho(j,i,ispec)*(kappa(j,i,ispec)+FOUR_THIRDS*mu(j,i,ispec)))
              rho_vs = dsqrt(rho(j,i,ispec)*mu(j,i,ispec))
            endif

            vx = veloc(1,iglob)
            vy = veloc(2,iglob)
            vz = veloc(3,iglob)

            vn = nx*vx+nz*vz

            tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
            ty = rho_vs*vy
            tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

            weight = jacobianb(ibb,j,ib)*wzgll(j)

            !if(LMX(i,j,ispec)/=0) then
              !rhs(1,LMX(i,j,ispec)) = rhs(1,LMX(i,j,ispec)) - tx*weight
              rhs(1,ID(1,iglob)) = rhs(1,ID(1,iglob)) - tx*weight
              rhs(1,ID(2,iglob)) = rhs(1,ID(2,iglob)) - ty*weight
            !endif
            !if(LMZ(i,j,ispec)/=0) then
              !rhs(1,iglob+NGLOB) = rhs(1,LMZ(i,j,ispec)) - tz*weight
              rhs(1,ID(3,iglob)) = rhs(1,ID(3,iglob)) - tz*weight
            !endif
          enddo
        enddo
  end subroutine add_absorb_bc
!---------------------------------------------

end module wave2d_solver

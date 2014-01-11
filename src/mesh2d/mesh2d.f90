program main
!mesh main program

  !use wave2d_variables
  use mesh2d_par
  use wave2d_mesh
  use build_mesh_interface_subs

  implicit none
  
  !local variables
  integer ::i,j

    
  !********* PROGRAM STARTS HERE *********************


  !==========
  !start mesh
  call read_parameters(nex,nez,model_x1,model_x2,model_z1,model_z2,nproc_x,nproc_z,&
        simul_type,debug,density0,incompressibility0,rigidity0)

!calculate more parameters based on the input parameters 
  length=dabs(model_x1-model_x2)/nproc_x
  height=dabs(model_z1-model_z2)/nproc_z
  nex_one_proc=NEX/nproc_x
  nez_one_proc=NEZ/nproc_z
  nspec=nex_one_proc*nez_one_proc
  nglob=NGLLX*NGLLZ+(nex_one_proc-1)*NGLLZ*(NGLLX-1)+&
        (NGLLX*NGLLZ+(nex_one_proc-1)*NGLLZ*(NGLLX-1)- &
         (1+(NGLLX-1)*nex_one_proc))*(nez_one_proc-1)
  max_interface_size=MAX(nex_one_proc,nez_one_proc)
  max_interface_size_node=MAX(nex_one_proc*(NGLLX-1)+1,nez_one_proc*(NGLLZ-1)+1)
!allocate the array
 allocate(x(nglob))
 allocate(z(nglob))
 allocate(ibool(NGLLX,NGLLZ,nspec))
 allocate(id(3,nglob))
 allocate(ibelm(max_interface_size,4))
 allocate(ibool_interfaces(max_interface_size_node,8))
 allocate(rho(NGLLX,NGLLZ,nspec))
 allocate(kappa(NGLLX,NGLLZ,nspec))
 allocate(mu(NGLLX,NGLLZ,nspec))

 print *,"max",max_interface_size


  call calcu_xigll_zigll(xigll,zigll)

!do each processor

  do i=1, nproc_z
    do j=1, nproc_x
       iproc=(i-1)*nproc_x+j-1
       mesh_anchor(1,1)=MODEL_X1+(j-1)*length
       mesh_anchor(2,1)=MODEL_Z1+(i-1)*height
       mesh_anchor(1,2)=mesh_anchor(1,1)+length
       mesh_anchor(2,2)=mesh_anchor(2,1)+height
!      print *,'anchor1',mesh_anchor(:,1),mesh_anchor(:,2)
       call init_variables(nglob,nspec,x,z,ibool,id,nspecb,max_interface_size,ibelm,my_neighbour,nibool_interfaces,&
              max_interface_size_node,ibool_interfaces,rho,kappa,mu)

       call mesh_one_proc(iproc,simul_type,mesh_anchor,nex_one_proc,nez_one_proc,nspec,length,height,&
              x,z,nglob,ibool,xigll,zigll,nspecb,ibelm,model_x1,model_x2,&
              model_z1,model_z2,id,debug)

       call build_mesh_interface(iproc,nex_one_proc,nez_one_proc,mesh_anchor,&
              ibool, max_interface_size_node, my_neighbour, &
              nibool_interfaces,ibool_interfaces,model_x1,model_x2,&
              model_z1,model_z2,nproc_x)
       call set_model_property(nspec,nglob,ibool,x,z,rho,kappa,mu)
       call save_mesh(iproc,nex_one_proc,nez_one_proc,nglob,ibool,x,z,nspec,nspecb,max_interface_size,ibelm, &
              MAX_INTERFACE,my_neighbour,nibool_interfaces,max_interface_size_node,ibool_interfaces,id,&
              rho,kappa,mu)
    end do
  enddo

end program main 
 

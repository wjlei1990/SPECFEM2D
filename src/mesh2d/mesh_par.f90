module constants
    include "constant.h"
end module constants


module mesh2d_par

use constants
integer ::nex,nez,nex_one_proc,nez_one_proc
double precision ::model_x1,model_x2,model_z1,model_z2
integer ::nproc_x,nproc_z
integer ::simul_type
logical ::debug

integer ::iproc
double precision::length,height
double precision,dimension(2,2)::mesh_anchor
integer ::nglob,nspec
double precision,dimension(:),allocatable ::x,z
integer,dimension(:,:,:),allocatable::ibool
integer,dimension(:,:),allocatable ::id
integer,dimension(4) ::nspecb
integer,dimension(:,:),allocatable::ibelm
integer,dimension(MAX_INTERFACE) ::my_neighbour,nibool_interfaces
integer,dimension(:,:),allocatable::ibool_interfaces
double precision:: density0,incompressibility0,rigidity0
double precision,dimension(:,:,:),allocatable::rho,kappa,mu


integer ::max_interface_size,max_interface_size_node
double precision,dimension(NGLLX) ::xigll
double precision,dimension(NGLLZ) ::zigll
end module

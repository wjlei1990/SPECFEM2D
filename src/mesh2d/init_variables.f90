subroutine    init_variables(nglob,nspec,x,z,ibool,id,nspecb,max_interface_size,ibelm,my_neighbour,nibool_interfaces,&
                  max_interface_size_node,ibool_interfaces,rho,kappa,mu)
use constants
implicit none
integer ::nglob,nspec
double precision,dimension(nglob)  ::x,z
integer,dimension(NGLLX,NGLLZ,nspec) ::ibool
integer,dimension(2,nglob) ::id
integer,dimension(4) :: nspecb
integer ::max_interface_size
integer,dimension(4,max_interface_size) ::ibelm
integer,dimension(MAX_INTERFACE) ::my_neighbour,nibool_interfaces
integer ::max_interface_size_node
integer,dimension(max_interface_size_node,8)  ::ibool_interfaces
double precision, dimension(NGLLX,NGLLZ,nspec) ::rho,kappa,mu

x(:)=HUGE
z(:)=HUGE
ibool(:,:,:)=HUGE_INT
id(:,:)=HUGE_INT
nspecb(:)=HUGE_INT
ibelm(:,:)=HUGE_INT
my_neighbour(:)=HUGE_INT
nibool_interfaces(:)=HUGE_INT
ibool_interfaces(:,:)=HUGE_INT
rho(:,:,:)=HUGE
kappa(:,:,:)=HUGE
mu(:,:,:)=HUGE
end subroutine

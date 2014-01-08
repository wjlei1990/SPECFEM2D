module wave2d_variables

 include "constants.h" 

!simulation type
!1)dynamic relaxation 2)wave simulation
  logical :: SOURCE_FLAG 
  logical :: FORCE_FLAG
  logical :: DAMP_FLAG
  logical :: ABSORB_FLAG
  logical :: IM_TRUE

  double precision :: S_ALPHA
  double precision :: S_BETA

! sound velocity
  double precision :: c

! material properties
  double precision :: kappal,mul,lambdal,lambdalplus2mul
  double precision, dimension(:,:,:), allocatable :: rho,kappa,mu
  !double precision, dimension(3,3) :: D


!ID array, IEN array, LM array
  integer, dimension(:,:), allocatable :: ID
  !integer, dimension(NGLLX,NGLLZ,NSPEC) :: IEN
  !integer, dimension(NGLLX,NGLLZ,NSPEC) :: LMX,LMZ
  integer :: N_EQ


!mesh info
  integer :: NEX, NEZ
  integer :: NELE, NGLLSQUARE, NGLOB

!some matrix
  !double precision, dimension(2*NGLLSQUARE,2*NGLLSQUARE) :: Stiff_e
  !double precision, dimension(2*NGLOB,2*NGLOB) :: Stiff
  

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz

!jacobian
  double precision, dimension(:,:,:) :: jacobian

! anchors
!  double precision, dimension(NSPEC) :: z1,z2
!  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(:), allocatable :: x,z

! local to global numbering
  integer, dimension(:,:,:), allocatable :: ibool

! boundary elements: 1-xmin, 2-xmax, 3-zmin, 4-zmaxs
  integer, dimension(:,:), allocatable     :: ibelm

!mesh info
  integer :: nglob, nspec, nspecb(4), ninterface
  integer, allocatable :: my_neighbour(:), nibool_interfaces(:)
  integer, allocatable :: ibool_interfaces(:,:)

! Jacobian matrix and Jacobian
!  double precision dxidxl,dxidzl,dgammadxl,dgammadzl,jacobianl
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: dxidx,dxidz,dgammadx,dgammadz,jacobian
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: dxdxi,dzdxi,dxdgamma,dzdgamma
!  double precision, dimension(4,NGLL,NELE) :: jacobianb
!  integer :: nspecb(4)

! absorbing boundary conditions
!  integer :: j1,j2, ib, ibb
!  double precision :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight

! Mass matrix
  !double precision mass_local
  !double precision, dimension(:), allocatable :: mass_global

! time marching
  !double precision :: DT
  double precision :: dh,time_step

! displacement, velocity and acceleration
 
! plotting
  double precision, dimension(:), allocatable :: norm

! space derivatives
  double precision tempx1l,tempx2l,tempy1l,tempy2l,tempz1l,tempz2l
  double precision fac1,fac2,hp1,hp2
  double precision dsxdxl,dszdxl,dsydxl,dsydzl,dsxdzl,dszdzl 
  double precision sigma_xx,sigma_xy,sigma_xz,sigma_zx,sigma_zy,sigma_zz
  double precision, dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempy1,tempy2,tempz1,tempz2
  double precision, dimension(3,3) :: ds

! half duration of the source
  double precision hdur

! file open or write status variable
  integer ios

! number of time steps to store wavefield
  integer NINT 


end module wave2d_variables

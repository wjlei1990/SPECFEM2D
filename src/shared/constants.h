  integer, parameter :: NELE = MAX(NEX,NEZ)
  integer, parameter :: NSPEC = NEX*NEZ

  integer, parameter :: NGLL = MAX(NGLLX,NGLLZ)

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLZ
  integer, parameter :: NGLLS = NGLLX * NGLLZ

! number of global points
  integer, parameter :: NGLOB = ((NGLLX-1)*NEX + 1)*((NGLLZ-1)*NEZ +1)




! number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = 5

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

!
! CONSTANTS
!
! pi
  double precision, parameter :: PI = 3.141592653589793d+00
! 4/3
  double precision, parameter :: FOUR_THIRDS = 4.d0/3.d0

  double precision, parameter :: ONE_THIRD = 1.d0/3.d0

  double precision, parameter :: ONEOVERTWO = 0.5d0

  double precision, parameter :: EPS = 1.0d-35


! normalization factor of force
  double precision, parameter :: FNORM = 1.0d12


!force term
  double precision :: force_x=0
  double precision :: force_z=10000

!damp term
  double precision :: damp_factor=0.5

end module wave2d_constants

  subroutine calcu_xigll_zigll(xigll,zigll)

  use constants
  implicit none
 
! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = 0.d0
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = 0.d0

end subroutine calcu_xigll_zigll


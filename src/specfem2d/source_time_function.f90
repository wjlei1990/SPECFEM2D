  double precision function source_time_function(t,hdur)

  implicit none

! source decay rate (also change in source spectrum if needed)
  double precision, parameter :: decay_rate = 2.628d0
  double precision, parameter :: PI = 3.141592653589793d0

  double precision t,hdur

  double precision alpha
  double precision, external :: erf

  alpha = decay_rate/hdur

! Error function
!  source_time_function = 0.5d0*(1.0d0+erf(decay_rate*t/hdur))

! Gaussian (this one causes static offset at stations)
!  source_time_function = alpha*dexp(-alpha*alpha*t*t)/dsqrt(PI)

! Ricker wavelet
  source_time_function = -2.*(alpha**3)*t*dexp(-alpha*alpha*t*t)/dsqrt(PI)
  

  end function source_time_function

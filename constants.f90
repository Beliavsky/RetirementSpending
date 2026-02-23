module constants_mod
use kind_mod, only: dp
implicit none
private
public :: pi, log_two_pi, log_two, sqrt_two
real(kind=dp), parameter :: &
   pi         = 3.141592653589793238462643_dp, &
   log_two_pi = 1.837877066409345483560659_dp, &
   log_two    = 0.69314718055994529_dp, &
   sqrt_two   = 1.4142135623730951_dp
end module constants_mod

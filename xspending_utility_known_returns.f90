program xspend_utility_known_returns
use kind_mod, only: dp
implicit none

! deterministic returns (sigma = 0):
! choose constant annual withdrawal w from savings to maximize
!   sum_{t=1..n_years} log(pension + withdraw_t),
! where withdraw_t = min(w, wealth_t),
! and wealth evolves deterministically:
!   wealth_{t+1} = max(0, (wealth_t - withdraw_t) * gross),
! with gross = 1 + mu.
!
! for mu >= 0 (gross >= 1): analytic optimum is the largest feasible constant
! withdrawal that lasts exactly n_years (annuity-due payment). computed below
! and printed alongside the grid-search estimate.

integer, parameter :: n_years = 30
integer, parameter :: ngrid  = 2001

real(kind=dp), parameter :: w0 = 1000000.0_dp    ! initial savings/portfolio wealth at t=0 (before any withdrawals)
real(kind=dp), parameter :: pension = 40000.0_dp ! constant annual pension income added to consumption each year
real(kind=dp), parameter :: mu = 0.05_dp         ! mean annual arithmetic return of the portfolio (e[r])
real(kind=dp), parameter :: gross = 1.0_dp + mu

real(kind=dp) :: w_lo, w_hi, w, u, u_best, w_best
real(kind=dp) :: step
real(kind=dp) :: w_formula
integer :: j

w_formula = w_opt_formula(w0, gross, n_years)

w_lo = 0.0_dp
w_hi = w0
step = (w_hi - w_lo) / real(ngrid - 1, kind=dp)

u_best = -huge(1.0_dp)
w_best = w_lo

do j=1, ngrid
   w = w_lo + step * real(j-1, kind=dp)
   u = utility_for_w(w)
   if (u > u_best) then
      u_best = u
      w_best = w
   end if
end do

print "('n_years     : ', i0)", n_years
print "('w0          : ', f0.2)", w0
print "('pension     : ', f0.2)", pension
print "('mu          : ', f0.6)", mu
print "('sigma       : ', f0.6)", 0.0_dp
print "('gross       : ', f0.8)", gross
print "('w_opt_grid  : ', f0.2)", w_best
print "('w_opt_form  : ', f0.2)", w_formula
print "('diff        : ', f0.6)", w_best - w_formula
print "('w_grid/w0   : ', f0.6)", w_best / w0
print "('w_form/w0   : ', f0.6)", w_formula / w0
print "('u_opt_grid  : ', f0.6)", u_best

call terminal_stats(w_best, 'grid')
call terminal_stats(w_formula, 'form')

contains

elemental real(kind=dp) function w_opt_formula(w0_in, g, n) result(w)
real(kind=dp), intent(in) :: w0_in, g
integer, intent(in) :: n
real(kind=dp) :: denom

! annuity-due payment that exactly exhausts wealth after n withdrawals:
! if g == 1: w = w0/n
! else:      w = w0 / sum_{t=0}^{n-1} g^{-t} = w0*(g-1)*g^{n-1}/(g^n - 1)
!
! note: this is the log-utility optimum for the deterministic case when g >= 1.
! for 0 < g < 1, the true optimum can be different.

if (n <= 0) then
   w = 0.0_dp
   return
end if

if (abs(g - 1.0_dp) <= 1.0e-12_dp) then
   w = w0_in / real(n, kind=dp)
else
   denom = g**n - 1.0_dp
   if (abs(denom) <= 1.0e-30_dp) then
      w = w0_in / real(n, kind=dp)
   else
      w = w0_in * (g - 1.0_dp) * g**(n-1) / denom
   end if
end if
end function w_opt_formula

real(kind=dp) function utility_for_w(w) result(u)
real(kind=dp), intent(in) :: w
integer :: t
real(kind=dp) :: wealth, withdraw, cons

wealth = w0
u = 0.0_dp

do t=1, n_years
   withdraw = min(w, wealth)
   cons = pension + withdraw
   if (cons <= 0.0_dp) then
      u = -huge(1.0_dp)
      return
   end if
   u = u + log(cons)

   wealth = (wealth - withdraw) * gross
   if (wealth < 0.0_dp) wealth = 0.0_dp
end do
end function utility_for_w

subroutine terminal_stats(w, label)
real(kind=dp), intent(in) :: w
character(len=*), intent(in) :: label
integer :: t
real(kind=dp) :: wealth, withdraw

wealth = w0
do t=1, n_years
   withdraw = min(w, wealth)
   wealth = (wealth - withdraw) * gross
   if (wealth < 0.0_dp) wealth = 0.0_dp
end do

print "('terminal wealth (', a, ') : ', f0.2)", trim(label), wealth
end subroutine terminal_stats

end program xspend_utility_known_returns
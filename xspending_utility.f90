program xspend_opt
use kind_mod,   only: dp
use random_mod, only: random_normal, random_seed_init
implicit none

! simulate log utility vs. constant dollar withdrawal rate w from savings each year
! optional output:
!   run with -p (or --print) to print the 6-line summary block
!   for each candidate spending rate evaluated in the grid searches.
!
! new:
!   at the end, print a table ONLY for the coarse grid values (ngrid1 rows)
!   with columns: w/w0, eu, mean terminal wealth, p(wealth=0 at end)

integer, parameter :: n_years = 30        ! planning horizon in years (number of annual withdrawals/consumption periods)
integer, parameter :: n_paths = 100000    ! number of monte carlo simulated return paths used to estimate expected utility
integer, parameter :: ngrid1  = 16        ! number of candidate spending values w evaluated on the coarse grid search
integer, parameter :: ngrid2  = 101       ! number of candidate w values evaluated per refinement grid search
integer, parameter :: nrefine = 2         ! number of refinement passes (each pass re-centers a finer grid around current best w)

real(kind=dp), parameter :: w0 = 1000000.0_dp   ! initial savings/portfolio wealth at t=0 (before any withdrawals)
real(kind=dp), parameter :: pension = 40000.0_dp ! constant annual pension income added to consumption each year

real(kind=dp), parameter :: mu = 0.05_dp        ! mean annual arithmetic return of the portfolio (e[r])
real(kind=dp), parameter :: sigma = 0.15_dp     ! annual standard deviation of arithmetic return of the portfolio (sd[r])

real(kind=dp), allocatable :: gross(:,:)

real(kind=dp) :: w_best, u_best
real(kind=dp) :: w_lo, w_hi, step
integer :: iseed

character(len=64) :: arg1

! store coarse-grid results only
real(kind=dp), allocatable :: w_coarse(:), eu_coarse(:), mean_wt_coarse(:), p_end0_coarse(:)


if (pension <= 0.0_dp) then
   print *, 'error: pension must be > 0 for log utility.'
   stop 1
end if

allocate(gross(n_years, n_paths))
allocate(w_coarse(ngrid1), eu_coarse(ngrid1), mean_wt_coarse(ngrid1), p_end0_coarse(ngrid1))

iseed = 0
call random_seed_init(iseed, nburn=1000)
call draw_returns(gross)

w_lo = 0.0_dp
w_hi = 0.15_dp * w0

call grid_search(w_lo, w_hi, ngrid1, w_best, u_best, print_all = .false., &
   store_coarse=.true.)

step = (w_hi - w_lo) / real(ngrid1 - 1, kind=dp)
call refine_search(w_best, step, ngrid2, nrefine, w_best, u_best, print_all = .false.)

call print_inputs()
call print_opt_block(w_best, u_best)

call print_coarse_table()

deallocate(w_coarse, eu_coarse, mean_wt_coarse, p_end0_coarse)
deallocate(gross)

contains

subroutine print_inputs()
print "('n_years  : ', i0)", n_years
print "('n_paths  : ', i0)", n_paths
print "('w0       : ', f0.2)", w0
print "('pension  : ', f0.2)", pension
print "('mu       : ', f0.6)", mu
print "('sigma    : ', f0.6)", sigma
end subroutine print_inputs

subroutine draw_returns(g)
real(kind=dp), intent(out) :: g(:,:)
integer :: t, i
real(kind=dp) :: r, gr

do i=1, size(g,2)
   do t=1, size(g,1)
      r = mu + sigma * random_normal()
      gr = 1.0_dp + r
      if (gr < 0.0_dp) gr = 0.0_dp
      g(t,i) = gr
   end do
end do
end subroutine draw_returns

subroutine eval_w(w, eu, mean_wt, p_end0, p_ever0)
real(kind=dp), intent(in)  :: w
real(kind=dp), intent(out) :: eu, mean_wt, p_end0, p_ever0

integer :: i, t
real(kind=dp) :: wealth, withdraw, cons
real(kind=dp) :: u_path, u_sum, wt_sum
integer :: ruin_by_end, ruin_ever
logical :: ruined

u_sum = 0.0_dp
wt_sum = 0.0_dp
ruin_by_end = 0
ruin_ever = 0

do i=1, n_paths
   wealth = w0
   ruined = .false.
   u_path = 0.0_dp

   do t=1, n_years
      withdraw = min(w, wealth)
      cons = pension + withdraw
      u_path = u_path + log(cons)

      wealth = (wealth - withdraw) * gross(t,i)
      if (wealth <= 0.0_dp) then
         wealth = 0.0_dp
         ruined = .true.
      end if
   end do

   if (wealth <= 0.0_dp) ruin_by_end = ruin_by_end + 1
   if (ruined) ruin_ever = ruin_ever + 1

   u_sum = u_sum + u_path
   wt_sum = wt_sum + wealth
end do

eu = u_sum / real(n_paths, kind=dp)
mean_wt = wt_sum / real(n_paths, kind=dp)
p_end0 = real(ruin_by_end, kind=dp) / real(n_paths, kind=dp)
p_ever0 = real(ruin_ever, kind=dp) / real(n_paths, kind=dp)
end subroutine eval_w

subroutine print_candidate_block(w, eu, mean_wt, p_end0, p_ever0)
real(kind=dp), intent(in) :: w, eu, mean_wt, p_end0, p_ever0
print "('w        : ', f0.2)", w
print "('w/w0     : ', f0.6)", w / w0
print "('eu       : ', f0.6)", eu
print "('mean terminal wealth : ', f0.2)", mean_wt
print "('p(wealth=0 at end)  : ', f0.6)", p_end0
print "('p(ruin by year n)   : ', f0.6)", p_ever0
end subroutine print_candidate_block

subroutine print_opt_block(w, eu0)
real(kind=dp), intent(in) :: w, eu0
real(kind=dp) :: eu, mean_wt, p_end0, p_ever0

call eval_w(w, eu, mean_wt, p_end0, p_ever0)

print "('w_opt    : ', f0.2)", w
print "('w_opt/w0 : ', f0.6)", w / w0
print "('eu_opt   : ', f0.6)", eu0
print "('mean terminal wealth : ', f0.2)", mean_wt
print "('p(wealth=0 at end)  : ', f0.6)", p_end0
print "('p(ruin by year n)   : ', f0.6)", p_ever0
end subroutine print_opt_block

subroutine print_coarse_table()
integer :: k
print *, '---'
print *, 'coarse grid table:'
print "(a)", "w/w0        eu            mean_terminal_wealth     p(wealth=0 at end)"
do k=1, ngrid1
   print "(f10.6,1x,f12.6,1x,f20.2,1x,f16.6)", w_coarse(k)/w0, eu_coarse(k), mean_wt_coarse(k), p_end0_coarse(k)
end do
end subroutine print_coarse_table

subroutine grid_search(a, b, ngrid, w_best, u_best, print_all, store_coarse)
real(kind=dp), intent(in) :: a, b
integer, intent(in) :: ngrid
real(kind=dp), intent(out) :: w_best, u_best
logical, intent(in) :: print_all
logical, intent(in) :: store_coarse

integer :: j
real(kind=dp) :: w, eu, mean_wt, p_end0, p_ever0

u_best = -huge(1.0_dp)
w_best = a

do j=1, ngrid
   w = a + (b-a) * real(j-1, kind=dp) / real(ngrid-1, kind=dp)

   call eval_w(w, eu, mean_wt, p_end0, p_ever0)

   if (store_coarse) then
      w_coarse(j) = w
      eu_coarse(j) = eu
      mean_wt_coarse(j) = mean_wt
      p_end0_coarse(j) = p_end0
   end if

   if (print_all) then
      print *, '---'
      call print_candidate_block(w, eu, mean_wt, p_end0, p_ever0)
   end if

   if (eu > u_best) then
      u_best = eu
      w_best = w
   end if
end do
end subroutine grid_search

subroutine refine_search(w0_guess, step0, ngrid, niter, w_best, u_best, print_all)
real(kind=dp), intent(in) :: w0_guess, step0
integer, intent(in) :: ngrid, niter
real(kind=dp), intent(out) :: w_best, u_best
logical, intent(in) :: print_all

integer :: k
real(kind=dp) :: a, b, step, w_guess
real(kind=dp) :: wb, ub

w_guess = w0_guess
step = step0

do k=1, niter
   a = max(0.0_dp, w_guess - 5.0_dp*step)
   b = min(w0,      w_guess + 5.0_dp*step)

   if (print_all) then
      print *, '---'
      print "('refine iter: ', i0, ', range: ', f0.2, ' to ', f0.2)", k, a, b
   end if

   call grid_search(a, b, ngrid, wb, ub, print_all, store_coarse=.false.)
   w_guess = wb
   step = (b-a) / real(ngrid-1, kind=dp)
end do

w_best = w_guess
u_best = ub
end subroutine refine_search

end program xspend_opt


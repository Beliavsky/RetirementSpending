program xbellman
use kind_mod,   only: dp
use random_mod, only: random_normal, random_seed_init
implicit none

! simulate log utility for three spending rules:
!   1) fixed annual dollar withdrawal (optimize w by grid search)
!   2) adaptive annuity withdrawal (recomputed each year from current wealth)
!   3) bellman dp spending rule (method b): reuse k simulated shocks per year to
!      approximate e[v_{t+1}(w_next)] inside the bellman equation, on a wealth grid
!
! modification:
!   all rules liquidate remaining wealth in the last year (withdraw = wealth when t == n_years)

integer, parameter :: n_years = 30        ! planning horizon in years (number of annual withdrawals/consumption periods)
integer, parameter :: n_paths = 100000    ! number of monte carlo simulated return paths used to estimate expected utility
integer, parameter :: ngrid1  = 16        ! number of candidate spending values w evaluated on the coarse grid search
integer, parameter :: ngrid2  = 101       ! number of candidate w values evaluated per refinement grid search
integer, parameter :: nrefine = 2         ! number of refinement passes (each pass re-centers a finer grid around current best w)

real(kind=dp), parameter :: w0 = 1000000.0_dp     ! initial savings/portfolio wealth at t=0 (before any withdrawals)
real(kind=dp), parameter :: pension = 10000.0_dp  ! constant annual pension income added to consumption each year

real(kind=dp), parameter :: mu = 0.06_dp          ! mean annual arithmetic return of the portfolio (e[r])
real(kind=dp), parameter :: sigma = 0.15_dp       ! annual standard deviation of arithmetic return of the portfolio (sd[r])

! bellman dp settings (m=200, a=50, k=500)
integer, parameter :: n_dp_wealth_grid   = 200          ! number of wealth grid points (state grid size m)
integer, parameter :: n_dp_withdraw_grid = 50           ! number of withdrawal choices per wealth (action grid size a)
integer, parameter :: n_dp_shocks        = 500          ! number of return shocks used for dp expectations (k)
real(kind=dp), parameter :: dp_wealth_max = 4.0_dp*w0   ! max wealth represented on dp grid (values above are capped)
logical, parameter :: print_terminal = .false.
logical, parameter :: call_pct = .false.
real(kind=dp), allocatable :: gross(:,:)               ! gross(t,i) = 1 + r, clipped at 0

real(kind=dp) :: w_best, u_best
real(kind=dp) :: w_lo, w_hi, step
integer :: iseed

character(len=64) :: arg1

! store coarse-grid results only (for fixed-dollar w grid)
real(kind=dp), allocatable :: w_coarse(:), eu_coarse(:), mean_wt_coarse(:), p_end0_coarse(:)

! adaptive annuity rule evaluation
real(kind=dp) :: eu_annuity, mean_wt_annuity, p_end0_annuity, p_ever0_annuity

! dp policy and dp evaluation
real(kind=dp), allocatable :: withdraw_dp(:,:)         ! (n_dp_wealth_grid, n_years)
real(kind=dp) :: eu_dp, mean_wt_dp, p_end0_dp, p_ever0_dp

if (pension <= 0.0_dp) then
   print *, 'error: pension must be > 0 for log utility.'
   stop 1
end if

allocate(gross(n_years, n_paths))
allocate(w_coarse(ngrid1), eu_coarse(ngrid1), mean_wt_coarse(ngrid1), p_end0_coarse(ngrid1))
allocate(withdraw_dp(n_dp_wealth_grid, n_years))

iseed = 0
call random_seed_init(iseed, nburn=1000)
call draw_returns( &
   gross) ! g: gross return matrix g(t,i)

w_lo = 0.0_dp
w_hi = 0.15_dp * w0

call grid_search( &
   w_lo, &               ! a: lower bound of withdrawal grid
   w_hi, &               ! b: upper bound of withdrawal grid
   ngrid1, &             ! ngrid: number of grid points
   w_best, &             ! w_best: output best w
   u_best, &             ! u_best: output best expected utility
   print_all=.false., &  ! print_all: print per-candidate blocks
   store_coarse=.true.)  ! store_coarse: store table rows for coarse grid only

step = (w_hi - w_lo) / real(ngrid1 - 1, kind=dp)
call refine_search( &
   w_best, &          ! w0_guess: starting guess for best w
   step, &            ! step0: initial step size around w0_guess
   ngrid2, &          ! ngrid: number of grid points in each refinement
   nrefine, &         ! niter: number of refinement iterations
   w_best, &          ! w_best: output refined best w
   u_best, &          ! u_best: output refined best expected utility
   print_all=.false.) ! print_all: print refinement diagnostics

call eval_rule_annuity( &
   eu_annuity, &       ! eu: expected sum of log consumption
   mean_wt_annuity, &  ! mean_wt: mean terminal wealth
   p_end0_annuity, &   ! p_end0: prob(terminal wealth = 0)
   p_ever0_annuity)    ! p_ever0: prob(wealth hits 0 at any time)

call solve_bellman_dp( &
   withdraw_dp) ! withdraw_policy: dp policy table (wealth grid x time)

call eval_policy_dp( &
   withdraw_dp, &  ! withdraw_policy: dp policy table (wealth grid x time)
   eu_dp, &        ! eu: expected sum of log consumption
   mean_wt_dp, &   ! mean_wt: mean terminal wealth
   p_end0_dp, &    ! p_end0: prob(terminal wealth = 0)
   p_ever0_dp)     ! p_ever0: prob(wealth hits 0 at any time)

call print_inputs()
call print_opt_block( &
   w_best, &  ! w: optimal fixed-dollar withdrawal (years 1..n_years-1)
   u_best)    ! eu0: expected utility (objective) from the optimizer

call print_rule_block( &
   "adaptive annuity", & ! name: rule label
   eu_annuity, &         ! eu: expected utility
   mean_wt_annuity, &    ! mean_wt: mean terminal wealth
   p_end0_annuity, &     ! p_end0: prob(wealth=0 at end)
   p_ever0_annuity)      ! p_ever0: prob(ruin ever)

call print_rule_block( &
   "bellman dp", & ! name: rule label
   eu_dp, &        ! eu: expected utility
   mean_wt_dp, &   ! mean_wt: mean terminal wealth
   p_end0_dp, &    ! p_end0: prob(wealth=0 at end)
   p_ever0_dp)     ! p_ever0: prob(ruin ever)

if (call_pct) call print_coarse_table()

deallocate(withdraw_dp)
deallocate(w_coarse, eu_coarse, mean_wt_coarse, p_end0_coarse)
deallocate(gross)

contains

subroutine print_inputs()
! print model and dp settings used in the run
print "('inputs:')"
print "('n_years  : ', i0)", n_years
print "('n_paths  : ', i0)", n_paths
print "('w0       : ', f0.2)", w0
print "('pension  : ', f0.2)", pension
print "('mu       : ', f0.6)", mu
print "('sigma    : ', f0.6)", sigma
print "('dp grid  : wealth=', i0, ', withdraw=', i0, ', shocks=', i0)", n_dp_wealth_grid, n_dp_withdraw_grid, n_dp_shocks
print "('dp wmax  : ', f0.2, /)", dp_wealth_max
end subroutine print_inputs

subroutine draw_returns( &
   g) ! g: gross return matrix g(t,i)
! draw gross returns for all years and paths
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

subroutine eval_w( &
   w, &        ! w: constant annual withdrawal (dollars) for years 1..n_years-1
   eu, &       ! eu: expected sum of log consumption
   mean_wt, &  ! mean_wt: mean terminal wealth
   p_end0, &   ! p_end0: prob(terminal wealth = 0)
   p_ever0)    ! p_ever0: prob(wealth hits 0 at any time)
! simulate the fixed-dollar withdrawal rule for a given w, with last-year liquidation
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
      if (t == n_years) then
         withdraw = wealth
      else
         withdraw = min(w, wealth)
      end if

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

subroutine eval_rule_annuity( &
   eu, &       ! eu: expected sum of log consumption
   mean_wt, &  ! mean_wt: mean terminal wealth
   p_end0, &   ! p_end0: prob(terminal wealth = 0)
   p_ever0)    ! p_ever0: prob(wealth hits 0 at any time)
! simulate the adaptive annuity withdrawal rule (recomputed each year), with last-year liquidation
real(kind=dp), intent(out) :: eu, mean_wt, p_end0, p_ever0

integer :: i, t, rem
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
      if (t == n_years) then
         withdraw = wealth
      else
         rem = n_years - t + 1
         withdraw = withdraw_annuity( &
            wealth, & ! wealth: current wealth before withdrawal
            rem)      ! rem: years remaining including this year
      end if

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
end subroutine eval_rule_annuity

real(kind=dp) function withdraw_annuity( &
   wealth, & ! wealth: current wealth before withdrawal
   rem) &    ! rem: years remaining including this year
   result(wd)
! compute annuity-style withdrawal based on wealth and remaining years
real(kind=dp), intent(in) :: wealth
integer, intent(in) :: rem

real(kind=dp) :: a, factor, gross_bar
integer :: k

if (wealth <= 0.0_dp) then
   wd = 0.0_dp
   return
end if

if (rem <= 0) then
   wd = 0.0_dp
   return
end if

gross_bar = 1.0_dp + mu
if (gross_bar <= 0.0_dp) then
   wd = wealth
   return
end if

a = 0.0_dp
factor = 1.0_dp
do k=1, rem
   a = a + factor
   factor = factor / gross_bar
end do

if (a <= 0.0_dp) then
   wd = wealth
else
   wd = wealth / a
end if

if (wd < 0.0_dp) wd = 0.0_dp
if (wd > wealth) wd = wealth
end function withdraw_annuity

subroutine solve_bellman_dp( &
   withdraw_policy) ! withdraw_policy: dp policy table (wealth grid x time)
! solve bellman equation on a wealth grid and store the optimal withdrawal policy, with last-year liquidation
real(kind=dp), intent(out) :: withdraw_policy(:,:)  ! (n_dp_wealth_grid, n_years)

integer :: t, iw, ia, k
integer :: k_use
real(kind=dp) :: dw, wealth, withdraw, remain, total, best, best_withdraw
real(kind=dp) :: cont, wnext, u_now

real(kind=dp), allocatable :: v_next(:), v_cur(:)

k_use = min(n_dp_shocks, n_paths)

allocate(v_next(n_dp_wealth_grid), v_cur(n_dp_wealth_grid))

dw = dp_wealth_max / real(n_dp_wealth_grid - 1, kind=dp)

v_next = 0.0_dp

do t=n_years, 1, -1

   if (t == n_years) then
      do iw=1, n_dp_wealth_grid
         wealth = dw * real(iw-1, kind=dp)
         v_cur(iw) = log(pension + wealth)
         withdraw_policy(iw,t) = wealth
      end do
      v_next = v_cur
      cycle
   end if

   do iw=1, n_dp_wealth_grid
      wealth = dw * real(iw-1, kind=dp)

      best = -huge(1.0_dp)
      best_withdraw = 0.0_dp

      do ia=1, n_dp_withdraw_grid
         withdraw = wealth * real(ia-1, kind=dp) / real(n_dp_withdraw_grid - 1, kind=dp)
         u_now = log(pension + withdraw)

         remain = wealth - withdraw
         cont = 0.0_dp

         do k=1, k_use
            wnext = remain * gross(t,k)
            cont = cont + interp_linear( &
               v_next, & ! v: value function on wealth grid at t+1
               wnext, &  ! x: next wealth value to evaluate
               dw)       ! dx: wealth grid spacing
         end do

         cont = cont / real(k_use, kind=dp)
         total = u_now + cont

         if (total > best) then
            best = total
            best_withdraw = withdraw
         end if
      end do

      v_cur(iw) = best
      withdraw_policy(iw,t) = best_withdraw
   end do

   v_next = v_cur
end do

deallocate(v_next, v_cur)
end subroutine solve_bellman_dp

subroutine eval_policy_dp( &
   withdraw_policy, & ! withdraw_policy: dp policy table (wealth grid x time)
   eu, &              ! eu: expected sum of log consumption
   mean_wt, &         ! mean_wt: mean terminal wealth
   p_end0, &          ! p_end0: prob(terminal wealth = 0)
   p_ever0)           ! p_ever0: prob(wealth hits 0 at any time)
! simulate the dp policy (interpolated in wealth) on monte carlo paths, with last-year liquidation
real(kind=dp), intent(in)  :: withdraw_policy(:,:) ! (n_dp_wealth_grid, n_years)
real(kind=dp), intent(out) :: eu, mean_wt, p_end0, p_ever0

integer :: i, t
real(kind=dp) :: dw, wealth, withdraw, cons
real(kind=dp) :: u_path, u_sum, wt_sum
integer :: ruin_by_end, ruin_ever
logical :: ruined

dw = dp_wealth_max / real(n_dp_wealth_grid - 1, kind=dp)

u_sum = 0.0_dp
wt_sum = 0.0_dp
ruin_by_end = 0
ruin_ever = 0

do i=1, n_paths
   wealth = w0
   ruined = .false.
   u_path = 0.0_dp

   do t=1, n_years
      if (t == n_years) then
         withdraw = wealth
      else
         withdraw = interp_policy( &
            withdraw_policy(:,t), & ! p: withdrawal policy values on wealth grid for year t
            wealth, &              ! wealth: current wealth
            dw)                    ! dx: wealth grid spacing
         if (withdraw > wealth) withdraw = wealth
      end if

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
end subroutine eval_policy_dp

real(kind=dp) function interp_linear( &
   v, &  ! v: values on an equally spaced grid
   x, &  ! x: query point
   dx) & ! dx: grid spacing
   result(y)
! linear interpolation of v on an equally spaced grid for query x
real(kind=dp), intent(in) :: v(:), x, dx
integer :: n, i
real(kind=dp) :: xmax, pos, frac

n = size(v)
xmax = dx * real(n-1, kind=dp)

if (x <= 0.0_dp) then
   y = v(1)
   return
end if

if (x >= xmax) then
   y = v(n)
   return
end if

pos = x / dx
i = int(pos) + 1
if (i < 1) i = 1
if (i >= n) i = n - 1

frac = pos - real(i-1, kind=dp)
y = v(i) + frac * (v(i+1) - v(i))
end function interp_linear

real(kind=dp) function interp_policy( &
   p, &      ! p: policy values on wealth grid
   wealth, & ! wealth: current wealth
   dx) &     ! dx: wealth grid spacing
   result(a)
! interpolate a policy p(wealth) on the wealth grid, with wealth capped to dp_wealth_max
real(kind=dp), intent(in) :: p(:), wealth, dx
real(kind=dp) :: wcap

wcap = wealth
if (wcap < 0.0_dp) wcap = 0.0_dp
if (wcap > dp_wealth_max) wcap = dp_wealth_max

a = interp_linear( &
   p, &    ! v: grid values (policy)
   wcap, & ! x: capped wealth value
   dx)     ! dx: grid spacing

if (a < 0.0_dp) a = 0.0_dp
if (a > wealth) a = wealth
end function interp_policy

subroutine print_opt_block( &
   w, &   ! w: optimal fixed-dollar withdrawal (years 1..n_years-1)
   eu0)   ! eu0: expected utility at that w (from optimizer)
! print the summary block for the optimized fixed-dollar withdrawal rule
real(kind=dp), intent(in) :: w, eu0
real(kind=dp) :: eu, mean_wt, p_end0, p_ever0

call eval_w( &
   w, &       ! w: constant annual withdrawal for years 1..n_years-1
   eu, &      ! eu: expected utility (recomputed here)
   mean_wt, & ! mean_wt: mean terminal wealth
   p_end0, &  ! p_end0: prob(wealth=0 at end)
   p_ever0)   ! p_ever0: prob(ruin ever)

print "('results:')"
print "('rule      : constant spending')"
print "('w_opt     : ', f0.2)", w
print "('w_opt/w0  : ', f0.6)", w / w0
print "('eu_opt    : ', f0.6)", eu0
if (print_terminal) then
   print "('mean terminal wealth : ', f0.2)", mean_wt
   print "('p(wealth=0 at end)  : ', f0.6)", p_end0
   print "('p(ruin by year n)   : ', f0.6)", p_ever0
end if
end subroutine print_opt_block

subroutine print_rule_block( &
   name, &    ! name: rule label
   eu, &      ! eu: expected utility
   mean_wt, & ! mean_wt: mean terminal wealth
   p_end0, &  ! p_end0: prob(wealth=0 at end)
   p_ever0)   ! p_ever0: prob(ruin ever)
! print a summary block for a named spending rule (annuity or dp)
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: eu, mean_wt, p_end0, p_ever0

print "(/,'rule      : ', a)", trim(name)
print "('eu        : ', f0.6)", eu
if (print_terminal) then
   print "('mean terminal wealth : ', f0.2)", mean_wt
   print "('p(wealth=0 at end)  : ', f0.6)", p_end0
   print "('p(ruin by year n)   : ', f0.6)", p_ever0
end if
end subroutine print_rule_block

subroutine print_coarse_table()
! print the coarse-grid table for fixed-dollar withdrawals
integer :: k
print *, '---'
print *, 'coarse grid table:'
print "(a)", "w/w0        eu            mean_terminal_wealth     p(wealth=0 at end)"
do k=1, ngrid1
   print "(f10.6,1x,f12.6,1x,f20.2,1x,f16.6)", w_coarse(k)/w0, eu_coarse(k), mean_wt_coarse(k), p_end0_coarse(k)
end do
end subroutine print_coarse_table

subroutine grid_search( &
   a, &            ! a: lower bound of withdrawal grid
   b, &            ! b: upper bound of withdrawal grid
   ngrid, &        ! ngrid: number of grid points
   w_best, &       ! w_best: output best w
   u_best, &       ! u_best: output best expected utility
   print_all, &    ! print_all: print per-candidate blocks (unused here)
   store_coarse)   ! store_coarse: store coarse table rows if true
! grid search over constant withdrawal w to maximize expected utility
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

   call eval_w( &
      w, &       ! w: constant annual withdrawal
      eu, &      ! eu: expected utility
      mean_wt, & ! mean_wt: mean terminal wealth
      p_end0, &  ! p_end0: prob(wealth=0 at end)
      p_ever0)   ! p_ever0: prob(ruin ever)

   if (store_coarse) then
      w_coarse(j) = w
      eu_coarse(j) = eu
      mean_wt_coarse(j) = mean_wt
      p_end0_coarse(j) = p_end0
   end if

   if (print_all) then
      print *, '---'
   end if

   if (eu > u_best) then
      u_best = eu
      w_best = w
   end if
end do
end subroutine grid_search

subroutine refine_search( &
   w0_guess, &  ! w0_guess: starting guess for best w
   step0, &     ! step0: initial step size around w0_guess
   ngrid, &     ! ngrid: number of grid points in each refinement
   niter, &     ! niter: number of refinement iterations
   w_best, &    ! w_best: output refined best w
   u_best, &    ! u_best: output refined best expected utility
   print_all)   ! print_all: print refinement diagnostics
! refine the best w by repeated local grid searches around the current best
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

   call grid_search( &
      a, &                 ! a: lower bound for local grid
      b, &                 ! b: upper bound for local grid
      ngrid, &             ! ngrid: number of grid points
      wb, &                ! w_best: best w found in this refinement
      ub, &                ! u_best: best utility found in this refinement
      print_all, &         ! print_all: pass through
      store_coarse=.false.)! store_coarse: do not store refinement rows

   w_guess = wb
   step = (b-a) / real(ngrid-1, kind=dp)
end do

w_best = w_guess
u_best = ub
end subroutine refine_search

end program xbellman

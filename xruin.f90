program xruin
! 11/17/2021 09:03 PM loop over spending rate
! 11/17/2021 05:00 PM loop over leverage
! 11/17/2021 10:41 AM simulate the probability of ruin given spending rate and return distribution
use kind_mod    , only: dp
use ziggurat_mod, only: rnor, jsrset
use stats_mod   , only: mean, median, mode
implicit none
integer, parameter :: nsteps = 40, nsim = 10000, nnfut = 4, h = nsteps/nnfut, nspend = 3, nlev = 3
real(kind=dp)      :: avg_return, sd_return, spend, wealth, wealth_init, wealth_final(nsim), &
                      spent(nsim), portfolio_return, leverage
integer            :: i, ifut, ilev, isim, ispend, tsurv(nsim), nfut(nnfut)
real(kind=dp), parameter :: xsim = real(nsim,kind=dp), leverage_min = 0.0_dp, leverage_h = 0.5_dp, &
                      spend_min = 0.02_dp, spend_h = 0.01_dp
character (len=20) :: label_fut(nnfut)
avg_return  = 0.06_dp ! average annual stock market return
sd_return   = 0.15_dp ! standard deviation of stock market returns
wealth_init = 1.00_dp
write (*,"(*(a10))") "#sim","avg_ret","sd_ret"
write (*,"(i10,2f10.2)") nsim,avg_return,sd_return
nfut = [(ifut*h,ifut=1,nnfut)] ! time horizons for which survival probabilities shown
do ifut=1,nnfut
   write (label_fut(ifut),"('p',i0)") nfut(ifut) ! labels for survival probabilities
end do
write (*,"(/,a6,10x,2a7,2a13)") "spend","spent","spent","years_surv","years_surv"
write (*,"(a6,a10,2a7,4a13,*(a7))") "rate","leverage","median","mean","median", &
       "mean","wealth_avg","wealth_surv",(trim(label_fut(ifut)),ifut=1,nnfut)
do ispend=1,nspend
   do ilev=1,nlev
      spend = spend_min + (ispend-1)*spend_h ! spending rate
      call jsrset(123456)    ! reset the seed to get the same normal variates for each level of spending and leverage
      leverage = leverage_min + (ilev-1)*leverage_h ! allocation to stock market
      spent = 0.0_dp
      do isim=1,nsim
         wealth = wealth_init
         do i=1,nsteps ! loop over time periods
            portfolio_return = leverage * (avg_return + sd_return*rnor())
            wealth = (wealth - spend) * (1 + portfolio_return) ! withdraw "spend" at beginning of period, invest the remainder
            spent(isim) = spent(isim) + spend ! cumulative amount spent
            if (wealth < spend) exit ! ran out of money
         end do
         wealth_final(isim) = wealth
         tsurv(isim) = i ! time that wealth survived
      end do
      write (*,"(f6.3,f10.4,2f7.3,4(1x,f12.4),*(1x,f6.4))") spend,leverage,median(spent),mean(spent), &
        median(tsurv),sum(tsurv)/xsim,mean(wealth_final),mean(pack(wealth_final,tsurv > nsteps)), &
        [(count(tsurv > nfut(ifut)),ifut=1,nnfut)]/xsim
   end do
   write (*,*)
end do
end program xruin

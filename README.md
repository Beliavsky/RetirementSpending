### Retirement Spending Simulation
Simulation in Fortran 2003 of spending and asset allocation rules in retirement, motivated by a Morningstar [study](https://www.morningstar.com/lp/the-state-of-retirement-income), discussed in a Wall Street Journal article [The 4% Retirement Rule Is in Doubt](https://www.wsj.com/articles/the-4-retirement-rule-is-in-doubt-will-your-nest-egg-last-11636713035). Suppose that

(1) A retiree will withdraw annually the same amount from savings, adjusted for inflation, as long as she lives.

(2) Annual after-inflation stock market returns are normally distributed with known mean and standard deviation.

(3) The investor rebalances annually to have a constant fraction of savings in stocks.

Then the two decisions for the investor to make are how much to spend annually and what fraction of savings to keep in stocks.
The program simulates the probability of savings lasting N years given the spending rule and stock market allocation. 
Parameters describing spending and asset allocation and stock market returns can easily be changed.

Compile with ```gfortran -std=f2003 kind.f90 stats.f90 ziggurat.f90 xruin.f90```

Output:

```
      #sim   avg_ret    sd_ret
     10000      0.06      0.15

 spend            spent  spent   years_surv   years_surv
  rate  leverage median   mean       median         mean   wealth_avg  wealth_surv    p10    p20    p30    p40
 0.020    0.0000  0.800  0.800      41.0000      41.0000       0.2000       0.2000 1.0000 1.0000 1.0000 1.0000
 0.020    0.5000  0.800  0.800      41.0000      40.9867       1.7227       1.7285 1.0000 1.0000 0.9999 0.9966
 0.020    1.0000  0.800  0.796      41.0000      40.7936       7.1893       7.3493 1.0000 0.9987 0.9917 0.9782

 0.030    0.0000  0.990  0.990      33.0000      33.0000       0.0100       0.0000 1.0000 1.0000 1.0000 0.0000
 0.030    0.5000  1.200  1.179      41.0000      40.1650       0.9740       1.1111 1.0000 1.0000 0.9763 0.8749
 0.030    1.0000  1.200  1.168      41.0000      39.8441       5.4993       6.1338 1.0000 0.9905 0.9469 0.8963

 0.040    0.0000  0.960  0.960      24.0000      24.0000       0.0400       0.0000 1.0000 1.0000 0.0000 0.0000
 0.040    0.5000  1.600  1.419      40.0000      35.9466       0.3888       0.7836 1.0000 0.9873 0.7713 0.4828
 0.040    1.0000  1.600  1.473      41.0000      37.5720       3.9787       5.3301 0.9997 0.9539 0.8361 0.7455
 ```
 
 Row 6 of the table above means that if the investor spends 3% of the initial portfolio value, adjusted for inflation, and 
 invests 100% of savings in the stock market, which has average after-inflation returns of 6% with standard deviation of 15%,
 that the probability of the savings lasting 30 years (column `p30`) is 94.7%. If the annual spending rate is 4%, the 30-year survival probability
 falls to 83.6%. The investor should decide whether spending 33% more per year is worth a higher risk of running out of money. Using
 different return assumptions, Morningstar recommends 3.3% as a safe withdrawal rate.

A second program, compiled with `gfortran kind.f90 constants.f90 random.f90 xspending_utility.f90`, looks at the optimal spending rate to maximize utility given a pension and log utility from consumption. Results are
```
n_years  : 30
n_paths  : 100000
w0       : 1000000.00
pension  : 40000.00
mu       : .060000
sigma    : .150000
w_opt    : 61700.00
w_opt/w0 : .061700
eu_opt   : 340.978462
mean terminal wealth : 1086372.63
p(wealth=0 at end)  : .528560
p(ruin by year n)   : .528560
 ---
 coarse grid table:
w/w0        eu            mean_terminal_wealth     p(wealth=0 at end)
  0.000000   317.899042           5750150.55         0.000000
  0.010000   324.593314           4912691.54         0.000080
  0.020000   330.051741           4075946.68         0.006370
  0.030000   334.537913           3248618.80         0.048290
  0.040000   337.970108           2459301.87         0.152350
  0.050000   340.138186           1752782.05         0.313240
  0.060000   340.961972           1171009.75         0.497630
  0.070000   340.637581            733113.55         0.666390
  0.080000   339.556105            429206.71         0.793710
  0.090000   338.097728            234757.69         0.883090
  0.100000   336.547079            120459.92         0.939040
  0.110000   335.089376             58099.06         0.969930
  0.120000   333.808004             26086.01         0.986400
  0.130000   332.713070             10998.80         0.994210
  0.140000   331.792868              4491.18         0.997590
  0.150000   331.018658              1654.03         0.999060
```
To maximize utility a spending rate of 6.17% of initial wealth is recommended.


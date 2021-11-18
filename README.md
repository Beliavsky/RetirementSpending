### Retirement Spending Simulation
Fortran simulation of spending rules in retirement. Suppose that

(1) A retiree will withdraw annually the same amount from savings, adjusted for inflation, as long as he lives

(2) Annual after-inflation stock market returns are normally distributed with known mean and standard deviation.

(3) The investor rebalances to have a constant fraction of savings in stocks.

Then the two decisions for the investor to make are how much to spend annually and what fraction of savings to keep in stocks.
The program simulates the probability of savings lasting N years given the spending rule and stock market allocation.

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

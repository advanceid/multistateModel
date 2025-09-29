This aims to test if increase number of requested cores beyond number of chains in rstan::sampling will still decrease time

md.stan is a function that will take some time, it will be cited in Rstan.

rstan-multicore.R cites md.stan and tries to run it on 4 chains but cores parameter in may differ (1, 2, 4, 8, 16).
For each running, the time will be logged.

When cores <= chains, time of execution is proportionate to 1/cores. 
When cores > chains, further increase in cores does not help with reducing time.
# Compare time consumption running Rstan
# Whether cores over chains reduces time



# Library ----
library(rstan)
library(parallel)
print(getwd())



# Model preparation ----
ms <- stan_model('multistate/md.stan')

N <- 100000
x <- rnorm(N, 0, 10)
df <- list(N = N, x = x)



# Model1 ----
start1 <- Sys.time()
sampling(ms,
         df,
         chains = 4,
         cores = 1,
         iter = 20000)
end1 <- Sys.time()



# Model 2 ----
start2 <- Sys.time()
sampling(ms,
         df,
         chains = 4,
         cores = 2,
         iter = 20000)
end2 <- Sys.time()



# Model 3 ----
start3 <- Sys.time()
sampling(ms,
         df,
         chains = 4,
         cores = 4,
         iter = 20000)
end3 <- Sys.time()



# Model 4 ----
start4 <- Sys.time()
sampling(ms,
         df,
         chains = 4,
         cores = 8,
         iter = 20000)
end4 <- Sys.time()



# Model 5 ----
start5 <- Sys.time()
sampling(ms,
         df, 
         chains = 4,
         cores = 16,
         iter = 20000)
end5 <- Sys.time()



# Time difference ----
print('cores = 1')
end1 - start1
print('cores = 2')
end2 - start2
print('cores = 4')
end3 - start3
print('cores = 8')
end4 - start4
print('cores = 16')
end5 - start5
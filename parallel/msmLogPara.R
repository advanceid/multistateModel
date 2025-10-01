# Initiation ----
# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/')
setwd('multistateModel/parallel')
rm(list = ls())
set.seed(123)
start <- Sys.time()

library(magrittr)
library(parallel)
library(readr)
library(rstan)
library(tidyverse)




# 1. Compile and expose Stan functions ----
n_cores <- min(2, detectCores())
print(n_cores)
cl <- makeCluster(n_cores)
print(getwd())

stan_code <- read_file('stan.txt')

clusterExport(cl, 'stan_code')
clusterEvalQ(cl, {
  library(rstan)
  sm <- stan_model(model_code = stan_code)
  rstan::expose_stan_functions(sm)
})
print('compiled')
# save.image('demonstration/faster/myenvironment.RData')
# load('demonstration/faster/myenvironment.RData')
sessionInfo()




# 2. Data preparation ----
cohort <- read.csv('dataMil/cohortInit.csv')
history <- read.csv('dataMil/history.csv')

cohort %<>% dplyr::select(-state0)
history %<>% filter(!is.na(id)) %>% 
  left_join(cohort, by = 'id') %>%
  mutate(age0 = age + time0/365,
         age1 = age + time1/365,
         stateE = ifelse(time1 > 365*5, state0, state1),
         timeE = pmin(time1, 365*5),
         time0 = floor(time0),
         timeE = floor(timeE)) %>%
  filter(time0 < 365*5)

followRS <- history %>% 
  mutate(timeE = timeE + 0.5,
         dur = timeE - time0,
         timeE2 = time0 + pmin(dur, 1825)) %>% 
  filter(id < 400000)


N <- nrow(followRS)
P <- 2
X <- as.matrix(followRS[,c('sex', 'age0')])
# t0 <- followRS$time0
# t1 <- followRS$timeE2
t <- followRS$timeE2 - followRS$time0
s0 <- followRS$state0
s1 <- followRS$stateE
# inten <- rep(0.001, N)

df_fl <- list(N = N, P = P, X = X, t = t, s0 = s0, s1 = s1)



# # 3. Split data into chunks for parallelization ----
# chunk_ids <- cut(seq_len(N), n_cores, labels=FALSE)
# chunks <- lapply(1:n_cores, function(i) list(
#   s0 = s0[chunk_ids==i],
#   s1 = s1[chunk_ids==i],
#   t = t[chunk_ids==i],
#   # t1 = t1[chunk_ids==i],
#   # N  = sum(chunk_ids==i),
#   
#   a12 = inten[chunk_ids==i],
#   a13 = inten[chunk_ids==i],
#   a15 = inten[chunk_ids==i],
#   a21 = inten[chunk_ids==i],
#   a24 = inten[chunk_ids==i],
#   a25 = inten[chunk_ids==i],
#   a31 = inten[chunk_ids==i],
#   a34 = inten[chunk_ids==i],
#   a35 = inten[chunk_ids==i],
#   a42 = inten[chunk_ids==i],
#   a43 = inten[chunk_ids==i],
#   a45 = inten[chunk_ids==i]
# ))
# print('check2')
# 
# 
# 
# # 4. Check function working ----
# print('check3')
# 
# sm <- stan_model(model_code = stan_code)
# rstan::expose_stan_functions(sm)
# a <- msm_loglik(chunks[[1]]$s0, chunks[[1]]$s1, chunks[[1]]$t,
#                 chunks[[1]]$a12, chunks[[1]]$a13, chunks[[1]]$a15,
#                 chunks[[1]]$a21, chunks[[1]]$a24, chunks[[1]]$a25,
#                 chunks[[1]]$a31, chunks[[1]]$a34, chunks[[1]]$a35,
#                 chunks[[1]]$a42, chunks[[1]]$a43, chunks[[1]]$a45)
# 
# print(a)
# print('check: rstan in global environment')
# 
# 
# 
# 5. Function:likelihood ----
func <- function(ch){
  print(exists('parLapply'))
  loglik <- msm_loglik(ch$s0, ch$s1, ch$t,
                       ch$a12, ch$a13, ch$a15,
                       ch$a21, ch$a24, ch$a25,
                       ch$a31, ch$a34, ch$a35,
                       ch$a42, ch$a43, ch$a45)
  return(loglik)
}


ready <- Sys.time()
# 6. Bayesian update ----
logpara <- rep(-10, 12)
best_logpara <- logpara 
step_size <- 0.1 
n_iter <- 2000
best_loglik <- -Inf
loglik_current <- -Inf
samples <- matrix(NA, nrow = n_iter, ncol = length(logpara))

print('check4')

for (i in 1:n_iter) {
  print(paste0('Iteration ', i))
  if (i %% 100 == 0){print(Sys.time())}
  print(best_logpara)
  proposal <- best_logpara + rnorm(length(logpara), mean = 0, sd = step_size)
  print(proposal)
  theta <- matrix(rep(exp(proposal), N), nrow = N, byrow = TRUE)
  
  chunk_ids <- cut(seq_len(N), n_cores, labels=FALSE)
  chunks <- lapply(1:n_cores, function(i) list(
    s0 = s0[chunk_ids==i],
    s1 = s1[chunk_ids==i],
    t  = t[chunk_ids==i],
    # N  = sum(chunk_ids==i),
    
    a12 = theta[chunk_ids==i, 1],
    a13 = theta[chunk_ids==i, 2],
    a15 = theta[chunk_ids==i, 3],
    a21 = theta[chunk_ids==i, 4],
    a24 = theta[chunk_ids==i, 5],
    a25 = theta[chunk_ids==i, 6],
    a31 = theta[chunk_ids==i, 7],
    a34 = theta[chunk_ids==i, 8],
    a35 = theta[chunk_ids==i, 9],
    a42 = theta[chunk_ids==i, 10],
    a43 = theta[chunk_ids==i, 11],
    a45 = theta[chunk_ids==i, 12]
  ))

  loglik <- parLapply(cl, chunks, func)
  print(loglik)
  loglik_prop <- sum(unlist(loglik))
  print(paste0('Sum of likelihood: ', loglik_prop))
  
  alpha <- exp(loglik_prop - loglik_current)
  if (runif(1) < alpha) {
    best_logpara <- proposal
    loglik_current <- loglik_prop
  }
  
  samples[i, ] <- best_logpara
}

stopCluster(cl)
write.csv(samples, 'mcmc-sample.csv', row.names = F)

finish <- Sys.time()
print(start)
print(ready)
print(finish)
q()


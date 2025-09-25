# informative observation will bias the results.
# continuous observation, even interval, and random interval will not bias the 
#   transition intensities



# RSTAN transition model ----
# Run multistate model with rstan ----
rm(list = ls())
# setwd('D:/ASPIRE/')
setwd('demonstration/generalData')
print(getwd())



# Library ----
library(expm)
library(magrittr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(msm)
library(rstan)



# Data Loading ----
follow <- read.csv('store/hkd.csv') %>%
  arrange(id, day) %>%
  filter(id <= 100)



# MSM model ----
Q <- rbind(c(0, 1, 1, 0, 1),
           c(1, 0, 0, 1, 1),
           c(1, 0, 0, 1, 1),
           c(0, 1, 1, 0, 1),
           c(0,0,0,0,0))

model <- msm(data = follow,
             stateNew ~ day,
             subject = id,
             qmatrix = Q,
             gen.inits = TRUE,
             control = list(fnscale = 100))
model



# RSTAN preparation ----
CTMT <- stan_model("msmModelLOG.stan")



# continuous observation ----
followRS <- follow %>%
  mutate(day0 = day - 1) %>%
  rename(day1 = day, state1 = stateNew) %>%
  dplyr::select(id, age, sex, day0, day1, state0, state1)

N <- nrow(followRS)
P <- 2

X <- as.matrix(followRS[,c('sex', 'age')])
t0 <- followRS$day0
t1 <- followRS$day1
y0 <- followRS$state0
y1 <- followRS$state1

df_fl <- list(N = N, P = P, X = X, t0 = t0, t1 = t1,
              s0 = y0, s1 = y1)

fit1 <- sampling(CTMT,
                 data = df_fl,
                 iter = 2000,
                 warmup = 500,
                 chains = 4,
                 cores = 4,
                 seed = 321)
saveRDS(fit1, file = "sampling/fit1.rds")
pars=c('loga12', 'loga13', 'loga15', "loga21", "loga24", 'loga25',
       'loga31', 'loga34', 'loga35', 'loga42', 'loga43', 'loga45')
print(fit1, pars = pars, digits = 4)
q()



# even-interval observation ----
followRS <- follow %>%
  rename(day1 = day, state1 = stateNew) %>%
  filter(day1 %% 2 == 1) %>%
  group_by(id) %>%
  mutate(state0 = lag(state1, 1), day0 = lag(day1, 1)) %>%
  filter(!is.na(day0)) %>%
  dplyr::select(id, age, sex, day0, day1, state0, state1) %>%
  ungroup

N <- nrow(followRS)
P <- 2

X <- as.matrix(followRS[,c('sex', 'age')])
t0 <- followRS$day0
t1 <- followRS$day1
y0 <- followRS$state0
y1 <- followRS$state1

df_fl <- list(N = N, P = P, X = X, t0 = t0, t1 = t1,
              s0 = y0, s1 = y1)

fit2 <- sampling(CTMT,
                data = df_fl,
                iter = 500,
                warmup = 200,
                chains = 4,
                cores = 4,
                seed = 321)
saveRDS(fit2, file = "sampling/fit2.rds")
pars=c('a12', 'a13', 'a15', "a21", "a24", 'a25',
       'a31', 'a34', 'a35', 'a42', 'a43', 'a45')
print(fit2, pars = pars, digits = 4)



# random-interval observation ----
followRS <- follow %>%
  rename(day1 = day, state1 = stateNew) %>%
  rowwise() %>%
  mutate(select = runif(1, min = 0, max = 1)) %>%
  filter(select < 0.6) %>%
  group_by(id) %>%
  mutate(state0 = lag(state1, 1), day0 = lag(day1, 1)) %>%
  filter(!is.na(day0)) %>%
  dplyr::select(id, age, sex, day0, day1, state0, state1) %>%
  ungroup

N <- nrow(followRS)
P <- 2

X <- as.matrix(followRS[,c('sex', 'age')])
t0 <- followRS$day0
t1 <- followRS$day1
y0 <- followRS$state0
y1 <- followRS$state1

df_fl <- list(N = N, P = P, X = X, t0 = t0, t1 = t1,
              s0 = y0, s1 = y1)

fit3 <- sampling(CTMT,
                 data = df_fl,
                 iter = 500,
                 warmup = 200,
                 chains = 4,
                 cores = 4,
                 seed = 321)
saveRDS(fit3, file = "sampling/fit3.rds")
pars=c('a12', 'a13', 'a15', "a21", "a24", 'a25',
       'a31', 'a34', 'a35', 'a42', 'a43', 'a45')
print(fit3, pars = pars, digits = 4)



# selective observation ----
followRS <- follow %>% 
  filter(visitNew) %>%
  rename(day1 = day, state1 = stateNew) %>%
  group_by(id) %>%
  mutate(state0 = lag(state1, 1), day0 = lag(day1, 1)) %>%
  filter(!is.na(day0)) %>%
  dplyr::select(id, age, sex, day0, day1, state0, state1) %>%
  ungroup

N <- nrow(followRS)
P <- 2

X <- as.matrix(followRS[,c('sex', 'age')])
t0 <- followRS$day0
t1 <- followRS$day1
y0 <- followRS$state0
y1 <- followRS$state1

df_fl <- list(N = N, P = P, X = X, t0 = t0, t1 = t1, 
              s0 = y0, s1 = y1)

fit4 <- sampling(CTMT,
                 data = df_fl,
                 iter = 500,
                 warmup = 200,
                 chains = 4,
                 cores = 4,
                 seed = 321)
saveRDS(fit4, file = "sampling/fit4.rds")
pars=c('a12', 'a13', 'a15', "a21", "a24", 'a25', 
       'a31', 'a34', 'a35', 'a42', 'a43', 'a45')
print(fit4, pars = pars, digits = 4)
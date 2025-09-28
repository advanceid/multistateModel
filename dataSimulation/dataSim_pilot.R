# Preparation ----
rm(list = ls())
# setwd('D:/github/multistateModel')
setwd('dataSimulation')
print(getwd())
set.seed(908)



# Library ----
library(expm)
library(magrittr)
library(tidyverse)
library(parallel)



# States ----
# 1: C0B0
# 2: C0B1
# 3: C1B0
# 4: C1B1
# 5: DEATH
# 1---2
# |\ /|
# | 5 |
# |/ \|
# 3---4



# Background ----
# annual death rate ~7/1000
# 200k infection every year
# annual BSI case 10k 

# visit probability if previously visited: 10, 80, 60, 90, 100
# visit probability if previously not visited: 1, 40, 30, 50, 90

# false negative: 90%



# Parameters ----
N <- 10^6
nDays <- 3650

a12 <- 1*10^-5
a13 <- 3*10^-5  
a15 <- 3*10^-5

a21 <- 1*10^-4
a24 <- 8*10^-4
a25 <- 3*10^-5

a31 <- 2*10^-5
a34 <- 8*10^-5
a35 <- 5*10^-5

a42 <- 1*10^-4
a43 <- 2*10^-4
a45 <- 5*10^-4

tmat <- matrix(ncol = 5, byrow = T,
               c(-a12-a13-a15, a12, a13, 0, a15,
                 a21, -a21-a24-a25, 0, a24, a25,
                 a31, 0, -a31-a34-a35, a34, a35,
                 0, a42, a43, -a42-a43-a45, a45,
                 0, 0, 0, 0, 0))

# probTrans <- expm(tmat)
# probTrans
# probTrans <- data.frame(probTrans) %>%
#   rename(prob1 = X1, prob2 = X2, prob3 = X3, prob4 = X4, prob5 = X5) %>%
#   mutate(state0 = 1:5)
# 
# probObs <- matrix(c(1, 0, 0, 0, 0,
#                   0.05, 0.95, 0, 0, 0,
#                   0.05, 0, 0.95, 0, 0,
#                   0.01, 0.05, 0.05, 0.89, 0,
#                   0, 0, 0, 0, 1), ncol = 5, byrow = T)
# probObs <- data.frame(probObs) %>%
#   rename(probO1 = X1, probO2 = X2, probO3 = X3, probO4 = X4, probO5 = X5) %>%
#   mutate(stateNew = 1:5)



# cohort initiation ----
cohort <- data.frame(id = 1:N)
cohort %<>% mutate(
  sex = rbinom(n = N, size = 1, prob = 0.4) + 1, # 40% male, 60% female
  age = rnorm(N, 70, 10),
  age2 = rnorm(N, 30, 5) + 
    rbinom(N, 1, 0.60)*rnorm(N, 40, 10),
  age = ifelse(sex == 2, age, age2),
  state0 = sample(c(1, 2, 3), size = N, 
                  replace = T, prob = c(0.99, 0.001, 0.009)),
  procedure = sample(0:1, size = nrow(cohort), replace = T, c(0.99, 0.01)),
  stay = 0, visit = 0) %>%
  filter(between(age, 20, 100)) %>%
  dplyr::select(!age2)

# cohort update ----
write.csv(cohort, 'cohortInit.csv', row.names = F)



# Data simulation ----
# for(subj in 1:nrow(cohort)){
history <- data.frame(id = NA, time0 = NA, time1 = NA, state0 = NA, state1 = NA)
for(subj in 1:nrow(cohort)){
  id <- cohort$id[subj]
  t0 <- 0
  # tTrans <- 0
  state0 <- cohort$state0[subj]
  
  while((t0 < nDays)&(state0 < 5)){
    lambda <- -tmat[state0,state0]
    tTrans <- rexp(1, rate = lambda)
    
    Ptrans <- expm(tmat*tTrans)
    Ptrans <- Ptrans[state0, ]
    Ptrans[state0] <- 0
    Ptrans <- Ptrans/sum(Ptrans)
    state1 <- sample(1:5, size = 1, prob = Ptrans)
    
    temp <- data.frame(
      id = id,
      time0 = t0, 
      time1 = t0 + tTrans,
      state0 = state0,
      state1 = state1)
      
    history <- rbind(history, temp)
    
    state0 <- state1
    t0 <- t0 + tTrans
  }
}

write.csv(history, 'history.csv', row.names = F)

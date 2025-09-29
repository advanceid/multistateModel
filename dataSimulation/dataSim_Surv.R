# Preparation ----
rm(list = ls())
# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/multistateModel')
setwd('multistateModel/dataSimulation')
print(getwd())
set.seed(908)



# Library ----
library(expm)
library(magrittr)
library(tidyverse)
library(doParallel)



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



# Parameters ----
N <- 2*10^6
nDays <- 3650

a12 <- 1*10^-3
a13 <- 5*10^-5  
a15 <- 5*10^-6

a21 <- 5*10^-3
a24 <- 2*10^-2
a25 <- 5*10^-4

a31 <- 5*10^-4
a34 <- 2*10^-3
a35 <- 5*10^-5

a42 <- 4*10^-2
a43 <- 4*10^-2
a45 <- 7*10^-4

tmat <- matrix(ncol = 5, byrow = T,
               c(-a12-a13-a15, a12, a13, 0, a15,
                 a21, -a21-a24-a25, 0, a24, a25,
                 a31, 0, -a31-a34-a35, a34, a35,
                 0, a42, a43, -a42-a43-a45, a45,
                 0, 0, 0, 0, 0))



# cohort initiation ----
cohort <- data.frame(id = 1:(N))
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
cohort <- cohort[1:(N/2), ]
write.csv(cohort, 'cohortInit.csv', row.names = F)



# History Function ----
func <- function(cohort){
  history <- data.frame(id = NA, time0 = NA, time1 = NA, state0 = NA, state1 = NA)
  for(subj in 1:nrow(cohort)){
    id <- cohort$id[subj]
    t0 <- 0
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
  history <- history[-1, ]
  return(history)
}



#### Parallelization ####
print('check')
start <- Sys.time()
ncores <- as.integer(Sys.getenv("PBS_NCPUS", unset = detectCores()))
ncores <- min(ncores, 30)
print(ncores)

cl <- makeCluster(ncores)
clusterExport(cl, varlist = c('nDays', 'tmat'))
clusterEvalQ(cl, {library(expm)})

chunks <-   cohorts <- split(cohort, cut(cohort$id, ncores, labels = FALSE))
history <- parLapply(cl, chunks, func)
history <- do.call(rbind, history)

end <- Sys.time()
end - start
stopCluster(cl)

write.csv(history, 'history.csv', row.names = F)

# Preparation ----
rm(list = ls())
# setwd('U:/dataSimulation')
# setwd('D:/ASPIRE/')
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

probTrans <- expm(tmat)
probTrans
probTrans <- data.frame(probTrans) %>%
  rename(prob1 = X1, prob2 = X2, prob3 = X3, prob4 = X4, prob5 = X5) %>%
  mutate(state0 = 1:5)

probObs <- matrix(c(1, 0, 0, 0, 0,
                  0.05, 0.95, 0, 0, 0,
                  0.05, 0, 0.95, 0, 0,
                  0.01, 0.05, 0.05, 0.89, 0,
                  0, 0, 0, 0, 1), ncol = 5, byrow = T)
probObs <- data.frame(probObs) %>%
  rename(probO1 = X1, probO2 = X2, probO3 = X3, probO4 = X4, probO5 = X5) %>%
  mutate(stateNew = 1:5)


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
write.csv(cohort, 'store/cohortDay000000.csv', row.names = F)


#### starting from middle ####
existing <- dir('store', full.names = T)
print(existing)
if (length(existing) > 2){
  existing <- existing[length(existing) - 1]
  breaking <- str_extract(existing, '[0-9]+')
  breaking <- as.numeric(breaking)
} else{
  breaking = 1
  existing <- existing[1]
}
cohort <- read.csv(existing)



#### Parallelization ####
start <- Sys.time()
ncores <- as.integer(Sys.getenv("PBS_NCPUS", unset = detectCores()))
ncores <- min(ncores, 100)
print(ncores)
cl <- makeCluster(ncores)



#### Update ####
updateC <- function(cohort = cohort){
  library(tidyverse)
  library(magrittr)
  
  cohort <- cohort %>% left_join(probTrans, by = 'state0') %>%
    rowwise() %>%
    mutate(day = i,
           age = age + 1/365,
           rollT = runif(1, 0, 1),
           stateNew = 1 + (rollT > prob1) + (rollT > prob1 + prob2)  + 
             (rollT > prob1 + prob2 + prob3) + 
             (rollT > prob1 + prob2 + prob3 + prob4),
           stay = ifelse(stateNew == state0, stay + 1, 0),
           
           probVisit = ifelse(visit == 0,
                              recode(stateNew, `1` = 0.1, `2` = 40, `3` = 20,
                                     `4` = 50, `5` = 90),
                              recode(stateO, `1` = 0.1, `2` = 60, `3` = 40,
                                     `4` = 80, `5` = 100)),
           rollV = runif(1, 0, 1),
           visitNew = rollV < probVisit/100) %>%
    left_join(probObs, by = 'stateNew') %>%
    mutate(rollO = runif(1, 0, 1),
      stateO = 1 + (rollO > probO1) + (rollO > probO1 + probO2)  +
             (rollO > probO1 + probO2 + probO3) +
             (rollO > probO1 + probO2 + probO3 + probO4)) %>%
    ungroup %>% dplyr::select(
      -c(prob1:prob5, probO1:probO5, rollT, rollV, rollO))
  
  return(cohort)
}

# Execute in parallel
for(i in breaking:nDays){
  print(i)
  clusterExport(cl, varlist = c("probTrans", 'probObs', 'i'))
  cohorts <- split(cohort, cut(cohort$id, ncores, labels = FALSE))
  cohorts <- parLapply(cl, cohorts, updateC)

  cohort <- do.call(rbind, cohorts)
  write.csv(cohort, row.names = F,
            paste0('store/cohortDay', str_pad(i, 7, pad = '0'), '.csv'))
  if (length(unique(cohort$id)) != nrow(cohort)) break
  cohort %<>% mutate(state0 = stateNew, visit = visitNew)
}
stopCluster(cl)
end <- Sys.time()

print(end - start)



# #### Comparison ####
# start2 <- Sys.time()
# for(i in 1:2){
# cohort <- cohort %>% left_join(probTrans, by = 'state0') %>%
#   rowwise() %>%
#   mutate(day = i,
#          age = age + 1/365,
#          rollT = runif(1, 0, 1),
#          stateNew = 1 + (rollT > prob1) + (rollT > prob1 + prob2)  + 
#            (rollT > prob1 + prob2 + prob3) + 
#            (rollT > prob1 + prob2 + prob3 + prob4),
#          stay = ifelse(stateNew == state0, stay + 1, 0),
#          
#          probVisit = ifelse(visit == 0,
#                             recode(stateNew, `1` = 1, `2` = 40, `3` = 30,
#                                    `4` = 50, `5` = 90),
#                             recode(stateNew, `1` = 10, `2` = 80, `3` = 60,
#                                    `4` = 90, `5` = 100)),
#          rollV = runif(1, 0, 1),
#          visitNew = rollV < probVisit/100) %>%
#   left_join(probObs, by = 'stateNew') %>%
#   mutate(rollO = runif(1, 0, 1),
#          stateO = 1 + (rollO > probO1) + (rollO > probO1 + probO2)  +
#            (rollO > probO1 + probO2 + probO3) +
#            (rollO > probO1 + probO2 + probO3 + probO4)) %>%
#   ungroup %>% dplyr::select(
#     -c(prob1:prob5, probO1:probO5, rollT, rollV, rollO))
# }
# end2 <- Sys.time()
# 
# 
# # cohort update ----
# sample_row <- function(prob_row){ sample(1:5, size = 1, prob = prob_row)}
# 
# for (i in 1:nDays){
#   cohort <- cohort %>% left_join(dfTrans, by = 'state0') %>%
#     rowwise() %>%
#     mutate(stateNew = sample_row(c_across(prob1:prob5)),
#            day = i,
#            age = age + 1/365,
#            stay = ifelse(stateNew == state0, stay + 1, 0),
#            probVisit = ifelse(visit == 0,
#                               recode(stateNew, `1` = 1, `2` = 40, `3` = 30, 
#                                      `4` = 50, `5` = 90),
#                               recode(stateNew, `1` = 10, `2` = 80, `3` = 60, 
#                                      `4` = 90, `5` = 100)),
#            probVisit = probVisit/100,
#            visitNew = rbinom(1, size = 1, prob = probVisit)) %>%
#     left_join(dfObs, by = 'stateNew') %>%
#     mutate(stateO = sample_row(c_across(probO1:probO5))) %>%
#     ungroup %>% dplyr::select(
#       -c(prob1:prob5, probO1:probO5))
#   write.csv(cohort, paste0('store/cohortDay', str_pad(i, 7, pad = '0'), '.csv'),
#             row.names = F)
#   cohort %<>% mutate(state0 = stateNew, visit = visitNew)
#   }
# 
# 
# 

# generate observation data from simulated data ----
rm(list = ls())
# setwd('U:/LS7')
# setwd('D:/ASPIRE/')
setwd('demonstration/generalData/')


# Library ----
library(tidyverse)
library(magrittr)



# Data loading ----
days <- dir('store', full.names = T)
days <- days[str_detect(days, '.csv')]
days <- days[str_detect(days, '[0-9]')]
days <- days[-1]

hkd <- data.frame()
for (file in days){
  temp <- read.csv(file)
#  temp <- temp %>% filter(visitNew)
  hkd <- rbind(hkd, temp)
}

print(getwd())
print(dir())

write.csv(hkd, 'store/hkd.csv', row.names = F)




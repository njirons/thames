library(MASS)
library(tidyverse)
library(future)
library(furrr)
library(mvtnorm)
source("functions/thames_gauss_multi")

no_cores <- availableCores() - 1
plan(multisession, workers =8)
Param <- expand_grid(simu =1:50, d=1,t =10000,s0=1, n = 20)

Res_d1 <- future_pmap(Param,
                      simu_gauss_mutli , .progress=TRUE
)

save(Res_d1, file="Res_d1.Rdata")
Res1 <- do.call(rbind, Res_d1)
summary(Res1[,1]+Res1[,2])
## on va le faire d par d

Param <- expand_grid(simu =1:50, d=20,t =10000,s0=1, n = 20)
Res_d20 <- future_pmap(Param,
                       simu_gauss_mutli , .progress=TRUE
)

save(Res_d20, file="Res_d20.Rdata")

## on va le faire d par d

Param <- expand_grid(simu =1:50, d=50,t =10000,s0=1, n = 20)

Res_d50 <- future_pmap(Param,
                       simu_gauss_mutli , .progress=TRUE
)

save(Res_d50, file="Res_d50.Rdata")


# d 100

Param <- expand_grid(simu =1:50, d=100,t =10000,s0=1, n = 20)
Res_d100 <- future_pmap(Param,
                        simu_gauss_mutli , .progress=TRUE
)

save(Res_d100, file="Res_d100.Rdata")
Res100 <- do.call(rbind, Res_d100)
summary(Res100[,1]+Res100[,2])

load("Res_d1.Rdata")
load("Res_d20.Rdata")
load("Res_d50.Rdata")
load("Res_d100.Rdata")

Res1 <- do.call(rbind, Res_d1) %>%
  as.data.frame() %>%
  mutate(d=1)
Res20 <- do.call(rbind, Res_d20) %>%
  as.data.frame() %>%
  mutate(d=20)
Res50 <- do.call(rbind, Res_d50) %>%
  as.data.frame() %>%
  mutate(d=50)
Res100 <- do.call(rbind, Res_d100) %>%
  as.data.frame() %>%
  mutate(d=100)
Res1


Res_tot <- bind_rows(Res1,Res20,Res50,Res100) %>%
  mutate(diff = V1+V2)
colnames(Res_tot)=NULL
save(Res_tot, file ="Res_changemu.Rdata")
write.matrix(Res_tot,'data/gauss_multi.csv')


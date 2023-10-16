library(genridge)
library(MASS)
library(mvtnorm)
library(expm)
library(tidyverse)
library(SimDesign)
library(bridgesampling)
library(invgamma)
library(thames)

source("functions/thames_zellner.R")

set.seed(2023)

data <- prostate # prostate dataset from Stamey et al. (1989)
x <- data[, 1:8]
y <- data[, 9]

X = cbind(1,x)# select variables and add (artificially add an intercept)
Xmat = as.matrix(X)

#hyperparameters
n = dim(X)[1]
d = dim(X)[2]+1
s0_2=1
g=sqrt(n)
nu0=4

#Zellner-g-prior on mu, Negative Gamma prior on sigma
bias_max=0
res = NULL
#repeat once for every model
for(i in (3:9)){
  reg = regression(y, Xmat[,1:i], s0_2, g, nu0, split_sample=1)
  real = reg[[2]]
  bias_max = max(reg[[3]],bias_max)
  res_0 = cbind(rep(i-1,dim(reg[[1]])[1]),reg[[1]],rep(real,dim(reg[[1]])[1]))
  if(is.null(res)){
    res = res_0
  } else{
    res = rbind(res,res_0)
  }
}
bias_max  # numerical 0
colnames(res)=NULL
write.matrix(res,'data/zellner.csv')

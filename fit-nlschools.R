rm (list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(invgamma)
library (MASS)
attach (nlschools)
library (lme4)
library (arm)
library(bridgesampling)
library(mvtnorm)
library(tidyr)
library(dplyr)
library(thames)

set.seed(2023)

######################################################### prepare nlschools data

n <- length (lang)
nj <- tapply (lang,class,length)
nj_cumsum <- cumsum(c(0,nj))
classmeans <- tapply (lang,class,mean)
nclass <- length (nj)
ses <- SES - mean(SES)

# Varying intercept model for Netherlands schools

# prior hyperparameters
mu_prec <- 1/(2*var(lang))	#Prior sd of mu is 1.4 sd(lang)
beta_prec <- var(ses) / (2*var(lang)) #Prior sd of beta is 1.4 sd(lang)/sd(ses)
tau_rate <- 0.5*var(lang) 	#Prior mean of tau is 1/var(lang)
tau_alpha_rate <- 0.5*var(classmeans) #Prior mean(tau.alpha) = 1/var(classmeans)

mu_hat <- mean(lang)
sigma_mu_hat <- sqrt(2) * sd(lang)
nu_epsilon_hat <- 0.5
beta_epsilon_hat <- 0.5*var(lang)
nu_alpha_hat <- 0.5
beta_alpha_hat <- 0.5*var(classmeans)

#################################################### source functions and models

source('functions/thames_lm.R')
source('functions/thames_lmm_full.R')
source('functions/thames_lmm_reduced.R')
source('functions/nlschools_stan.R')

# compile stan models
lmm_full_model <- stan_model(model_code = lmm_full_stan)
lmm_reduced_model <- stan_model(model_code = lmm_reduced_stan)
lm_model <- stan_model(model_code = lm_stan)
lm_analytical <- stan_model(model_code = lm_stan_analytical)
lm_analytic0 <- stan_model(model_code = lm_stan_analytic0)

# define data as named list for stan
standata <- list(
  n = n,
  nclass = nclass,
  nj = as.integer(nj),
  nj_cumsum = as.integer(nj_cumsum),
  classmeans = as.numeric(classmeans),
  ses = ses,
  lang = lang,
  mu_prec = mu_prec,
  beta_prec = beta_prec,
  tau_rate = tau_rate,
  tau_alpha_rate = tau_alpha_rate,
  mu_hat = mu_hat,
  sigma_mu_hat = sigma_mu_hat,
  nu_epsilon_hat = nu_epsilon_hat,
  beta_epsilon_hat = beta_epsilon_hat,
  nu_alpha_hat = nu_alpha_hat,
  beta_alpha_hat = beta_alpha_hat
)

################################################################## fit the model

iters <- round(seq(5,5000,length.out=20)); # posterior sample sizes

for(iter in iters){
  lmm_full_fit <- sampling(lmm_full_model,
                           data = standata,
                           iter = 2*iter,
                           pars = c('mu','alpha',
                                    'sigma_alpha_squared','sigma_epsilon_squared',
                                    'alpha_long'))
  save(lmm_full_fit,file=paste0('data/lmm_full_fit_',4*iter,'.Rda'))
  rm(lmm_full_fit)
}

for(iter in iters){
  lmm_reduced_fit <- sampling(lmm_reduced_model,
                              data = standata,
                              iter = 2*iter,
                              pars = c('mu','sigma_alpha_squared','sigma_epsilon_squared'))
  save(lmm_reduced_fit,file=paste0('data/lmm_reduced_fit_',4*iter,'.Rda'))
  rm(lmm_reduced_fit)
}

for(iter in iters){
  lm_fit <- sampling(lm_model,
                     data = standata,
                     iter = 2*iter,
                     pars = c('mu','sigma_epsilon_squared'))
  save(lm_fit,file=paste0('data/lm_fit_',4*iter,'.Rda'))
  rm(lm_fit)
}

############################ compute log marginal likelihood via bridge sampling

# full lmm
lml_bridge_full <- matrix(NA, nrow = length(iters), ncol = 3)
new_mod <- stan(model_code = lmm_full_stan, data = standata, chains = 0)
for(j in 1:length(iters)){
  load(paste0('lmm_full_fit_',4*iters[j],'.Rda'))
  lml_bridge <- bridge_sampler(lmm_full_fit,new_mod,silent = TRUE)
  lml_err <- error_measures(lml_bridge)

  # estimate
  lml_bridge_full[j,1] <- lml_bridge$logml

  # lower 95% bound
  lml_bridge_full[j,2] <- lml_bridge$logml + log(trunc_quantile(0.025,lml_err$cv))

  # upper 95% bound
  lml_bridge_full[j,3] <- lml_bridge$logml + log(trunc_quantile(0.975,lml_err$cv))

  rm(lmm_full_fit)
  rm(lml_bridge)
  rm(lml_err)
}
write.matrix(lml_bridge_full,'data/lml_bridge_full.csv')

# reduced lmm
lml_bridge_reduced <- matrix(NA, nrow = length(iters), ncol = 3)
new_mod <- stan(model_code = lmm_reduced_stan, data = standata, chains = 0)
for(j in 1:length(iters)){
  load(paste0('lmm_reduced_fit_',4*iters[j],'.Rda'))
  lml_bridge <- bridge_sampler(lmm_reduced_fit,new_mod,silent = TRUE)
  lml_err <- error_measures(lml_bridge)

  # estimate
  lml_bridge_reduced[j,1] <- lml_bridge$logml

  # lower 95% bound
  lml_bridge_reduced[j,2] <- lml_bridge$logml + log(trunc_quantile(0.025,lml_err$cv))

  # upper 95% bound
  lml_bridge_reduced[j,3] <- lml_bridge$logml + log(trunc_quantile(0.975,lml_err$cv))

  rm(lmm_reduced_fit)
  rm(lml_bridge)
  rm(lml_err)
}
write.matrix(lml_bridge_reduced,'data/lml_bridge_reduced.csv')

# lm
lml_bridge_lm <- matrix(NA, nrow = length(iters), ncol = 3)
new_mod <- stan(model_code = lm_stan, data = standata, chains = 0)
for(j in 1:length(iters)){
  load(paste0('lm_fit_',4*iters[j],'.Rda'))
  lml_bridge <- bridge_sampler(lm_fit,new_mod,silent = TRUE)
  lml_err <- error_measures(lml_bridge)

  # estimate
  lml_bridge_lm[j,1] <- lml_bridge$logml

  # lower 95% bound
  lml_bridge_lm[j,2] <- lml_bridge$logml + log(trunc_quantile(0.025,lml_err$cv))

  # upper 95% bound
  lml_bridge_lm[j,3] <- lml_bridge$logml + log(trunc_quantile(0.975,lml_err$cv))

  rm(lm_fit)
  rm(lml_bridge)
  rm(lml_err)
}
write.matrix(lml_bridge_lm,'data/lml_bridge_lm.csv')

##################################### compute log marginal likelihood via THAMES

# full lmm
lml_thames_full <- matrix(NA, nrow = length(iters), ncol = 3)
for(k in 1:length(iters)){
  load(paste0('lmm_full_fit_',4*iters[k],'.Rda'))
  sims <- extract(lmm_full_fit)
  params <- cbind(sims$mu,sims$alpha,sims$sigma_alpha_squared,sims$sigma_epsilon_squared)
  lb <- c(rep(-Inf,nclass+1),0,0)
  ub <- rep(Inf,nclass+3)
  bound <- function(x){prod(x >= lb)*prod(x <= ub)}

  n_samples = dim(params)[1]
  d_par = dim(params)[2]
  c_opt = sqrt(d_par+1)

  lps <- rep(0,n_samples)
  for(j in 1:n_samples){
    lps[j] <- logpost_lmm_full(lang,sims$mu[j],mu_hat,sigma_mu_hat,
                               sims$alpha[j,],sims$alpha_long[j,],
                               sqrt(sims$sigma_alpha_squared[j]),
                               nu_alpha_hat,beta_alpha_hat,
                               sqrt(sims$sigma_epsilon_squared[j]),
                               nu_epsilon_hat,beta_epsilon_hat)
  }

  # estimate THAMES
  lml_thames <- thames(lps=lps,params=params,bound=bound)

  lml_thames_full[k,1] <- -lml_thames$log_zhat_inv

  # lower 95% bound
  lml_thames_full[k,2] <- -lml_thames$log_zhat_inv_U

  # upper 95% bound
  lml_thames_full[k,3] <- -lml_thames$log_zhat_inv_L

  rm(lmm_full_fit)
  rm(lml_thames)
  rm(sims)
  rm(params)
  rm(lps)
}
write.matrix(lml_thames_full,'data/lml_thames_full.csv')

# reduced lmm
lml_thames_reduced <- matrix(NA, nrow = length(iters), ncol = 3)
for(k in 1:length(iters)){
  load(paste0('lmm_reduced_fit_',4*iters[k],'.Rda'))
  sims <- extract(lmm_reduced_fit)
  params <- cbind(sims$mu,sims$sigma_alpha_squared,sims$sigma_epsilon_squared)
  lb <- c(-Inf,0,0)
  ub <- rep(Inf,3)
  bound <- function(x){prod(x >= lb)*prod(x <= ub)}

  n_samples = dim(params)[1]
  d_par = dim(params)[2]
  c_opt = sqrt(d_par+1)

  lps <- rep(0,n_samples)
  for(j in 1:n_samples){
    lps[j] <- logpost_lmm_reduced(lang,nj,nj_cumsum,nclass,
                                  sims$mu[j],mu_hat,sigma_mu_hat,
                                  sqrt(sims$sigma_alpha_squared[j]),
                                  nu_alpha_hat,beta_alpha_hat,
                                  sqrt(sims$sigma_epsilon_squared[j]),
                                  nu_epsilon_hat,beta_epsilon_hat)
  }

  # estimate THAMES
  lml_thames <- thames(lps=lps,params=params,bound=bound)

  lml_thames_reduced[k,1] <- -lml_thames$log_zhat_inv

  # lower 95% bound
  lml_thames_reduced[k,2] <- -lml_thames$log_zhat_inv_U

  # upper 95% bound
  lml_thames_reduced[k,3] <- -lml_thames$log_zhat_inv_L

  rm(lmm_reduced_fit)
  rm(lml_thames)
  rm(sims)
  rm(params)
  rm(lps)
}
write.matrix(lml_thames_reduced,'data/lml_thames_reduced.csv')

# lm
lml_thames_lm <- matrix(NA, nrow = length(iters), ncol = 3)
for(k in 1:length(iters)){
  load(paste0('lm_fit_',4*iters[k],'.Rda'))
  sims <- extract(lm_fit)
  params <- cbind(sims$mu,sims$sigma_epsilon_squared)
  lb <- c(-Inf,0)
  ub <- rep(Inf,2)
  bound <- function(x){prod(x >= lb)*prod(x <= ub)}

  n_samples = dim(params)[1]
  d_par = dim(params)[2]
  c_opt = sqrt(d_par+1)

  lps <- rep(0,n_samples)
  for(j in 1:n_samples){
    lps[j] <- logpost_lm(lang,sims$mu[j],mu_hat,sigma_mu_hat,
                         sqrt(sims$sigma_epsilon_squared[j]),
                         nu_epsilon_hat,beta_epsilon_hat)
  }

  # estimate THAMES
  lml_thames <- thames(lps=lps,params=params,bound=bound)

  lml_thames_lm[k,1] <- -lml_thames$log_zhat_inv

  # lower 95% bound
  lml_thames_lm[k,2] <- -lml_thames$log_zhat_inv_U

  # upper 95% bound
  lml_thames_lm[k,3] <- -lml_thames$log_zhat_inv_L

  rm(lm_fit)
  rm(lml_thames)
  rm(sims)
  rm(params)
  rm(lps)
}
write.matrix(lml_thames_lm,'data/lml_thames_lm.csv')

##################################### compute log marginal likelihood via simple Monte Carlo

# linear model
lml_mc_lm <- matrix(NA, nrow = length(iters), ncol = 3)
for(j in 1:length(iters)){
  # sample from prior
  mus <- rnorm(iters[j],mu_hat,sigma_mu_hat)
  sigma_epsilon <- sqrt(invgamma::rinvgamma(iters[j],shape=nu_epsilon_hat,rate=beta_epsilon_hat))
  samps <- cbind(mus,sigma_epsilon)

  # evaluate likelihood
  lls <- apply(samps,1,function(x){loglik_lm(lang,x[1],x[2])})

  # monte carlo approximation of log marginal likelihood
  lml <- log(mean(exp(lls - max(lls)))) + max(lls)
  lml_mc_lm[j,1] <- lml

  # calculate coefficient of variation
  cv <- sd(exp(lls-max(lls)))/
    (sqrt(iters[j])*mean(exp(lls-max(lls))))

  # calculate 95% lower bound
  lml_mc_lm[j,2] <- lml + log(trunc_quantile(0.025,cv))

  # calculate 95% upper bound
  lml_mc_lm[j,3] <- lml + log(trunc_quantile(0.975,cv))
}
write.matrix(lml_mc_lm,'data/lml_mc_lm.csv')

# LMM reduced
lml_mc_reduced <- matrix(NA, nrow = length(iters), ncol = 3)
for(j in 1:length(iters)){
  # sample from prior
  mus <- rnorm(iters[j],mu_hat,sigma_mu_hat)
  sigma_alpha <- sqrt(invgamma::rinvgamma(iters[j],shape=nu_alpha_hat,rate=beta_alpha_hat))
  sigma_epsilon <- sqrt(invgamma::rinvgamma(iters[j],shape=nu_epsilon_hat,rate=beta_epsilon_hat))
  samps <- cbind(mus,sigma_alpha,sigma_epsilon)

  # evaluate likelihood
  lls <- apply(samps,1,function(x){loglik_lmm_reduced(lang,nj,nj_cumsum,nclass,x[1],x[2],x[3])})

  # monte carlo approximation of log marginal likelihood
  lml <- log(mean(exp(lls - max(lls)))) + max(lls)
  lml_mc_reduced[j,1] <- lml

  # calculate coefficient of variation
  cv <- sd(exp(lls-max(lls)))/
    (sqrt(iters[j])*mean(exp(lls-max(lls))))

  # calculate 95% lower bound
  lml_mc_reduced[j,2] <- lml + log(trunc_quantile(0.025,cv))

  # calculate 95% upper bound
  lml_mc_reduced[j,3] <- lml + log(trunc_quantile(0.975,cv))
}
write.matrix(lml_mc_reduced,'data/lml_mc_reduced.csv')


# LMM full
lml_mc_full <- matrix(NA, nrow = length(iters), ncol = 3)
for(j in 1:length(iters)){
  # sample from prior
  mus <- rnorm(iters[j],mu_hat,sigma_mu_hat)
  sigma_alpha <- sqrt(invgamma::rinvgamma(iters[j],shape=nu_alpha_hat,rate=beta_alpha_hat))
  sigma_epsilon <- sqrt(invgamma::rinvgamma(iters[j],shape=nu_epsilon_hat,rate=beta_epsilon_hat))
  samps <- cbind(mus,sigma_alpha,sigma_epsilon)
  alphas <- t(apply(samps,1,function(x){rnorm(nclass,x[1],x[2])}))
  samps <- cbind(samps,alphas)

  # evaluate likelihood
  lls <- apply(samps,1,function(x){loglik_lmm_full(lang,tail(x,nclass),x[3])})

  # monte carlo approximation of log marginal likelihood
  lml <- log(mean(exp(lls - max(lls)))) + max(lls)
  lml_mc_full[j,1] <- lml

  # calculate coefficient of variation
  cv <- sd(exp(lls-max(lls)))/
    (sqrt(iters[j])*mean(exp(lls-max(lls))))

  # calculate 95% lower bound
  lml_mc_full[j,2] <- lml + log(trunc_quantile(0.025,cv))

  # calculate 95% upper bound
  lml_mc_full[j,3] <- lml + log(trunc_quantile(0.975,cv))
}
write.matrix(lml_mc_full,'data/lml_mc_full.csv')


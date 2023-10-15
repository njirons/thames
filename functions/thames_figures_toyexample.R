library(tidyverse)
library(EnvStats)
library(pracma)

### Simulations for thames_open ###

# Estimator in the simple gaussian case (with T=1)
# See f.e. draft_02 from the dropbox reports
# So far the mean vector and n are required as input: 
# Instead one could require full_data as input. The user needs to check that
# full_data matches with mean and n
thames_open_simple <- function(est_mu, n, c_squ, sigma_squ, 
                               type="inv_likeli", 
                               preset_random_numbers=FALSE){
  #browser()
  d = length(est_mu)
  sigma_n_squ = 1/(1/sigma_squ+n)
  m_n = n*est_mu*sigma_n_squ
  # starting values from draft_02
  
  if (isFALSE(preset_random_numbers)){
    param = rnorm(d)
  } else {
    param = preset_random_numbers
  }
  #Need this if I want to compare the same run for different c_squ
  
  # if (isFALSE(isFALSE(full_data))){
  #   log_constant_01 = sum(apply(full_data, 1,function(x) sum((x-mean(x))^2)/2))
  #   log_constant = log(2*pi/n)/2-n*log(2*pi)/2- log_constant_01
  # }
  #Since the arithmetic mean is a sufficient statistic, 
  #the results always only differ by a multiplicative constant w.r.t. theta
  
  if(type=="harmonic_mean"){
    param_sr <- param * sqrt(sigma_n_squ) + est_mu
    sum_likelihood = sum(n*(est_mu-param_sr)^2)
    return(exp(log(2*pi/n)*(d/2)+sum_likelihood/2))
  } 
  # else if(type=="harmonic_mean_full_data"){
  #   return(exp(log(2*pi/n)*(d/2)+sum_likelihood/2 + log_constant))
  # }
  
  if(sum(param^2)>c_squ){
    return(0)
  }
  # Same as in thames_open
  
  param_sr <- param * sqrt(sigma_n_squ) + est_mu
  sum_likelihood = sum(param_sr^2/sigma_squ + n*(est_mu-param_sr)^2)
  log_V_A = (d/2) * log(pi) + (d/2)*(log(c_squ)+log(sigma_n_squ))-lgamma(d/2+1)
  log_sol = sum_likelihood/2-(log_V_A + d*log(sqrt(n/sigma_squ)/(2*pi)))
  # Easy Gaussian case: A lot of stuff simplifies
  
  if(type=="inv_likeli"){   # usual THAMES
    return(exp(log_sol))
  } else if(type=="inv_likeli_scaled"){   # scaled THAMES
    marg_var = sigma_squ+1/n
    log_true_sol = (d/2)*log(2*pi*marg_var)+sum(est_mu^2)/(2*marg_var)
    return(exp(log_sol - log_true_sol))
  } 
  # else if(type=="inv_likeli_full_data"){
  #   return(exp(log_sol + log_constant))
  # }
    
} 

# the thames when T>1
mean_thames_open_simple = function(est_mu, n, c_squ, sigma_squ,T=10000,type="inv_likeli"){
  mean(sapply(1:T, function(t) thames_open_simple(est_mu, n, c_squ, sigma_squ,type=type)))
}

# mse estimation using a sample of size T
mse_thames_open_simple = function(T, est_mu, n, c_squ, sigma_squ, 
                                              type="inv_likeli_scaled"){
  if(type=="inv_likeli_scaled"){
    estimators = sapply(1:T, function(t) 
      thames_open_simple(est_mu, n, c_squ, sigma_squ, type))
    mean((estimators - 1)^2)  # due to unbiasedness
  } else if(type=="log_mean_likeli"){
    estimators = sapply(1:T, function(t) 
      mean_thames_open_simple(est_mu, n, c_squ, sigma_squ))
    d = length(est_mu)
    marg_var = sigma_squ+1/n
    log_true_sol = (d/2)*log(2*pi*marg_var)+sum(est_mu^2)/(2*marg_var)
    mean((log(estimators)-log_true_sol)^2)
  }
}

# comparison of mse with different radius functions c_d^2
sim_mse_wrt_d = function(Td, T, est_mu, n, c_squ_d, sigma_squ, 
                                type="inv_likeli_scaled"){
  sapply(1:Td, 
         function(d) mse_thames_open_simple(T, est_mu[1:d], n, 
                                            c_squ_d(d), sigma_squ, type))

}

# arithmetic mean with varying t in 1:T
sim_mean_thames_open_simple = function(T, est_mu, n, c_squ, sigma_squ, 
                                                   type="inv_likeli", 
                                                   preset_random_numbers = FALSE, full_data=FALSE){
  estimators = sapply(1:T, function(t) 
    thames_open_simple(est_mu, n, c_squ, sigma_squ, 
                       type, preset_random_numbers[t,]))
  if(isFALSE(isFALSE(full_data))){
    log_constant_01 = sum(apply(full_data, 1,function(x) sum((x-mean(x))^2)/2))
    log_constant = log(2*pi/n)/2-n*log(2*pi)/2- log_constant_01
    cumsum(estimators)*exp(-log_constant)/(1:T) 
    #See rapport_14 "Changing the marginal of the sufficient statistic"
  } else{
    cumsum(estimators)/(1:T)
  }
}


loglik_lmm_reduced <- function(y,nj,nj_cumsum,nclass,mu,sigma_alpha,sigma_epsilon){
  ll <- 0
  for(j in 1:nclass){
    ll <- ll + dmvnorm(lang[(nj_cumsum[j]+1):nj_cumsum[j+1]], 
                       mean = rep(mu,nj[j]), 
                       sigma = diag(sigma_epsilon^2,nj[j]) + 
                         (sigma_alpha^2)*matrix(data=1,nrow=nj[j],ncol=nj[j]), 
                       log = TRUE)
  }
  return(ll)
}

logprior_lmm_reduced <- function(mu,mu_hat,sigma_mu_hat,
                                 sigma_alpha,nu_alpha_hat,beta_alpha_hat,
                                 sigma_epsilon,nu_epsilon_hat,beta_epsilon_hat){
  lp <- dnorm(mu,mu_hat,sigma_mu_hat,log=TRUE) + 
    invgamma::dinvgamma(
      sigma_alpha^2,
      shape = nu_alpha_hat, 
      rate = beta_alpha_hat,
      log=TRUE) + 
    invgamma::dinvgamma(
      sigma_epsilon^2,
      shape = nu_epsilon_hat, 
      rate = beta_epsilon_hat,
      log=TRUE)
  return(lp)
}

logpost_lmm_reduced <- function(y,nj,nj_cumsum,nclass,
                                mu,mu_hat,sigma_mu_hat,
                                sigma_alpha,nu_alpha_hat,beta_alpha_hat,
                                sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat){
  ll <- loglik_lmm_reduced(y,nj,nj_cumsum,nclass,mu,sigma_alpha,sigma_epsilon)
  lp <- logprior_lmm_reduced(mu,mu_hat,sigma_mu_hat,
                             sigma_alpha,nu_alpha_hat,beta_alpha_hat,
                             sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat)
  return(ll+lp)
}


# calculate likelihood L(y|mu,beta,tau,tau_alpha) after integrating out random intercepts
# using prior alpha_i ~ N(mu, tau^2) (iid)
nlschools_ll <- function(mu,beta,tau,tau_alpha){
  ll <- 0
  for(j in 1:nclass){
    ll <- ll + dmvnorm(lang[(nj_cumsum[j]+1):nj_cumsum[j+1]], 
                       mean = mu+beta*ses[(nj_cumsum[j]+1):nj_cumsum[j+1]], 
                       sigma = diag(1/tau,nj[j]) + 
                         (1/tau_alpha)*matrix(data=1,nrow=nj[j],ncol=nj[j]), 
                       log = TRUE)
  }
  return(ll)
}

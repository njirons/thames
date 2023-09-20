loglik_lmm_full <- function(y,alpha_long,sigma_epsilon){
  ll <- sum(dnorm(
    x = y, 
    mean = alpha_long,
    sd = sigma_epsilon,
    log = TRUE
  ))
  return(ll)
}

logprior_lmm_full <- function(
    mu,mu_hat,sigma_mu_hat,
    alpha,sigma_alpha,nu_alpha_hat,beta_alpha_hat,
    sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat){
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
      log=TRUE) +
    sum(dnorm(alpha,mu,sigma_alpha,log=TRUE))
  
  return(lp)
}

logpost_lmm_full <- function(y,mu,mu_hat,sigma_mu_hat,
                             alpha,alpha_long,sigma_alpha,nu_alpha_hat,beta_alpha_hat,
                             sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat){
  ll <- loglik_lmm_full(y,alpha_long,sigma_epsilon)
  lp <- logprior_lmm_full(mu,mu_hat,sigma_mu_hat,
                          alpha,sigma_alpha,nu_alpha_hat,beta_alpha_hat,
                          sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat)
  return(ll+lp)
}

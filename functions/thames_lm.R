loglik_lm <- function(y,mu,sigma_epsilon){
  return(sum(dnorm(y,mu,sigma_epsilon,log=TRUE)))
}

logprior_lm <- function(mu,mu_hat,sigma_mu_hat,
                        sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat){
  lp <- dnorm(mu,mu_hat,sigma_mu_hat,log=TRUE) + 
    dinvgamma(
      sigma_epsilon^2,
      shape = nu_epsilon_hat, 
      rate = beta_epsilon_hat,
      log=TRUE)
  return(lp)
}

logpost_lm <- function(y,mu,mu_hat,sigma_mu_hat,
                       sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat){
  ll <- loglik_lm(y,mu,sigma_epsilon)
  lp <- logprior_lm(mu,mu_hat,sigma_mu_hat,
                    sigma_epsilon, nu_epsilon_hat,beta_epsilon_hat)
  return(ll+lp)
}

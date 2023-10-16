reg_log_likelihood <- function(param, X, sig2,...) {
  n <- nrow(X)
  d = length(param)
  sum(sapply(1:n,function(i){mvtnorm::dmvnorm(X[i,], param, diag(d), 
                                              log = TRUE)}))
} 

reg_log_prior <- function(param, sig2=1,...) {
  d <- length(param)
  p = dmvnorm(param, rep(0,d), (sig2)*diag(d), log = TRUE)
  return(p)
}

simu_gauss_mutli <- function(n, t,d, s0 =1, simu){
  set.seed(3785)
  mu_star = rnorm(d,0,1)
  
  set.seed(simu*18+36)
  Y = rmvnorm(n, mu_star, diag(rep(1,d)))
  
  sn = 1/(n+1/s0)
  mn = apply(Y,2,sum)/(n + 1/s0)
  
  true_marginal = sum(unlist(apply(Y,2,function(i){
    mvtnorm::dmvnorm(i, rep(0, n),  s0 *matrix(1, n, n)+diag(n), log = TRUE)})))
  
  mu_sample <- rmvnorm(t, mean = mn,sn*diag(d))
  log_prior <- apply(mu_sample, 1, reg_log_prior, s0)
  log_likeliood <- apply(mu_sample, 1, reg_log_likelihood, Y)
  ap <- thames(log_prior + log_likeliood, as.matrix(mu_sample))
  return(c(ap$log_zhat_inv, true_marginal))
}
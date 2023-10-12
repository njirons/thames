log_reg_likelihood <- function(param, Xmat, Y) {
  if(!is.matrix(param)){
    param=matrix(param,nrow=1)
  }
  n = dim(Xmat)[1]
  p = numeric(dim(param)[1])
  for(i in seq_along(p)){
    beta = param[i,1:dim(Xmat)[2]]
    sig2 = param[i,dim(Xmat)[2]+1]
    p[i] = mvtnorm::dmvnorm(Y, Xmat%*%beta, sig2*diag(n),log=TRUE)
  }
  return(p)
}

log_reg_prior <- function(param, Xmat, s0_2, g, nu0){
  if(!is.matrix(param)){
    param=matrix(param,nrow=1)
  }
  d = dim(param)[2]

  p_sigma2 = dinvgamma(param[,d],nu0/2,nu0*s0_2/2,log=TRUE)
  p_beta = numeric(dim(param)[1])
  for (i in seq_along(param[,1])){
    p_beta[i] = mvtnorm::dmvnorm(t(param[i,1:(d-1)]),rep(0,d-1),param[i,d]*solve(t(Xmat)%*%Xmat)*g,log=TRUE)
  }
  p = p_sigma2 + p_beta
  return(p)
}

log_reg_posterior <- function(param, Xmat, y, s0_2, g, nu0){
  if(!is.matrix(param)){
    param=matrix(param,nrow=1)
  }
  d = dim(param)[2]

  m_n = solve(t(Xmat)%*%Xmat) %*% t(Xmat) %*% y
  s_n = t(y)%*%y - (g/(g+1))*t(y)%*%Xmat%*%m_n
  p_sigma2 = dinvgamma(param[,d],(n+nu0)/2,(nu0*s0_2+s_n)/2,log = TRUE)
  p_beta = numeric(dim(param)[1])
  for (i in seq_along(param[,1])){
    p_beta[i] = dmvnorm(param[i,(1:(d-1))], mean=m_n*(g/(g+1)), sigma = param[i,d]*solve(t(Xmat)%*%Xmat)*(g/(g+1)),log=TRUE)
  }
  p = p_sigma2 + p_beta
  return(p)
}

regression = function(y, Xmat, s0_2, g, nu0, split_sample=0){
  n = dim(Xmat)[1] # number of observations
  d = dim(Xmat)[2] + 1 #sigma is also unknown
  #browser()
  nb_sample = 10000
  m_n = solve(t(Xmat)%*%Xmat) %*% t(Xmat) %*% y
  s_n = t(y)%*%y - (g/(g+1))*t(y)%*%Xmat %*% m_n
  sig2_sample = rinvgamma(nb_sample,(n+nu0)/2,(nu0*s0_2+s_n)/2)
  beta_sample = matrix(nrow=nb_sample,ncol=d-1)
  for(i in seq_along(beta_sample[,1])){
    beta_sample[i,] = rmvnorm(1, mean=m_n*(g/(g+1)), sigma = sig2_sample[i]*solve(t(Xmat)%*%Xmat)*(g/(g+1)))
  }
  param = cbind(beta_sample,t(t(sig2_sample)))
  log_unnorm_posterior = log_reg_prior(param, Xmat, s0_2, g, nu0) + log_reg_likelihood(param, Xmat, y)

  real = -((d-1)/2)*log(g+1)-(n/2)*log(pi)+lgamma((nu0+n)/2)-lgamma(nu0/2)+(nu0/2)*log(nu0*s0_2)-((nu0+n)/2)*log(nu0*s0_2+s_n)
  #real2=log_unnorm_posterior[10]-log_reg_posterior(param, Xmat, y, s0_2, g, nu0)[10]

  x<- c(50,100,200,500,seq(1000,10000,1000))

  ub = rep(Inf,d)
  lb = c(rep(-Inf,d-1),0)
  names(lb) <- names (ub) <- colnames(param) <- names(as.data.frame(param))
  bridge=bridge_sampler(samples=param, log_posterior=function(param.row,data) log_reg_prior(param.row,Xmat, s0_2, g, nu0) + log_reg_likelihood(param.row, Xmat, y),lb=lb,ub=ub,data=NULL,silent=TRUE)

  bridge_logml = sapply(x, function(nb_sample) bridge_sampler(samples=param[1:nb_sample,], log_posterior=function(param.row,data) log_reg_prior(param.row,Xmat, s0_2, g, nu0) + log_reg_likelihood(param.row, Xmat, y),
                                                              lb=lb,ub=ub,data=NULL,silent=TRUE)$logml)
  bridge_cv = sapply(x, function(nb_sample) error_measures(bridge_sampler(samples=param[1:nb_sample,], log_posterior=function(param.row,data) log_reg_prior(param.row,Xmat, s0_2, g, nu0) + log_reg_likelihood(param.row, Xmat, y),
                                                                          lb=lb,ub=ub,data=NULL,silent=TRUE))$cv)

  bridge_L = bridge_logml + log(trunc_quantile(0.025,bridge_cv))
  bridge_mid = bridge_logml
  bridge_U = bridge_logml + log(trunc_quantile(0.975,bridge_cv))

  tic_U = sapply(x, function(nb_sample) thames(lps=log_unnorm_posterior[1:nb_sample],params=param[1:nb_sample,])$log_zhat_inv_U)
  tic_mid = sapply(x, function(nb_sample) thames(lps=log_unnorm_posterior[1:nb_sample],params=param[1:nb_sample,])$log_zhat_inv)
  tic_L = sapply(x, function(nb_sample) thames(lps=log_unnorm_posterior[1:nb_sample],params=param[1:nb_sample,])$log_zhat_inv_L)

  naive_mc_sample = matrix(ncol=d,nrow=nb_sample)
  for(i in (1:nb_sample)){
    sig2_prior_sample = rinvgamma(1,nu0/2,s0_2*nu0/2)
    beta_prior_sample = rmvnorm(1, mean=rep(0,d-1),
                                sigma=g*sig2_prior_sample*solve(t(Xmat)%*%Xmat))
    naive_mc_sample[i,] = c(beta_prior_sample,sig2_prior_sample)
  }
  naive_mc_likelihood = apply(naive_mc_sample,1,
                              function(param) log_reg_likelihood(param, Xmat, y))

  naive_mc_point_est = log(cumsum(exp(naive_mc_likelihood - max(naive_mc_likelihood)))[x]/x) + max(naive_mc_likelihood)
  naive_mc_log_se = sapply(x, function(t) log(sd(exp(naive_mc_likelihood[1:t] - max(naive_mc_likelihood)))/sqrt(t)) + max(naive_mc_likelihood))
  naive_mc_cv = exp(naive_mc_log_se - naive_mc_point_est)

  naive_mc_U = naive_mc_point_est + log(trunc_quantile(0.025,naive_mc_cv))
  naive_mc_mid = naive_mc_point_est
  naive_mc_L = naive_mc_point_est + log(trunc_quantile(0.975,naive_mc_cv))

  bias =  thames(lps=log_unnorm_posterior[1:nb_sample],params=param[1:nb_sample,])
  bias_max = abs(log(1-(1/50)*(bias$cv^2)*exp(-bias$log_zhat_inv)))

  return(list(cbind(-tic_U, -tic_mid, -tic_L,bridge_U,bridge_mid,bridge_L,
                    naive_mc_U, naive_mc_mid, naive_mc_L), real, bias_max))
}

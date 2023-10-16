library(pracma)
library(gtools)
library(stats)
library(corpcor)
library(gridExtra)
library(thames)

d_samples = c(1,20,50,100)
y_sims = 50
res_unif_nobound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_unif_bound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_random_bound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_random_nobound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)

set.seed(1)
for(c in (1:2)){
  for(j in (1:length(d_samples))){
    #browser()
    d = d_samples[j]
    K = d+1
    alpha = rep(1,K)
    # Since the last parameter is entirely determined by the first K-1 parameters,
    # the number of dimensions is one lower than the vector
    n = 400
    l = 150
    if(c==1){
      prob = alpha/K#uniform weights
    } else {
      prob = rdirichlet(1,alpha) # 1 fixed probability for each d
    }
    for(i in (1:y_sims)){
      T=10000
      print(c*i*j)
      Y = rmultinom(n,size=l,prob=prob) #The data is simulated from the model

      log_prior = function(p) -sum(lgamma(alpha)) + lgamma(sum(alpha)) + sum(log(c(p,(1-sum(p))))*(alpha-1))

      log_likelihood = function(p) sum(apply(Y,2,function(y) dmultinom(y,prob=c(p,1-sum(p)),log=TRUE)))

      post_alpha = apply(Y,1,sum) + matrix(alpha,nrow=1)
      log_coeff_1 = lgamma(sum(alpha))-sum(lgamma(alpha))
      log_coeff_2 = n*lgamma(l+1) - sum(lgamma(Y+1))
      log_coeff_3 = sum(lgamma(post_alpha)) - lgamma(sum(post_alpha))
      log_marginal_likelihood = log_coeff_1 + log_coeff_2 + log_coeff_3

      param_full = rdirichlet(T,post_alpha)[,(1:d)]

      log_joint = function(param) -(log_prior(param)+log_likelihood(param))

      thames_nobound = thames(params=t(t(param_full)),
                              lp_func=function(param) -apply(param,1, log_joint))
      log_thames_nobound = thames_nobound$log_zhat_inv

      thames_bound = thames(params=t(t(param_full)),
                            lp_func=function(param) -apply(param,1, log_joint),
                            bound=function(x) (x>0) * (sum(x) <= 1))
      log_thames_bound = thames_bound$log_zhat_inv

      if(c==1){
        res_unif_nobound[i,j] = -log_thames_nobound - log_marginal_likelihood
        res_unif_bound[i,j] = -log_thames_bound - log_marginal_likelihood
      } else{
        res_random_nobound[i,j] = -log_thames_nobound - log_marginal_likelihood
        res_random_bound[i,j] = -log_thames_bound - log_marginal_likelihood
      }
    }
  }
}

res = cbind(res_unif_bound, res_unif_nobound, res_random_bound, res_random_nobound)
write.matrix(res,'data/multinomial.csv')

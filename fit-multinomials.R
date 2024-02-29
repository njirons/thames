source("multinom-logit_stan.R")

PACKAGES=c("pracma","gtools","stats","corpcor","gridExtra",
           "thames","bridgesampling","ohenery","MASS",
           "foreach","parallel", "doParallel","rstan")
lapply(PACKAGES,require, character.only = TRUE)

d_samples = c(1,20,50,100)
y_sims = 50
res_unif_nobound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_unif_bound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_random_bound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)
res_random_nobound = matrix(rep(0,y_sims*length(d_samples)),ncol=length(d_samples),nrow=y_sims)

res_thames_transformed = matrix(ncol=length(d_samples),nrow=y_sims)
res_bridge = matrix(ncol=length(d_samples),nrow=y_sims)
res_mcest = matrix(ncol=length(d_samples),nrow=y_sims)

res_systemtimes_bridge = matrix(ncol=length(d_samples),nrow=y_sims)
res_systemtimes_thames = matrix(ncol=length(d_samples),nrow=y_sims)
res_systemtimes_mcest = matrix(ncol=length(d_samples),nrow=y_sims)

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
    ncores=8

    cores=detectCores()
    cl <- makeCluster(min(cores[1]-1,ncores)) #not to overload your computer
    registerDoParallel(cl)
    #browser()
    this <- foreach(i=1:y_sims, .combine=cbind,
                    .packages=PACKAGES,
                    .export=c(names(as.list(.GlobalEnv)),ls())) %dopar% {
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

                      # Generating the sample from the transformed posterior and calculating
                      # THAMES, bridge_sampling and mc_est on it
                      if(c==1){

                        samples= cbind(1-rowSums(t(t(param_full))),t(t(param_full)))
                        transformed_samples = t(apply(samples,1,inv_smax))
                        transformed_samples = t(t(transformed_samples[,-1]))
                        if(d==1){
                          colnames(transformed_samples) = c("1")
                        } else{
                          colnames(transformed_samples) = c("1",2:dim(transformed_samples)[2])
                        }
                        ub = rep(Inf,d)
                        lb = rep(-Inf,d)
                        names(lb) <- names(ub) <- colnames(transformed_samples)

                        # Jacobian of the softmax
                        soft_max_der = function(par) diag(par)-
                          matrix(rep(par,length(par)),ncol=length(par),nrow=length(par))*t(matrix(rep(par,length(par)),ncol=length(par),nrow=length(par)))

                        # multivariate chain rule component
                        sum_der = function(par) rbind(diag(length(par)),-1)
                        #browser()
                        log_bridge = bridge_sampler(samples=transformed_samples,
                                                    log_posterior=function(param.row,data) -log_joint(smax(c(-sum(param.row),param.row))[-1])
                                                    +sum(log(abs(eigen(cbind(diag(length(param.row)),0)%*%soft_max_der(smax(c(-sum(param.row),param.row)))%*%
                                                                         sum_der(param.row))$values))),# chain rule component
                                                    data=NULL,ub=ub,lb=lb, silent=TRUE)$logml
                        log_thames_transformed = thames(params=transformed_samples,
                                                        lp_func=function(param) apply(param,1,function(param.row)
                                                          -log_joint(smax(c(-sum(param.row),param.row))[-1])
                                                          +sum(log(abs(eigen((cbind(diag(length(param.row)),0)%*%soft_max_der(smax(c(-sum(param.row),param.row)))%*%
                                                                                sum_der(param.row)))$values)))))$log_zhat_inv
                        prior_sample = rdirichlet(T,alpha)
                        log_likelihood_prior_sample = sapply(1:T,function(s) log_likelihood(prior_sample[s,-1]))

                        # define naive Monte Carlo estimator
                        mc_est = function(log_likelihood_prior_sample){
                          res = list(log_est = log(mean(exp(log_likelihood_prior_sample-max(log_likelihood_prior_sample))))+max(log_likelihood_prior_sample),
                                     log_sd=sd(exp(log_likelihood_prior_sample-max(log_likelihood_prior_sample))))
                          return(res)
                        }

                        log_mc_est = mc_est(log_likelihood_prior_sample)$log_est

                        # Precompute lps and calculate computation times on a stan fit object
                        lps = apply(transformed_samples ,1,function(param.row)
                          -log_joint(smax(c(-sum(param.row),param.row))[-1])
                          +sum(log(abs(eigen((soft_max_der(smax(c(-sum(param.row),param.row)))%*%
                                                sum_der(param.row))[-1,])$values))))

                        lmm_full_test <- stan_model(model_code = lmm_test)
                        standata = list(d=K,n_samp=n,y=t(Y),a=alpha)
                        lmm_full_fit <- sampling(lmm_full_test,
                                                 data = standata,
                                                 iter = T,
                                                 chains=1,
                                                 warmup=0,
                                                 algorithm="Fixed_param",
                                                 pars = c('theta_tilde'))

                        # Set samples and lps to exact simulation from posterior
                        #browser()
                        for(i in 1:d){
                          lmm_full_fit@sim$samples[[1]][[i]] = transformed_samples[,i]
                        }
                        #browser()
                        #lmm_full_fit@sim$samples[[1]]$`theta_tilde[1]`=transformed_samples
                        lmm_full_fit@sim$samples[[1]]$lp__=lps[lmm_full_fit@sim$permutation[[1]]]
                        new_mod <- stan(model_code = lmm_test, data = standata, chains = 0)

                        # Calculate system times
                        systemtimes_bridge = system.time(bridge_sampler(lmm_full_fit,new_mod,silent = TRUE)$logml)[3]
                        systemtimes_thames = system.time(thames(params=transformed_samples,lps=lps)$log_zhat_inv)[3]
                        systemtimes_mcest = system.time(mc_est(log_likelihood_prior_sample)$log_est)[3]
                      }

                      if(c==1){
                        c(-log_thames_nobound-log_marginal_likelihood,
                          -log_thames_bound-log_marginal_likelihood,
                          log_mc_est-log_marginal_likelihood,
                          log_bridge-log_marginal_likelihood,
                          -log_thames_transformed-log_marginal_likelihood,
                          systemtimes_bridge, systemtimes_thames, systemtimes_mcest)
                      } else{
                        c(-log_thames_nobound, -log_thames_bound)-log_marginal_likelihood
                      }
                    }
    stopCluster(cl)
    #browser()
    if(c==1){
      res_unif_nobound[,j] = this[1,]#res_thames_thames_nobound[,j] = this[1,]#-log_thames_transformed - log_marginal_likelihood
      res_unif_bound[,j] = this[2,]#res_bridge[,j] = this[2,]#log_bridge - log_marginal_likelihood
      res_mcest[,j] = this[3,]
      res_bridge[,j] = this[4,]#-log_thames_nobound - log_marginal_likelihood
      res_thames_transformed[,j] = this[5,]#-log_thames_bound - log_marginal_likelihood
      res_systemtimes_bridge[,j] = this[6,]
      res_systemtimes_thames[,j] = this[7,]
      res_systemtimes_mcest[,j] = this[8,]
    } else{
      res_random_nobound[,j] = this[1,]#-log_thames_nobound - log_marginal_likelihood
      res_random_bound[,j] = this[2,]#-log_thames_bound - log_marginal_likelihood
    }
    #browser()

  }
}

res = cbind(res_unif_nobound, res_unif_bound,
            res_random_nobound, res_random_bound,
            res_mcest, res_bridge, res_thames_transformed,
            res_systemtimes_bridge, res_systemtimes_thames, res_systemtimes_mcest)
write.matrix(res,'data/multinomial.csv')

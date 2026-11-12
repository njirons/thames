# full lmm
lmm_full_stan <- 'data {
  int<lower=1> n; // total number of students
  int<lower=1> nclass; // number of classes
  int<lower=0> nj[nclass]; // number of students in each class
  int<lower=0> nj_cumsum[nclass+1];
  vector[nclass] classmeans; // class means of lang
  vector[n] lang; // language scores
  
  real mu_hat;
  real<lower=0> sigma_mu_hat;
  real<lower=0> nu_epsilon_hat;
  real<lower=0> beta_epsilon_hat;
  real<lower=0> nu_alpha_hat;
  real<lower=0> beta_alpha_hat;
}
parameters {
  real mu;
  vector[nclass] alpha;
  real<lower=0> sigma_alpha_squared;
  real<lower=0> sigma_epsilon_squared;
}
transformed parameters {
  vector[n] alpha_long;
  for(j in 1:nclass){
    alpha_long[(nj_cumsum[j]+1):nj_cumsum[j+1]] = rep_vector(alpha[j],nj[j]);
  }
}
model {
  // log prior
  target += normal_lpdf(mu | mu_hat, sigma_mu_hat);
  target += inv_gamma_lpdf(sigma_alpha_squared | nu_alpha_hat, beta_alpha_hat);
  target += inv_gamma_lpdf(sigma_epsilon_squared | nu_epsilon_hat, beta_epsilon_hat);
  target += normal_lpdf(alpha | mu, sqrt(sigma_alpha_squared));
  
  // log likelihood
  target += normal_lpdf(lang | alpha_long, sqrt(sigma_epsilon_squared));
}
'

# reduced lmm (random effects integrated out)
lmm_reduced_stan <- 'data {
  int<lower=1> n; // total number of students
  int<lower=1> nclass; // number of classes
  int<lower=0> nj[nclass]; // number of students in each class
  int<lower=0> nj_cumsum[nclass+1];
  vector[nclass] classmeans; // class means of lang
  vector[n] lang; // language scores
  
  real mu_hat;
  real<lower=0> sigma_mu_hat;
  real<lower=0> nu_epsilon_hat;
  real<lower=0> beta_epsilon_hat;
  real<lower=0> nu_alpha_hat;
  real<lower=0> beta_alpha_hat;
}
parameters {
  real mu;
  real<lower=0> sigma_alpha_squared;
  real<lower=0> sigma_epsilon_squared;
}
model {
  // log prior
  target += normal_lpdf(mu | mu_hat, sigma_mu_hat);
  target += inv_gamma_lpdf(sigma_alpha_squared | nu_alpha_hat, beta_alpha_hat);
  target += inv_gamma_lpdf(sigma_epsilon_squared | nu_epsilon_hat, beta_epsilon_hat);
  
  // log likelihood
  for(j in 1:nclass){
    target += multi_normal_lpdf(lang[(nj_cumsum[j]+1):nj_cumsum[j+1]] | rep_vector(mu,nj[j]), diag_matrix(rep_vector(sigma_epsilon_squared,nj[j]))+sigma_alpha_squared*rep_matrix(1,nj[j],nj[j]));
  }
}
'

# simple mean model
lm_stan <- 'data {
  int<lower=1> n; // total number of students
  int<lower=1> nclass; // number of classes
  int<lower=0> nj[nclass]; // number of students in each class
  int<lower=0> nj_cumsum[nclass+1];
  vector[nclass] classmeans; // class means of lang
  vector[n] lang; // language scores
  
  real mu_hat;
  real<lower=0> sigma_mu_hat;
  real<lower=0> nu_epsilon_hat;
  real<lower=0> beta_epsilon_hat;
}
parameters {
  real mu;
  real<lower=0> sigma_epsilon_squared;
}
model {
  // log prior
  target += normal_lpdf(mu | mu_hat, sigma_mu_hat);
  target += inv_gamma_lpdf(sigma_epsilon_squared | nu_epsilon_hat, beta_epsilon_hat);
  
  // log likelihood
  //target += normal_lpdf(lang | rep_vector(mu,n), sqrt(sigma_epsilon_squared));
  target += normal_lpdf(lang | mu, sqrt(sigma_epsilon_squared));
}
'

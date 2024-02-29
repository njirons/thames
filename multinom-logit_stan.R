lmm_test <- 'data {
  int<lower=1> d;
  int<lower=1> n_samp;
  //int<lower=0> y[d];
  int y[n_samp,d];
  vector[d] a;
}
parameters {
  //simplex[d] theta;
  vector[d-1] theta_tilde;
}
transformed parameters {
  simplex[d] theta;
  theta = append_row(exp(-sum(theta_tilde)),exp(theta_tilde))/sum(append_row(exp(-sum(theta_tilde)),exp(theta_tilde)));
  //theta = append_row(exp(theta_tilde)./(1+exp(theta_tilde)),1-sum(exp(theta_tilde)./(1+exp(theta_tilde))));
}

model {

  // chain rule component
  matrix[d, d] softmaxder;

  for(i in 1:d){
    for(j in i:d){
      softmaxder[i,j] = -theta[i]*theta[j];
      if(i==j){
        softmaxder[i,i] += theta[i];
      }
      softmaxder[j,i] = softmaxder[i,j];
    }
  }

  matrix[d,d-1] sum_der;

  for(i in 1:d){
    for(j in 1:(d-1)){
      if(i<d){
        if(i==j){
          sum_der[i,j] = 1;
        } else{
          sum_der[i,j] = 0;
        }
      } else{
        sum_der[i,j] = -1;
      }
    }
  }

  matrix[d-1,d] restriction_der;

  for(i in 1:(d-1)){
    for(j in 1:d){
      if(i==j){
        restriction_der[i,j] = 1;
      } else{
        restriction_der[i,j] = 0;
      }
    }
  }

  matrix[d-1,d-1] res;
  res = (restriction_der*softmaxder)*sum_der;

  target += log_determinant(res);

  // priors
  target += dirichlet_lpdf(theta | a);

  for(j in 1:n_samp){
    target += multinomial_lpmf(y[j] | theta);
  }

}'

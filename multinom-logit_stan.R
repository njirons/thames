lmm_test <- '
data {
  int<lower=1> d;
  int<lower=1> n_samp;
  array[n_samp, d] int y;
  vector[d] a;
}

parameters {
  vector[d - 1] theta_tilde;
}

transformed parameters {
  simplex[d] theta;

  theta = append_row(exp(-sum(theta_tilde)), exp(theta_tilde))
          / sum(append_row(exp(-sum(theta_tilde)), exp(theta_tilde)));
}

model {
  matrix[d, d] softmaxder;
  matrix[d, d - 1] sum_der;
  matrix[d - 1, d] restriction_der;
  matrix[d - 1, d - 1] res;

  // Chain rule component
  for (i in 1:d) {
    for (j in i:d) {
      softmaxder[i, j] = -theta[i] * theta[j];
      if (i == j) {
        softmaxder[i, i] += theta[i];
      }
      softmaxder[j, i] = softmaxder[i, j];
    }
  }

  for (i in 1:d) {
    for (j in 1:(d - 1)) {
      if (i < d) {
        if (i == j) {
          sum_der[i, j] = 1;
        } else {
          sum_der[i, j] = 0;
        }
      } else {
        sum_der[i, j] = -1;
      }
    }
  }

  for (i in 1:(d - 1)) {
    for (j in 1:d) {
      if (i == j) {
        restriction_der[i, j] = 1;
      } else {
        restriction_der[i, j] = 0;
      }
    }
  }

  res = (restriction_der * softmaxder) * sum_der;
  target += log_determinant(res);

  // Prior
  target += dirichlet_lpdf(theta | a);

  // Likelihood
  for (j in 1:n_samp) {
    target += multinomial_lpmf(y[j] | theta);
  }
}
'

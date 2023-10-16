#' Bounded parameter correction
#'
#' Computes Monte Carlo estimate of the bounded parameter correction to THAMES
#'
#' @param theta_hat center of the ellipsoid A.
#' @param sigma_svd Singular value decomposition (SVD) of the posterior
#'     covariance matrix used to define A. If the dimension of the parameters space
#'     is 1, sigma_svd shall be equal to the variance.
#' @param bound Function calculating membership of a point in the posterior support.
#' @param radius positive number, radius determining the ellipsoid A
#' @param n_simuls Integer, number of Monte Carlo simulations to use in the calculation.
#'
#' @return A number in [0,1] estimating of the proportion of the volume of A
#'     contained in the posterior support.
#'
#' @keywords internal
bound_par_cor = function(theta_hat, sigma_svd, bound, radius, n_simuls=1e5){
  #browser()
  if(length(sigma_svd)==1){
    simuls = runif(n_simuls,min=theta_hat-radius*sqrt(sigma_svd),
                   max=theta_hat+radius*sqrt(sigma_svd))
    return(mean(Vectorize(bound)(simuls)))
  } else{
    #browser()
    # calculate inverse of sigma_hat
    sigma_inv <- sigma_svd$v %*% diag(sigma_svd$d^(-1)) %*% t(sigma_svd$v)
    
    # sample uniformly from the ellipsoid A
    simuls = runif_in_ellipsoid(n_simuls, sigma_inv, radius) +
      t(matrix(rep(theta_hat,n_simuls),ncol=n_simuls))
    # calculate proportion of Monte Carlo samples contained in the posterior support.
    return(mean(apply(simuls,1,bound)))
  }
}

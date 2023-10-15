library(pracma)
library(TreeTools)  # Only need this for "LnDoubleFactorial()"

# the recursive expression of f to compare with the explicit
get_rec <- function(d,x){
  if(d%%2 == 0){
    rec <- exp(0.5 * x^2) -1 ## on initialise  pour k = 2
    deb <- 4
  }else{
    rec <- x * sqrt(pi / 2) * erfi(x / sqrt(2))
    deb <- 3
  }
  if(d >= 3){
    for(k in seq(deb, d, 2)){
      rec <- exp(0.5 * x^2) - (k - 2) / x^2 * rec
    }
  }
  return(rec)  
}

#Explicit solution from lemma 1, used for 
get_rec_explicit = function(d,x){
  if(d%%2 == 0){
    rec <- exp(0.5 * x^2) - 1 ## on initialise  pour k = 2
    if(d == 2){
      return(rec)
    }
  }else{
    rec <- x * sqrt(pi / 2) * erfi(x / sqrt(2))
    if (d == 1){
      return(rec)
    }
  }
  coeff_1 = factorial2(d-2)
  coeff_2 = (-1/(x^2))^ceiling((d-2)/2)
  
  #using Adrian's trick for the sum
  log_series_1 = function(r){
    r*log(x^2) - LnDoubleFactorial(2*r+2*floor(d/2)-d)
  }# add signs later
  total = Vectorize(log_series_1)(1:ceiling((d-2)/2))
  max_of_total = max(total)
  total = total - max_of_total
  
  coeff_3 = sum((-1)^(1:ceiling((d-2)/2)) * exp(total)) * exp(max_of_total + x^2/2)
  # added the signs here
  # The sum can be problematic but it works fine if c_d=sqrt(d+L) for some L
  
  coeff_1 * coeff_2 * (rec + coeff_3)
  
}

# the log of the explicit solution: may avoid some numerical issues
log_get_rec_explicit = function(d,x){
  if(d%%2 == 0){
    rec <- exp(0.5 * x^2) - 1 ## on initialise  pour k = 2
    if(d == 2){
      return(log(abs(rec)))
    }
  }else{
    rec <- x * sqrt(pi / 2) * erfi(x / sqrt(2))
    if (d == 1){
      return(log(abs(rec)))
    }
  }
  log_coeff_1 = LnDoubleFactorial(d-2)
  log_coeff_2 = log(abs((-1/(x^2))))*ceiling((d-2)/2)
  series_1 = function(r){
    (-1/(x^2))^{-r}/factorial2(2*r+2*floor(d/2)-d)
  }
  
  # The sum can be problematic but it works fine if c_d=sqrt(d+L) for some L
  coeff_3 = exp(x^2/2) * sum(Vectorize(series_1)(1:ceiling((d-2)/2)))
  
  # We know for a fact that f is positive. 
  # That is why I can take the absolute value here.
  log_coeff_1 + log_coeff_2 + log(abs(rec + coeff_3))
  
}

#explicit expression of the scv for the plots
scv_rec_explicit = function(d,x){
  log_kappa_d_div_p_y = log(d)+(d/2)*log(2)+lgamma(d/2+1)
  exp(log_kappa_d_div_p_y - (d+2)*log(x) + log_get_rec_explicit(d,x))-1
}

#necessary condition for the 1st order condition to find c_d
get_deriv <- function(d, x, f){
  -2*d*f(d,x)+exp(0.5*x^2)*x^2
}

get_deriv_alt <- function(d, x, log_f){
  log(2*d)+log_f(d,x)-0.5*x^2-log(x^2)
}
# use this instead of the above function, much more stable
# also a necessary condition for the optimal value that this is 0

#not used so far: rough upper bound on the SCV in c=sqrt(d+1)
upper_bound = function(d){
  (2/(d+1))^{d/2}*gamma(d/2+1)*exp((d+1)/2)-1
}

#lower bound on the SCV for any choice of A; to be plotted in the figures
lower_bound = function(d){
  (2/d)^{d/2}*gamma(d/2+1)*exp(d/2)/2-1
}

#not used so far: the limiting function of the lower bound;
limit_lower_bound = function(d){
  sqrt((d+2)*pi/4)-1
}
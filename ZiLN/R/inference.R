
log_likelihood = function(par, a = a,  log_x = log_x) {
  mu = par[1]
  sigma = par[2]
  #log_x = log(x)
  n = length(log_x)
  nz_log_x = log_x[is.finite(log_x)]
  nz = sum(is.finite(log_x))
  ll = (n - nz) * log(pnorm(q = a, mean = mu, sd = sigma)) - 0.5 * sum((nz_log_x - mu)^2 / sigma^2) - nz * log(sigma)
  return(-ll)
}

z_d = function(z, mean = 0, sd = 1){
  return(z * dnorm(z,  mean, sd))
}

#' Hi, I'm not exported!
#'
#' @param txt       what to say
#'
#' @export
#' 
infer_Z = function(X, seq_depth = "TS") {
  
  d = dim(X)[2]
  n = dim(X)[1]

  X_clr = matrix(0, n, d)
  nz = X != 0
  
  if( seq_depth == "TS"){
    X_n = X / apply(X, 1, sum)
    for(i in 1:n) {
      X_clr[i,] = log(X_n[i,]) - mean(log(X_n[i,X_n[i,] != 0]))
    }
  }
  if(seq_depth == "unif"){
    for(k in 1:d){
      nz_X_k = X[nz[,k],k]
      emp = ecdf(c(nz_X_k, max(nz_X_k) + 1))
      X_clr[,k] = qnorm(emp(X[,k]))
    }
    X_clr[!nz] = -Inf
  }
  
  log_ratios = X_clr
  a = double(n)
  
  for(i in 1:d){
    a[i] = min(log_ratios[nz[,i],i])
  }
  
  Z = matrix(0, n, d)
  Z[] = log_ratios
  
  for(i in 1:d){
    if(sum(nz[,i]) > 1){
      mu = mean(log_ratios[is.finite(log_ratios)[,i],i])
      o = optim(c(0, 1), log_likelihood, log_x = X_clr[,i], a  =  a[i])
      par = o$par
      mu = par[1]
      sigma = par[2]
      exp_Z_0 = integrate(z_d, lower = -Inf, upper = a[i], mean = mu, sd =  sigma)
      Z_0 =  integrate(dnorm, lower = -Inf, upper = a[i], mean = mu, sd =  sigma)

      Z[X[,i] == 0, i] = exp_Z_0$value / Z_0$value
      # Handle small sigma
      if(sigma < 1e-5){
        Z[X[,i] == 0, i] = -10
      }
    }
    else if(sum(nz[,i]) == 1){
      Z[X[,i] == 0, i] = -10
    }
    else{
      Z[,i] = 0
    }
  }
  return(Z)
}


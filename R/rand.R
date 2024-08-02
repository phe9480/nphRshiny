#' Random Number Generator for Customized Cumulative Distribution Function
#' 
#' Suppose the CDF of a survival time T is provided, which has a domain of (min, max).
#' Then t = inverse.CDF (u) is a random variable following the CDF. Note: If a survival function S(t)
#' is provided, the CDF = 1 - S(t) is the CDF.
#' 
#' @param  n Sample size
#' @param  CDF Cumulative distribution function, e.g. function(t){1-exp(-log(2)/12*t)}
#' @param  min Lower bound of CDF's domain
#' @param  max Upper bound of CDF's domain
#' 
#' @return Generate a random sample of size n
#' 
#' @examples
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#' set.seed(2022)
#' 
#' #(1) Exponential distribution
#' t = rand(n=1000, CDF=function(t){1-exp(-log(2)/12*t)})
#' km <- survival::survfit(survival::Surv(t, rep(1, 1000)) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(2) Standard normal distribution
#' x = rand(n=1000, CDF=pnorm)
#' hist(x)
#' 
#' #(3) Uniform enrollment CDF = t/24 in 24 months.
#' x = rand(n=500, CDF=function(t){(t/24)*as.numeric(t <= 24)}, min=0, max=24)
#' hist(x)
#' 
#' @export
#' 
rand = function(n=1, CDF=pnorm, min=0, max=1e8) {

  #Inverse function of CDF
  inv.CDF = function(u){
      f = function(x){CDF(x)-u}
      t = uniroot(f=f, interval = c(min, max), tol=1e-3)$root
      return(t)
  }
  u = runif(n)
  ans = sapply(u, inv.CDF)
  
  return(ans)
}


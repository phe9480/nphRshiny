#' Generate Survival Data for Given Distributions or Customized Survival Function
#' 
#' Suppose the CDF of a survival time T is provided, which has a domain of (min, max).
#' Then t = inverse.CDF (u) is a random variable following the CDF. Note: If a survival function S(t)
#' is provided, the CDF = 1 - S(t) is the CDF.
#' 
#' @param  n Sample size
#' @param dist The following options are available. 
#' (1) "exponential": lam required;
#' (2) "weibull": shape and scale required. When shape = 1, it reduces to exponential distribution with hazard rate = 1/scale.
#' (3) "piecewise exponential": lamb and cuts required;
#' (4) "mixture cure rate of exponential": p1 and lam required; 
#' (5) "mixture cure rate of weibull": p1, shape, and scale required
#' (6) "customized": S as a survival function is required.
#'
#' @param lam hazard rate for exponential distribution
#' @param shape shape parameter for weibull distribution. Refer to rweibull() for details.
#' @param scale scale parameter for weibull distribution. Refer to rweibull() for details. 
#' @param p1 cure rate parameter for mixture cure rate distribution
#' @param S Survival function for customized distribution.
#' @param cuts Cut points for piecewise exponential distribution
#' 
#' @return Generate a random sample of size n
#' 
#' @examples
#' 
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#' set.seed(2022)
#' 
#' #(1) Exponential distribution
#' t = genSurv(dist = "exponential", n = 100, lam=log(2)/12) 
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(2) Weibull
#' t = genSurv(dist = "weibull", n = 100, shape=1.5, scale=12/log(2)) 
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(3) Piecewise exponential
#' t = genSurv(dist = "piecewise exponential", n = 100, lam=c(log(2)/10, log(2)/15), cuts = c(6))
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(4) Mixture cure rate of exponential
#' t = genSurv(dist = "mixture cure rate of exponential", n = 100, lam=log(2)/10, p1=0.15)
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(5) Mixture cure rate of weibull
#' t = genSurv(dist = "mixture cure rate of weibull", n = 100, shape=1.5, scale=12/log(2), p1=0.15)
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' #(6) Customized survival distribution
#' t = genSurv(dist = "customized", n = 100, S=function(t){exp(-log(2)/12*t)})
#' km <- survival::survfit(survival::Surv(t, rep(1, length(t))) ~ 1)
#' plot(km, xlab="Month Since Randomization",ylab="Survival")
#' 
#' 
#' @export
#' 
genSurv = function(dist = "exponential", n = 1, lam=log(2)/10, shape=NULL, scale=NULL,
                   p1=NULL, S=NULL, cuts=NULL) {

  #(1) Exponential
  if (dist == "exponential") {t = rexp(n, rate = lam)}
  
  #(2) Weibull
  if (dist == "weibull") {t = rweibull(n, shape = shape, scale = scale)}

  #(3) piecewise exponential
  if (dist == "piecewise exponential"){
    intervals = rep(NA, length(cuts)); intervals[1] = cuts[1]
    if (length(cuts) > 1){
      for (i in 2:length(cuts)){intervals[i] = cuts[i] - cuts[i-1]}
    }
    t = nphsim::rpwexp(n, rate = lam, intervals = intervals)
  }
  
  #(4) Mixture cure rate of exponential
  if (dist == "mixture cure rate of exponential") {
    t = rmcr(n=n, p=p1, alpha = lam, beta=1, gamma=1, lambda=0, tau=0, psi=1)
  }

  #(5) Mixture cure rate of weibull
  if (dist == "mixture cure rate of weibull") {
    t = rmcr(n=n, p=p1, alpha = scale^(-shape), beta=1, gamma=shape, lambda=0, tau=0, psi=1)
  }
  
  #(5) customized distribution
  if (dist == "customized") {
    CDF = function(t){1-S(t)}
    t = rand(n=n, CDF=CDF)
  }
    
  return(t)
}


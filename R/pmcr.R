#' Cumulative Distribution Function for Mixture Cure Rate Distribution With Delayed Effect
#' 
#' Consider the survival function of mixture cure rate (MCR) distribution: 
#' S(t) = p + (1-p)S0(t), where S0(t) is a survival distribution 
#' for susceptible subject, i.e., S0(0)=1 and S0(t) -> 0 as t -> Infinity.
#' S0(t) can be any proper survival function. For generality of S0(t), consider 
#' the generalized modified Weibull (GMW) distribution with parameters 
#' (alpha, beta, gamma, lambda) Martinez et al(2013). A delayed effect tau is also incorporated.
#' The delayed effect will modify the survival function after tau corresponding to the hazard ratio.
#' i.e. if psi = 0.6, and tau = 6 months, then after tau, the survival function is assumed having HR 0.6 compared to mcr(p, alpha, beta, gamma, lambda).
#' This feature is convenient when modeling delayed effect on top of MCR.
#' 
#' alpha: scale parameter
#' beta and gamma: shape parameters
#' lambda: acceleration parameter
#' 
#' S0(t) = 1-(1-exp(-alpha*t^gamma*exp(lambda*t)))^beta
#' 
#' Special cases:
#' (1) Weibull dist: lambda = 0, beta = 1. For traditional weibull parameterization (shape, scale), then alpha = scale^(-shape), gamma = shape. 
#' (2) Exponential dist: lambda = 0, beta = 1, gamma = 1. The hazard rate = alpha = 1/scale; and gamma=shape=1.
#' (3) Rayleigh dist: lambda = 0, beta = 1, gamma = 2.
#' (4) Exponentiated Weibull dist (EW): lambda = 0
#' (5) Exponentiated exponential dist (EE): lambda = 0 and gamma = 1
#' (6) Generalized Rayleigh dist (GR): lambda = 0, gamma = 2
#' (7) Modified Weibull dist (MW): beta = 1
#' 
#' Let T be the survival time according to survival function S1(t) below.
#' Denote T's distribution as MCR(p, alpha, beta, gamma, lambda, tau, psi), 
#' where tau is the delayed effect and psi is the hazard ratio after delayed effect,
#' ie. proportional hazards to S(t) after delay tau.
#' S(t) = p + (1-p)*S0(t)
#' S1(t) = S(t)I(t<tau) + S(tau)^(1-psi)*S(t)^psi
#' In brief, S0(t) is the proper GMW distribution (alpha, beta, gamma, lambda);
#' S(t) is MCR with additional cure rate parameter p;
#' In reference to S(t), S1(t) is a delayed effect distribution and proportional hazards after delay.
#' 
#' MCR(p, alpha, beta, gamma, lambda, tau, psi=1) reduces to S(t)
#' MCR(p, alpha, beta, gamma, lambda, tau=0, psi=0.6) means a distribution having proportional hazard HR=0.6 compared to MCR(p, alpha, beta, gamma, lambda).
#' MCR(p, alpha, beta, gamma, lambda, tau=6, psi=0.6) means a distribution having a delayed effect 6 months and afterwards proportional hazard HR=0.6 compared to MCR(p, alpha, beta, gamma, lambda).
#' MCR(p=0, alpha, beta, gamma, lambda, tau=0, psi=1) reduces to S0(t)
#' MCR(p=0, alpha, beta=1, gamma=1, lambda=0, tau=0, psi=1) reduces to exponential dist.
#' 
#' @param  p Cure rate parameter. When p = 0, it reduces to GMW(alpha, beta, gamma, lambda) distribution.
#' @param  alpha Generalized modified Weibull (GMW) distribution parameters. alpha > 0
#' @param  beta Generalized modified Weibull (GMW) distribution parameters. beta > 0
#' @param  gamma Generalized modified Weibull (GMW) distribution parameters. gamma >= 0 but gamma and lambda cannot be both 0.
#' @param  lambda Generalized modified Weibull (GMW) distribution parameters. lambda >= 0 but gamma and lambda cannot be both 0.
#' @param  tau  Threshold for delayed effect period. tau = 0 reduces no delayed effect.
#' @param  psi  Hazard ratio after delayed effect. psi = 1 reduces to the survival function without proportional hazards after delayed period (tau).
#' 
#' @return Cumulative distribution function
#' 
#' @examples
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#' #(1) Exponential distribution mixture with cure rate 0.15
#' t=seq(1,100)
#' c1 = sapply(t, pmcr, p=0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0)
#' plot(t, 1-c1, type="l", ylab="", main="Survival Function", ylim=c(0,1))
#' 
#' #(2) Weibull distribution mixture with cure rate 0.15
#' c2 = sapply(t, pmcr, p=0.15, alpha = log(2)/12, beta=1, gamma=0.8, lambda=0)
#' plot(t, 1-c2, type="l", ylab="", main="Survival Function", ylim=c(0,1))
#' 
#' #(3) Rayleigh distribution mixture with cure rate 0.15
#' c3 = sapply(t, pmcr, p=0.15, alpha = log(2)/12, beta=1, gamma=2, lambda=0)
#' plot(t, 1-c3, type="l", ylab="", main="Survival Function", ylim=c(0,1))
#' 
#' #(4) GMW distribution mixture with cure rate 0.15
#' c4 = sapply(t, pmcr, p = 0.15, alpha = log(2)/12, beta=1, gamma=1, lambda=0.2)
#' plot(t, 1-c4, type="l", ylab="", main="Survival Function", ylim=c(0,1))
#' 
#' #(5) Weibull distribution mixture with cure rate 0.15 with delayed effect 6 months
#' c5 = sapply(t, pmcr, p=0.15, alpha = log(2)/12, beta=1, gamma=0.8, lambda=0, tau=6, psi=0.6)
#' plot(t, 1-c5, type="l", ylab="", main="Survival Function", ylim=c(0,1))
#' lines(t, 1-c2, col=2)
#' 
#' @export
#' 
pmcr = function(t, p=0.3, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1) {
  
  S0 = function(t){1-(1-exp(-alpha*t^gamma*exp(lambda*t)))^beta}
  S  = function(t){p + (1-p)*S0(t)}
  
  F0 = function(t){1 - S0(t)}
  Ft = (1-p)*F0(t)
  
  f = 1-S(t); f2 = 1 - S(tau)^(1-psi)*S(t)^psi
  ans = f*as.numeric(t <= tau) + f2*as.numeric(t > tau)
  
  return(ans)
}


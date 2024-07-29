#' Power Calculation Using Haybittle-Peto Boundary and Logrank Test in Group Sequential Design
#' 
#' This function calculates the power for the given number of events and rejection boundary
#' at each analysis for group sequential design
#' 
#' @param p1      A fixed p value boundary for IAs
#' @param overall.alpha Overall alpha for the test in the group sequential design. 
#' Default, overall.alpha = 0.025. overall.alpha must be one-sided. 
#' @param events  A vector of target events for IAs and FA
#' @param events0 A vector of target events for IAs and FA for control arm. 
#' Required when variance option is H1.
#' @param events1 A vector of target events for IAs and FA for experimental arm
#' Required when variance option is H1.
#' @param hr      Hazard ratio (smaller better)
#' @param r       Randomization ratio. r=1 means 1:1; r=2 means 2:1.
#' @param variance Option for variance estimate. "H1" or "H0". Default H1, which is usually more conservative than H0.
#'  
#' @return Objects of two dataframes alpha and design
#'  \itemize{
#'  \item alpha      Incremental alpha
#'  \item cum.alpha  Cumulative alpha
#'  \item overall.alpha    Overall alpha
#'  \item side       Always one-sided
#'  \item events     Number of events by each analysis
#'  \item information Information fraction by each analysis
#'  \item p          Rejection boundary in p value by each analysis
#'  \item z          Rejection boundary in z value by each analysis
#'  \item cv         Critical value (minimum detectable difference) in hazard ratio
#'  \item marg.power Marginal power, defined as the power at each analysis assuming as a separate study. Note: This is for reference only and it is not consistent with group sequential design.
#'  \item incr.power Incremental power by each analysis
#'  \item cum.power Cumulative power by each analysis
#'  \item overall.power Overall Power
#'  }
#' @examples
#' 
#' #(1) Haybittle-Peto method for analyses performed at 96, 153, and 250 events with 1-sided p-value rejection boundary 0.0003 at IAs
#' lr.power.HaybittlePeto(p1=0.0003, overall.alpha=0.025, hr = 0.7, r = 1, events=c(96, 153, 250))
#' 
#' #(2) Haybittle-Peto method for analyses performed at 229 and 340 events with 1-sided p-value rejection boundary 0.0003 at IA
#' lr.power.HaybittlePeto(p1=0.0003, overall.alpha=0.025, hr = 0.7, r = 1, events=c(229, 340))
#' 
#' #(3) Haybittle-Peto method for analyses performed at 210 and 280 events with 1-sided p-value rejection boundary 0.0003 at IA, variance under H1
#' lr.power.HaybittlePeto(p1=0.0003, overall.alpha=0.025, hr = 14/22, r = 1, events=c(210, 280))
#'                      
#' @export
#'
lr.power.HaybittlePeto <- function(p1=0.0003, overall.alpha=0.025, 
                          events=c(96, 153, 250), 
                          events0=NULL, events1=NULL, 
                          hr = 0.6, r = 1, variance="H0"){
  #p1: 1-sided boundary for all IAs. 
  #overall.alpha: 1-sided overall type I error, default 0.025
  #timing: information fractions for performing the analyses
  
  r1 = r/(r+1) #Proportion of experimental arm subjects among all subjects
  
  M = length(events) #number of analyses
  frac = events/events[M]
  
  #Correlation structure
  corr = matrix(1, nrow = M, ncol = M)
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      corr[i, j] = corr[j, i] = sqrt(frac[i]/frac[j])
    }
  }
  
  if (M == 1) {
    pf = p1 = overall.alpha
    zf = z1 = qnorm(1 - overall.alpha)
  }
  else {
    z1 = qnorm(1-p1)
    a = rep(NA, M) #incremental alpha
    if (M == 2) {
      a[1] = p1; a[2] = overall.alpha-a[1]
    }
    if (M > 2){
      a[1] = p1
      for (i in 2:(M-1)){
        a[i] = mvtnorm::pmvnorm(lower = c(rep(-Inf, i - 1), z1), 
                                upper = c(rep(z1, i - 1), Inf), 
                                corr = corr[1:i, 1:i], 
                                abseps = 1e-08, 
                                maxpts = 1e+05)[1]
      }
      a[M] = overall.alpha - sum(a[1:(M-1)])
    }
    f.x = function(x) {
      I = mvtnorm::pmvnorm(lower = c(rep(-Inf, M - 1), x), 
                           upper = c(rep(z1, M - 1), Inf), 
                           corr = corr, 
                           abseps = 1e-08, 
                           maxpts = 1e+05)[1]
      return(I - a[M])
    }
    zf = uniroot(f = f.x, interval = c(1, 10), tol = 1e-08)$root
    pf = 1 - pnorm(zf)
  }
  p = c(rep(p1, M-1), pf)
  
  # CV calculation
  if (variance=="H1"){ 
    events = events1 + events0
    re = events1/events
  }
  #(5) Critical Values in HR
  #######By default, use variance under H1 for more conservative estimate of power
  z = c(rep(z1, M-1), zf)
  if (variance == "H0"){
    cv = exp(-z/sqrt(r1*(1-r1)*events))
  } else {   
    cv = exp(-z/sqrt(re*(1-re)*events))
  }
  
  #(1) Mean of z statistics
  #######By default, use variance under H1 for more conservative estimate of power
  if (variance == "H0"){
    mu = -log(hr) * sqrt(r1*(1-r1)*events)
  } else {
    #Proportion of events
    mu = -log(hr) * sqrt(re*(1-re)*events)
  }
  
  #(4)Marginal Power
  marg.power = rep(NA, M)
  
  for (k in 1:M){
    marg.power[k] = 1-pnorm(qnorm(1-p[k]), mean=mu[k])
  }
  
  #(4c)Incremental Power
  incr.power = rep(NA, M)
  cum.power =  rep(NA, M)
  cum.power[1] = incr.power[1] = marg.power[1]
  
  
  if(M > 1){
    for (k in 2:M){
      incr.power[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), z[k]), 
                                       upper = c(z[1:(k-1)], Inf), mean=mu[1:k], 
                                       corr = corr[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      cum.power[k] = cum.power[k-1] + incr.power[k]
    }
    overall.power = cum.power[M]
  } else{
    cum.power = incr.power = marg.power = overall.power = 1-pnorm(z, mean=mu)
  }
  side = 1; alpha = a; information = frac
  cum.alpha = rep(NA, M)
  for (i in 1:M){cum.alpha[i] = sum(alpha[1:i])}
  
  o <- list()
  o$alpha = data.frame(cbind(alpha, cum.alpha, overall.alpha, side))
  o$design = data.frame(cbind(events, information, p, z, cv, marg.power,incr.power,cum.power,overall.power))
  return(o)
}


#' Determine the Expected Date of Cutoff (DCO) According to the Expected Numbers of Events
#' 
#' This function calculates the expected DCO according to the number of events, where the DCO is calculated from first subject randomized. 
#' The function returns the expected number of events for each arm at time t, based on the provided
#' enrollment distribution function and random lost-to-followup distribution if applicable.
#' 
#' @param n Total sample size for two arms.
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param Lambda Distribution function of enrollment. For uniform enrollment, 
#' Lambda(t) = (t/A) where A is the enrollment period, i.e., Lambda(t) = t/A for 0<=t<=A, and 
#' Lambda(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' Lambda(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default Lambda is uniform distribution function.
#' @param G0 Cumulative distribution function of drop-off for control arm. 
#' For example, 3 percent drop-off in 12 months of followup means then the hazard rate for unit time 
#' is eta0 = -log(1-0.03/12), so G0=function(t){1-exp(-eta0*t)}.
#' @param G1 Cumulative distribution function of drop-off for experimental arm. 
#' Similarly, G1=function(t){1-exp(-eta1*t)}, where eta1=-log(1-0.03/12) is the 
#' hazard rate for 3 percent drop-off in 12 months of followup.  
#' 
#' @return Date of Cutoff (DCO)
#'  
#' @examples 
#' #Example (1) Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' #6 months after last patient randomized.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' #Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' max.DCO = 50
#' nE = matrix(NA, nrow = max.DCO, ncol=3)
#' for (DCO in 1:max.DCO) {
#'   o = fe(DCO = DCO, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#'   nE[DCO, 1] = o$e0; nE[DCO, 2] = o$e1; nE[DCO, 3] = o$e
#' }
#' plot(1:max.DCO, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:max.DCO, nE[, 1], lty = 1, col=1)
#' lines(1:max.DCO, nE[, 2], lty = 2, col=2)
#' lines(1:max.DCO, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' e=fe(DCO = 30, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#' e  
#' fDCO(events = e$e, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#' 
#' fDCO(events = c(250, 300), r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450)
#' 
#' #Example (2) Same trial set up as example (1) but assuming delayed effect for 
#' #experimental arm. The delayed period is assumed 6 months, and after delay the
#' #hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' #Hazard function and survival function for experimental arm
#' h1 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' #dropout probability: 0.05 in 12 months for control arm. 0.03 for experimental arm
#' eta0 = -log(1-0.05/12)
#' eta1 = -log(1-0.03/12)
#' G0=function(t){1-exp(-eta0*t)}
#' G1=function(t){1-exp(-eta1*t)}
#' 
#' max.DCO = 30
#' nE = matrix(NA, nrow = max.DCO, ncol=3)
#' for (DCO in 1:max.DCO) {
#'   o = fe(DCO = DCO, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450, G0=G0, G1=G1)
#'   nE[DCO, 1] = o$e0; nE[DCO, 2] = o$e1; nE[DCO, 3] = o$e
#' }
#' plot(1:max.DCO, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:max.DCO, nE[, 1], lty = 1, col=1)
#' lines(1:max.DCO, nE[, 2], lty = 2, col=2)
#' lines(1:max.DCO, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' e=fe(DCO = 30, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450, G0=G0, G1=G1)
#' e
#' fDCO(events = e$e, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, n = 450, G0=G0, G1=G1)
#' 
#' @export
fDCO = function(n = 450, events = c(250, 300), r = 1, 
              h0 = function(t){log(2)/12}, 
              S0=function(t){exp(-log(2)/12 * t)}, 
              h1 = function(t){log(2)/12*0.70}, 
              S1= function(t){exp(-log(2)/12 * 0.7 * t)},
              Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
              G0 = function(t){0}, G1 = function(t){0}){

  K = length(events) #number of data cutoffs
  dco = rep(NA, K)
  
  for (i in 1:K){
    f.e = function(x){
      e = fe(DCO = x, r=r, h0=h0, S0=S0, h1=h1, S1=S1, Lambda=Lambda, n=n, G0=G0, G1=G1)$e
      return(e - events[i])
    }
    dco[i] = uniroot(f=f.e, lower=0, upper=1000)$root  
  }
  return(dco)
}

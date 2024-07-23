#'  Expected Average Weighted Hazard Ratio from Weighted Cox-Regression at a Calendar Time
#' 
#'  This function calculates the expected average weighted HR under H1
#'  evaluated at a calendar time using (He et al 2021) method and Kalbfleisch and Prentice (1981) method.
#'  Both methods are very close and the difference is usually negligible.
#'  It's shown that the expected average weighted HR is approximately the weighted
#'  average of piecewise HRs if the true HR(t) is piecewise under H1, where the weight
#'  at each subinterval is the weighted probability of event at the subinterval 
#'  divided by the cumulative probability of event at a calendar time. For general form of H1 in terms
#'  of HR(t), the expected average weighted HR is 
#'  the expectation of weighted hazard ratio function divided by the weighted cumulative
#'  probability of event at a calendar time. For the special case of Cox
#'  regression, the weight is 1, and the expected average hazard ratio reduces 
#'  to the expected hazard ratio function evaluated at the calendar time. For proportional hazards, it further reduces to
#'  the constant hazard ratio under H1 regardless of weight functions.
#'  
#'  This function allows flexible weight function as used in the weighted log-rank test. This 
#'  function also allows flexible alternative hypothesis in terms of HR(t), the 
#'  instaneous hazard ratio function over time. For delayed effect scenario under H1,
#'  one can define HR(t) as a piecewise constant function of survival time t.
#'  In addition, the function can handle user-defined flexible non-uniform enrollment 
#'  distribution function and independent time to lost-to-followup process which 
#'  is user-defined function of any lost-to-followup pattern such as constant
#'  lost-to-followup rate or Weibull distribution. For most common setting
#'  in practice, assuming the same lost-to-followup pattern in both arms.
#'  If the total number of subjects n is not provided, the function returns
#'  the non-centrality parameter of n^(-1/2)*Z where Z is the normalized 
#'  weighted log-rank test statistic Z = U/sqrt(var(U)) with U as the weighted
#'  log-rank score statistic. 
#'  
#' @param DCO Calendar times, calculated from first subject randomization date, 
#'           when the weighted log-rank test will be evaluated at the same calendar time. 
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param n Total sample size for two arms. Default is NULL. 
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param  rho Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#' @param  gamma Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  s.tau  Survival rate cut S(tau) at t = tau1; default 0.5, ie. cut at median.
#'       s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  f.ws  Self-defined weight function of survival rate. 
#'         For example, f.ws = function(s){1/max(s, 0.25)}
#'         When f.ws is specified, the weight function takes them as priority.
#' @param Lambda Distribution function of enrollment. For uniform enrollment, 
#' Lambda(t) = (t/A) where A is the enrollment period, i.e., Lambda(t) = t/A for 0<=t<=A, and 
#' Lambda(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' Lambda(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default Lambda is uniform distribution function.
#' @param G Distribution function of lost-to-follow-up censoring process. The observed
#' survival time is min(survival time, lost-to-follow-up time). Default G = 0 (no lost-to-followup)
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  AHR: Expected average hazard ratio
#'  \item  AHR.KP Kalbfleisch and Prentice (1981) method
#'  \item  var.logHR: Asymptotic variance of the log expected average hazard ratio
#'  \item  var.AHR: Asymptotic variance of the expected average hazard ratio using delta method
#'  \item  var.AHR.KP: Asymptotic variance of the expected average hazard ratio of 
#'  Kalbfleisch and Prentice (1981) method using delta method
#'  \item  wt: Weight function used
#'  }
#'  
#' @examples 
#' #############
#' #Example (1) 
#' #############
#' #Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' #6 months after last patient randomized.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' lambda1 = lambda0 * HR
#' f.logHR.PH = function(t){log(HR)}
#' h1.PH = function(t){lambda1}; S1.PH= function(t){exp(-lambda1 * t)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' rho = 0; gamma = 0; tau = NULL; s.tau = 0; f.ws = NULL; G = function(t){0}
#' DCO = 24; r = 1; n = 450
#' 
#' wlr.AHR(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G)
#'      
#' #If proportional hazards, then the average HR is independent of any weight functions       
#' wlr.AHR(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'      rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G)
#'      
#' #############      
#' #Example (2) 
#' #############
#' #Same trial set up as example (1) but assuming delayed effect for 
#' experimental arm. The delayed period is assumed 6 months, and after delay the
#' hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' h1.D6 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1.D6 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HR)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G = function(t){0}
#' 
#' wlr.AHR(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G)
#'      
#' wlr.AHR(DCO = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G)
#'      
#' #############
#' #Example (3). 
#' #############
#' #Draw a plot of average weighted HR over time    
#' plot.AHR = function(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1, S1=S1, s.tau=0.5){
#'   maxT = 48; AHR.lr = AHR.fh = AHR.fh11 = AHR.sfh = rep(NA, maxT)
#'   for(i in 1:maxT){
#'     AHR.lr[i] = wlr.AHR(DCO = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)$AHR
#'     AHR.fh[i] = wlr.AHR(DCO = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)$AHR      
#'     AHR.sfh[i] = wlr.AHR(DCO = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = s.tau, f.ws = NULL,
#'        Lambda = Lambda, G = G)$AHR    
#'     AHR.fh11[i] = wlr.AHR(DCO = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 1, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)$AHR          
#'   }
#'
#'   maxY = max(c(AHR.lr, AHR.fh, AHR.sfh))
#'   plot(1:maxT, AHR.lr, type="n", ylim=c(0.4, 1), 
#'   xlab="Months", ylab = "Average Weighted HR")   
#'   lines(1:maxT, AHR.lr, lty = 1, col=1)
#'   lines(1:maxT, AHR.fh, lty = 2, col=2)
#'   lines(1:maxT, AHR.sfh, lty = 3, col=3)
#'   lines(1:maxT, AHR.fh11, lty = 4, col=4)
#'   abline(h=0.8, col="gray", lty=2)
#'   
#'   legend(0, 0.65, c("Log-rank", "FH(0,1)", paste("sFH(0,1,",s.tau,")", sep=""),"FH(1,1)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.8)
#' }
#' 
#' #(a)Proportional Hazards
#' #########################
#' plot.AHR(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1.PH, S1=S1.PH)
#' 
#' #(b)Delayed effect of 6 months
#' #########################
#' plot.AHR(f.logHR = f.logHR.D6, h0=h0, S0=S0, h1=h1.D6, S1=S1.D6, s.tau=0.5)
#' 
#' #(c)Crossover to effective subsequent IO therapy and delay 6 months
#' #########################
#' crossDCO = 24; HRx = 0.9; #after crossover, assuming the tail piecewise HR = 0.9.
#' h0x = function(t){lambda0*as.numeric(t < crossT) + HR/HRx*lambda0*as.numeric(t >= crossT)}; 
#' c0 = exp(-crossT*lambda0*(1-HR/HRx)); 
#' S0x = function(t){exp(-lambda0*t)*as.numeric(t<crossT) + c0*exp(-HR/HRx*lambda0*t)*as.numeric(t>=crossT)}
#' h1.D6x = function(t){lambda0*as.numeric(t < delay) + HR*lambda0*as.numeric(t >= delay)}
#' c1 = exp(-delay*lambda0*(1-HR)); 
#' S1.D6x = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c1*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' f.logHR.D6x = function(t){log(as.numeric(t < delay) 
#' + as.numeric(t >= delay & t < crossT)*HR + as.numeric(t >= crossT)*HRx)}
#' 
#' plot.AHR(f.logHR = f.logHR.D6x, h0=h0x, S0=S0x, h1=h1.D6x, S1=S1.D6x, s.tau=0.5)
#' 
#' #############
#' #Example (4)
#' #Comparison of KP method and (He et al) method 
#' #Delayed effect 6 mo
#' #############
#' HRseq = seq(0.1, 0.9, by=0.01); L = length(HRseq)
#' AHR.lr = AHR.KP.lr = AHR.fh = AHR.KP.fh = rep(NA, L)
#' for (i in 1:L){
#'   HR = HRseq[i]
#'   h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#'   h1.D6 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#'   c = exp(-delay*lambda0*(1-HR)); 
#'   S1.D6 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#'   f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HR)}
#'   
#'   lr = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   fh = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   AHR.lr[i] = lr$AHR; AHR.KP.lr[i] = lr$AHR.KP
#'   AHR.fh[i] = fh$AHR; AHR.KP.fh[i] = fh$AHR.KP
#' }
#'   plot(HRseq, AHR.lr, type="n", ylim=c(0, 1), 
#'   xlab="Hazard Ratio", ylab = "Average Weighted HR")   
#'   lines(HRseq, AHR.lr, lty = 1, col=1)
#'   lines(HRseq, AHR.KP.lr, lty = 2, col=2)
#'   
#'   lines(HRseq, AHR.fh, lty = 3, col=3)
#'   lines(HRseq, AHR.KP.fh, lty = 4, col=4)
#'   
#'   legend(0.6, 0.5, c("Log-rank: AHR", "Log-rank: AHR(KP)", "FH01: AHR", "FH01: AHR(KP)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.7)
#' ####################
#' #Example (5)
#' #Comparison of KP method and (He et al) method 
#' #Proportional Hazards. Both methods equivalent
#' ####################
#' HRseq = seq(0.1, 0.9, by=0.01); L = length(HRseq)
#' AHR.lr = AHR.KP.lr = AHR.fh = AHR.KP.fh = rep(NA, L)
#' for (i in 1:L){
#'   HR = HRseq[i]
#'   lambda1 = lambda0 * HR
#'   f.logHR.PH = function(t){log(HR)}
#'   h1.PH = function(t){lambda1}; 
#'   S1.PH= function(t){exp(-lambda1 * t)}
#'   
#'   lr = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   fh = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   AHR.lr[i] = lr$AHR; AHR.KP.lr[i] = lr$AHR.KP
#'   AHR.fh[i] = fh$AHR; AHR.KP.fh[i] = fh$AHR.KP
#' }
#'   plot(HRseq, AHR.lr, type="n", ylim=c(0, 1), 
#'   xlab="Hazard Ratio", ylab = "Average Weighted HR")   
#'   lines(HRseq, AHR.lr, lty = 1, col=1)
#'   lines(HRseq, AHR.KP.lr, lty = 2, col=2)
#'   
#'   lines(HRseq, AHR.fh, lty = 3, col=3)
#'   lines(HRseq, AHR.KP.fh, lty = 4, col=4)
#'   
#'   legend(0.6, 0.5, c("Log-rank: AHR", "Log-rank: AHR(KP)", "FH01: AHR", "FH01: AHR(KP)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.7)
#' 
#' ####################
#' #Example (6)
#' #Comparison of KP method and (He et al) method 
#' #Crossover to subsequent therapy
#' ####################
#' HRseq = seq(0.1, 0.9, by=0.01); L = length(HRseq)
#' AHR.lr = AHR.KP.lr = AHR.fh = AHR.KP.fh = rep(NA, L)
#' crossT = 30; HRx = 0.9; #after crossover, assuming the tail piecewise HR = 0.9.
#' for (i in 1:L){
#'   HR = HRseq[i]
#'   lambda1 = lambda0 * HR
#'   h0x = function(t){lambda0*as.numeric(t < crossT) + HR/HRx*lambda0*as.numeric(t >= crossT)}; 
#'   c0 = exp(-crossT*lambda0*(1-HR/HRx)); 
#'   S0x = function(t){exp(-lambda0*t)*as.numeric(t<crossT) + c0*exp(-HR/HRx*lambda0*t)*as.numeric(t>=crossT)}
#'   h1.D6x = function(t){lambda0*as.numeric(t < delay) + HR*lambda0*as.numeric(t >= delay)}
#'   c1 = exp(-delay*lambda0*(1-HR)); 
#'   S1.D6x = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c1*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#'   f.logHR.D6x = function(t){log(as.numeric(t < delay) 
#'     + as.numeric(t >= delay & t < crossT)*HR + as.numeric(t >= crossT)*HRx)}
#'   lr = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   fh = wlr.AHR(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        Lambda = Lambda, G = G)
#'   AHR.lr[i] = lr$AHR; AHR.KP.lr[i] = lr$AHR.KP
#'   AHR.fh[i] = fh$AHR; AHR.KP.fh[i] = fh$AHR.KP
#' }
#'   plot(HRseq, AHR.lr, type="n", ylim=c(0, 1), 
#'   xlab="Hazard Ratio", ylab = "Average Weighted HR")   
#'   lines(HRseq, AHR.lr, lty = 1, col=1)
#'   lines(HRseq, AHR.KP.lr, lty = 2, col=2)
#'   
#'   lines(HRseq, AHR.fh, lty = 3, col=3)
#'   lines(HRseq, AHR.KP.fh, lty = 4, col=4)
#'   
#'   legend(0.5, 0.5, c("Log-rank: AHR", "Log-rank: AHR(KP)", "FH01: AHR", "FH01: AHR(KP)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.7)
#' 
#' @export
#' 
wlr.AHR = function(DCO = 24, r = 1, n = 450, 
       h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
       h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
       f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
       rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
       Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
       G = function(t){0}){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 

  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}  

  #UNDER H1
  f.bar = function(t){r0 * f0(t) + r1 * f1(t)}
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)}
  
  #Weight function
  f.w = function(t, f.S = S0, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma){
    s = f.S(t)
    #First priority: f.ws
    if(!is.null(f.ws)){
      w = f.ws(s)
    }else {
      #Second priority: s.tau
      if (!is.null(s.tau)){
        s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max);
      } else {
        s.til = apply(cbind(s, f.S(tau)), MARGIN=1,FUN=max);        
      }
      w = s.til^rho*(1-s.til)^gamma
    }
    return(w)
  }
  
  #eta: weighted average of HR(t) by prob. of event
  I.eta = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * f.logHR(t)*Lambda(DCO-t) * (1 - G(t)) * f.bar(t))
  }
  eta = integrate(I.eta, lower=0, upper=T, abs.tol=1e-8)$value
  
  #d: prob. of event
  I.d = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * Lambda(DCO-t) * (1 - G(t)) * f.bar(t))
  }
  d = integrate(I.d, lower=0, upper=T, abs.tol=1e-8)$value
  AHR = exp(eta / d)  

  #Kalbfleisch and Prentice (1981) Approach
  I.KP1 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    HRt = exp(f.logHR(t))
    return(w *HRt/(1 + HRt)*Lambda(DCO-t) * (1 - G(t)) * f.bar(t))
  }
  I.KP2 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    HRt = exp(f.logHR(t))
    return(w *1/(1 + HRt)*Lambda(DCO-t) * (1 - G(t)) * f.bar(t))
  }
  AHR.KP = integrate(I.KP1, lower=0, upper=T, abs.tol=1e-8)$value/integrate(I.KP2, lower=0, upper=T, abs.tol=1e-8)$value
  
  #Asymptotic variance of AHR
  I.I0 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w^2 * Lambda(DCO-t) * (1 - G(t)) * f.bar(t))
  }
  I0 = integrate(I.I0, lower=0, upper=T, abs.tol=1e-8)$value

  I.d0 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (exp(f.logHR(t))-1)*Lambda(DCO-t) * (1 - G(t)) * f0(t))
  }
  d0 = integrate(I.d0, lower=0, upper=T, abs.tol=1e-8)$value
  I.d1 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (1-exp(-f.logHR(t)))*Lambda(DCO-t) * (1 - G(t)) * f1(t))
  }
  d1 = integrate(I.d1, lower=0, upper=T, abs.tol=1e-8)$value
  AHR2 = exp((r1*d1 + r0*d0)/d)

  V = I0 / (n*Lambda(T)*r0*r1*d^2)
  
  o = list()
  o$AHR = AHR
  o$AHR.KP= AHR.KP
  o$AHR2 = AHR2
  o$var.logHR = V
  o$var.AHR = exp(2*AHR)*V
  o$var.AHR.KP = exp(2*AHR.KP)*V
  if(!is.null(f.ws)){wt = f.ws} else{
    wt = data.frame(cbind(rho, gamma, tau, s.tau))
  }
  o$wt = wt
  return(o)
}

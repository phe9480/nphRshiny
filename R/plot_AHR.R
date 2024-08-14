#'  Display Expected Average Weighted Hazard Ratio Over Time
#' 
#'  This function calculates the expected average weighted HR under H1
#'  evaluated at a time since 1st subject randomized. Three methods are implemented: 
#'  (1) Geometric Schoenfeld; (2) Kalbfleisch and Prentice (1981) method; (3) Geometric.
#'  These methods are very close and the difference is usually negligible.
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
#' @param Tmax  Maximum range of calendar times. Default Tmax = 50.
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
#' @param G Cumulative distribution function of drop-off. 
#' For example, 3 percent drop-off in 12 months of followup means then the hazard rate for unit time 
#' is eta = -log(1-0.03/12), so G=function(t){1-exp(-eta*t)}. If control and experimental arms have separate dropouts, then create G(t)=(1-r)*G0(t)+r*G1(t)
#' @param  method Methods to calculate the AHR: Options include "Geometric Schoenfeld", "Geometric",
#'                and "Kalbfleisch and Prentice". Default is "Geometric Schoenfeld".
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
#' #Assuming delayed effect 6 months, and after delay the hazard ratio is assumed 0.65.
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
#' ahr1=plot_AHR(n = 450, Tmax = 50, r = 1, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G, method="Kalbfleisch and Prentice")
#' 
#' ahr2=plot_AHR(n = 450, Tmax = 50, r = 1, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G, method="Geometric")
#'      
#' ahr3=plot_AHR(n = 450, Tmax = 50, r = 1, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G, method="Geometric Schoenfeld")
#' 
#' ahr4=plot_AHR(n = 450, Tmax = 50, r = 1, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = function(s){1/s},
#'      Lambda = Lambda, G = G, method="Kalbfleisch and Prentice")
#'      
#' #control: Weibull (shape = 2, scale = 20); experimental arm: exp(rate = 0.03)
#' ahr5=plot_AHR(n = 450, Tmax = 30, r = 1, h0 = function(t){0.005*t}, S0=function(t){exp(-(t/20)^2)},
#'      h1 = function(t){0.03}, S1=function(t){exp(-0.03*t)}, 
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      Lambda = Lambda, G = G, method="Kalbfleisch and Prentice")
#'      
#' @export
#' 
plot_AHR = function(n = 450, Tmax = 50, r = 1,  
                   h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
                   h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                   rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
                   Lambda = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                   G = function(t){0}, method="Geometric Schoenfeld", ...){
  
  t = seq(1, Tmax, by = 1); ahr = rep(NA, length(t))
  col.seq = c("seagreen3","blue3","turquoise4","deeppink3","orange")
  
  f.logHR = function(u){log(h1(u)/h0(u))}
  
  for (i in 1:length(t)){
    #print(i)
    ahr.i=wlr.AHR(DCO=t[i], r=r, n = n, h0=h0, S0=S0, h1=h1, S1 = S1, f.logHR = f.logHR,
                  rho=rho, gamma=gamma, tau=tau, s.tau=s.tau, f.ws=f.ws,
                  Lambda = Lambda, G = G)
    
    if (method == "Geometric Schoenfeld") {ahr[i] = ahr.i$AHR}
    if (method == "Geometric") {ahr[i] = ahr.i$AHR2}
    if (method == "Kalbfleisch and Prentice") {ahr[i] = ahr.i$AHR.KP}
  }
  ahr = round(ahr, 3)
  plot(t, ahr, type="n", bty = "l", xlab="Time Since 1st Subject Randomized (mo)", 
       ylab="Expected Average Hazard Ratio", ...)
  lines(t, ahr, lwd=3, col=col.seq[1], lty=1)
  abline(h = seq(0, 1, 0.05), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  return(ahr)
}


#' Display of Cumulative Number of Events Over Time per Study Design
#' 
#' This function plots the cumulative number of events for each arm and total per study design
#' 
#' @param n Total sample size from 2 arms
#' @param Tmax Maximum range of time since 1st subject randomized for plot
#' @param DCO  A sequence of calendar time for the plot (x-axis)
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
#' @param G0 Distribution function of lost-to-follow-up censoring process for control arm.
#' @param G1 Distribution function of lost-to-follow-up censoring process for experimental arm.
#' @param ... Other graphic parameters passed to the plot
#'
#' @return Display of the graph
#'  
#' @examples 
#' 
#' HR = 0.65; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' 
#' #Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' G0 = function(t){1-exp(-0.03/12*t)} #exponential distribution CDF for 3% dropoff every 12 months
#' G1 = function(t){1-exp(-0.10/12*t)} #exponential distribution CDF for 10% dropoff every 12 months
#' 
#' leg = list(x=0, y=80, txt=c("Control", "Exp. Arm", "Total"))
#' par = list(main="OS: Number of Events Over time",
#'            xlab="Calendar Time (mo)", 
#'            ylab="Events")
#'            
#' plot_events(n = 450, Tmax = 50, r=1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, G0 = G0, G1 = G1)
#' 
#' #Specify graphic parameters 
#' plot_events(n = 450, Tmax = 50, r=1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, Lambda = Lambda, G0 = G0, G1 = G1, par=par)
#' 
#' #Understanding the drop off effect from two arms
#' plot_events(n = 450, Tmax = 50, r=1, h0 = h0, S0 = S0, h1 = h0, S1 = S0, Lambda = Lambda, G0 = G0, G1 = G1)
#' 
#' @export
#' 
plot_events = function(n = 450, Tmax = 50, r=1, 
                       h0=function(t){log(2)/12}, S0=function(t){exp(-log(2)/12 * t)}, 
                       h1=function(t){log(2)/12*0.7}, S1=function(t){exp(-log(2)/12 *0.7* t)},
                Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                G0 = function(t){1-exp(-0.03/12*t)}, 
                G1 = function(t){1-exp(-0.10/12*t)}, 
                leg=list(x=0, y=n/2,txt=c("Control", "Exp. Arm", "Total")), 
                param=list(xlab = "Time Since 1st Subject Randomized (mo)", 
                         ylab = "Cumulative Events",
                         main = "")){
  t = seq(1, Tmax, by = 1)
  nE = matrix(NA, nrow = length(t), ncol=3)
  for (i in 1:length(t)) {
    o = fe(n = n, DCO = t[i], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                  Lambda = Lambda, G0 = G0, G1 = G1)
    nE[i, 1] = o$e0
    nE[i, 2] = o$e1
    nE[i, 3] = o$e
  }

  plot(t, nE[, 3], type = "n", xlab=param$xlab, ylab=param$ylab, main=param$main) 
  lines(t, nE[, 1], lty = 1, col=1, lwd=3)
  lines(t, nE[, 2], lty = 2, col=2, lwd=3)
  lines(t, nE[, 3], lty = 3, col=3, lwd=3)
  abline(h = seq(0, max(nE), 25), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  
  legend(0, max(nE[, 3]), leg$txt, col=1:3, lty=1:3, bty="n", cex=0.8)
  
}


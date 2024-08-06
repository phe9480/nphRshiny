#' Simulation of Two-Arm Trial Data Following Customized Survival Distribution with Non-uniform Accrual Pattern and Lost to Follow-up
#'
#' Simulate Randomized two-arm trial data with the following characteristics:  
#' (1) randomization time (entry time) is generated according to the specified non-uniform accrual pattern, 
#' i.e. the cumulative recruitment at calendar time t is (t/A)^w with weight w and enrollment complete in A months.
#' w = 1 means uniform enrollment, which is usually not realistic due to graduate sites activation process.
#' (2) Survival time follows a customized survival distribution for each arm.
#' (3) N total patients with r:1 randomization ratio
#' (4) Random drop off can be incorporated into the censoring process.
#' (5) Data cutoff dates are determined by specified vector of target events for all analyses.
#' (6) A dataset is generated for each analysis according to the specified number of target events. 
#'     Multiple analyses can be specified according to the vector of targetEvents, eg, targetEvents = c(100, 200, 300) 
#'     defines 3 analyses at 100, 200, and 300 events separately.
#'
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param S0 Survival function for control arm
#' @param S1 Survival function for experimental arm
#' @param drop Drop Off rate per month for each arm, eg, 3% every year for each arm, then drop=c(0.03/12, 0.03/12)
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @param DCO   A vector of data cut-off time in months, calculated from first subject in. 
#'              Default NULL. The date cut-off times will be determined by targetEvents.
#'              If provided, then the targetEvents will be ignored.
#'              
#' @param seed seed for random sampling              
#' 
#' @return An object with a dataframe for each analysis including the following variables:
#' \describe{
#' \item{sim}{sequence number of simulated dataset;}
#' \item{treatment}{treatment group with values of "control" and "experimental"}
#' \item{enterTime}{Time of randomization in calendar time}
#' \item{calendarTime}{the time when event/censoring occurred in calendar time}
#' \item{survTime}{Survival time for analysis, = calendarTime - enterTime}
#' \item{cnsr}{censor status (0=event; 1=censor) before administrative censoring due to data cut}
#' \item{calendarCutOff}{Data CutOff Time (DCO);}
#' \item{survTimeCut}{Survival time after cut}
#' \item{cnsrCut}{Censor status after cut}
#' }
#' @examples
#' #Example (1): Simulate 10 samples from proportional hazards scenario. 
#' #Total 600 pts, 1:1 randomization, control median OS 12 mo; 
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' m0 = 12; #median OS for control arm
#' lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' S1=function(t){exp(-lambda0*0.65*t)}
#' 
#' sim.ph = simulation.rand(nSim=10, N = 600, A = 21, w=1.5, r=1, S0 = S0, S1=S1, drop=c(0, 0), targetEvents = c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,50))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,50))
#' 
#' #Example (2): Simulate 10 samples with delayed effect at month 6;
#' #Total 600 pts, 1:1 randomization, control median OS 12 mo; 
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' HRd = 0.60 #hazard ratio after delay
#' 
#' h.D6=function(t){lambda0*as.numeric(t<6)+HRd*lambda0*as.numeric(t>=6)}
#' c6 = exp(-6*lambda0*(1-HRd)); 
#' S.D6 = function(t){S0(t)*as.numeric(t<6)+c6*exp(-HRd*lambda0*t)*as.numeric(t>=6)}
#' 
#' sim.D6 = simulation.rand(nSim=10, N = 600, A = 24, w=1.5, r=1, S0=S0, S1=S.D6, drop=c(0, 0), targetEvents = c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.D6[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.D6[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' 
#' @export 
simulation.rand = function(nSim=1, N = 600, A = 21, w=1.5, Lambda=NULL, r=1, S0=function(t){exp(-log(2)/12*t)}, S1=function(t){exp(-log(2)/12*0.65*t)}, drop=c(0,0), targetEvents = c(400, 500), DCO = NULL) {

  #number of analyses
  if(is.null(DCO)){
    L = length(targetEvents)
  } else{
    L = length(DCO)
  }
  
  nEachMonth = f.nEachMonth(N=N, A=A, w=w, r=r, Lambda=Lambda)
  n0 = sum(nEachMonth$n0); n1 = sum(nEachMonth$n1)
  CDF0 = function(t){1-S0(t)}
  CDF1 = function(t){1-S1(t)}
  
  out = list(NULL)
  if (L > 1){for (k in 2:L){out = c(out, list(NULL))}}
  
  ########################
  for (j in 1:nSim) {
    
    ########################
    #Sampling for each arm
    ########################
    T0 = rand(n=n0, CDF=CDF0)
    T1 = rand(n=n1, CDF=CDF1)
    
    #Permutation of the original ordered samples
    T0 = sample(T0); T1 = sample(T1)
    
    ############################
    #Drop Off data for each arm
    ############################
    if (is.na(drop[2])){drop = rep(drop[1], 2)}
    if (drop[1] > 0) {W0 = rexp(n0, rate=-log(1-drop[1]))} else {W0 = rep(Inf, n0)}
    if (drop[2] > 0) {W1 = rexp(n1, rate=-log(1-drop[2]))} else {W1 = rep(Inf, n1)}

    ############################
    #Censor data from Drop-off
    ############################
    Y0 = apply(cbind(T0, W0), 1, min)
    Y1 = apply(cbind(T1, W1), 1, min)
    
    event0 = as.numeric(T0 < Inf)
    event0[W0 < T0] = 0
    event1 = as.numeric(T1 < Inf)
    event1[W1 < T1] = 0
    
    ############################
    #EnterTime, CalendarTime
    ############################
    enterTime0 = rep(NA, n0)
    enterTime0[1:nEachMonth$n0[1]] = runif(nEachMonth$n0[1], min=0, max=1)
    if (A > 1) {for (m in 2:A){
      LL = sum(nEachMonth$n0[1:(m-1)])+1
      UU = sum(nEachMonth$n0[1:m])
      enterTime0[LL:UU] = runif(nEachMonth$n0[m], min=m-1, max=m)
    }}
    
    enterTime1 = rep(NA, n1)
    enterTime1[1:nEachMonth$n1[1]] = runif(nEachMonth$n1[1], min=0, max=1)
    if (A > 1) {for (m in 2:A){
      LL = sum(nEachMonth$n1[1:(m-1)])+1
      UU = sum(nEachMonth$n1[1:m])
      enterTime1[LL:UU] = runif(nEachMonth$n1[m], min=m-1, max=m)
    }}
    
    ############################
    #Assemble data before cut
    ############################
    sim = rep(j, N)
    treatment = c(rep("control", n0), rep("experimental", n1))
    enterTime = c(enterTime0, enterTime1)
    survTime = as.numeric(c(Y0, Y1))
    
    #trick infinity
    survTime[survTime > 1e6] = 1e6
    calendarTime = as.numeric(enterTime) + as.numeric(survTime)
    cnsr = c(1-event0, 1-event1)
    
    dati = data.frame(cbind(sim,enterTime, calendarTime, survTime, cnsr))
    dati$treatment = treatment
    
    ############################
    #Cut data
    ############################
    dati.cut = NULL
    for (ii in 1:L){
      dati.cut[[ii]] = f.dataCut(data=dati, targetEvents=targetEvents[ii], DCO = DCO[ii])
      out[[ii]] = rbind(out[[ii]], dati.cut[[ii]])
    }   
  }
  return(out)
}



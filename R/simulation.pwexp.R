#' Simulate data following piece-wise exponential distribution with non-uniform accrual pattern and lost to follow-up
#'
#' Simulate Randomized two-arm trial data with the following characteristics:  
#' (1) randomization time (entry time) is generated according to the specified non-uniform accrual pattern, 
#' i.e. the cumulative recruitment at calendar time t is (t/A)^w with weight w and enrollment complete in A months.
#' w = 1 means uniform enrollment, which is usually not realistic due to graduate sites activation process.
#' (2) Survival time follows piece-wise exponential distribution for each arm.
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
#' @param Lambda Cumulative distribution function (CDF) for enrollment on (0, infinity). For example, uniform enrollment of 20 patients / month for 24 months has Lambda = function(t){t/24*as.numeric(t<= 24) + as.numeric(t>24)}.   
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lam0 Hazard rates for control arm of intervals defined by cuts; for exponential(lam0) distribution,
#'         lam0 = log(2) / median;
#' @param  lam1 Hazard rates for experimental arm for intervals; for exponential(lam1) distribution,
#'         lam1 = log(2) / median; For delayed effect under H1, lam1 is a vector (below).
#' @param  cuts0 Timepoints to form intervals for piecewise exponential distribution for control arm. 
#' For example, cuts0 = c(6, 12) and lam0 = c(0.05, 0.07, 0.10) mean hazard 0.05 on (0, 6), 0.07 on (6, 12) and 0.1 on (12, infinity). 
#' length(cuts0) must be 1 less than length(lam0).
#' @param  cuts1 Timepoints to form intervals for piecewise exponential distribution for experimental arm. See cuts0.
#' @param drop0 Drop Off rate per month, eg, 3% every year for control arm, then drop0=0.03/12
#' @param drop1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              100 target events are used to determine IA DCO; and 150 events are used 
#'              to determine the FA DCO. Must be integers.
#' @param DCO   A vector of data cut-off time in months, calculated from first subject in. 
#'              Default NULL. The date cut-off times will be determined by targetEvents.
#'              If provided, then the targetEvents will be ignored.
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
#' sim.ph = simulation.pwexp(nSim=10, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, lam1= log(2)/12*0.65, cuts0=NULL, cuts1=NULL, drop0= 0, drop1= 0, targetEvents = c(400, 500))
#' #same as above
#' sim.ph = simulation.pwexp(nSim=10, N = 600, Lambda= function(t){(t/24)^1.5*as.numeric(t<= 24) + as.numeric(t>24)}, r=1, lam0=log(2)/12, lam1= log(2)/12*0.65, cuts0=NULL, cuts1=NULL, drop0= 0, drop1= 0, targetEvents = c(400, 500))
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
#' sim.delay6 = simulation.pwexp(nSim=10, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, lam1= c(log(2)/12,log(2)/12*0.65), cuts0=NULL, cuts1=6, drop0= 0, drop1= 0, targetEvents = c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Example (3): Simulate 10 samples with delayed effect at month 6 
#' #Control arm has crossover to subsequent IO after 24 mo, so its hazard decreases 20%.
#' #control arm has constant hazard (median 11.7 mo) and experimental arm has 
#' #hr = 1 and 0.65 at intervals (0, 6) and (6, 24) respectively.
#' #HR = 0.65, enrollment 24 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 400 and 500 events respectively.
#' 
#' lam0 = log(2)/12*c(1, 1, 0.8); lam1 = log(2)/12*c(1, 0.65, 0.65)
#' sim.delay6crs=simulation.pwexp(nSim=100,N=600,A=24,w=1.5,r=1,lam0=lam0, lam1=lam1,cuts0=c(6, 24),cuts1=c(6, 24),drop0=0,drop1=0, targetEvents=c(400, 500))
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6crs[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6crs[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' @export 
simulation.pwexp = function(nSim=100, N = 600, A = 21, w=1.5, Lambda = NULL, 
                            r=1, lam0=log(2)/12, lam1=log(2)/12*0.65, 
                            cuts0=NULL, cuts1=NULL, drop0=0, drop1=0, 
                            targetEvents = c(400, 500), DCO = NULL) {

  nEachMonth = nphRshiny::f.nEachMonth(N=N, A=A, w=w, r=r, Lambda=Lambda)
  
  gamma = nEachMonth$n0 + nEachMonth$n1
  eta0 = -log(1-drop0)
  eta1 = -log(1-drop1)
  
  #G0 Cumulative distribution function of drop-off for control arm. 
  #For example, 3 percent drop-off in 12 months of followup means
  #the hazard rate per month is 
  #eta0 = -log(1-0.03/12), so G0=function(t){1-exp(-eta0*t)}.
  
  #Determine intervals, lambdaC, lambdaE
  if (is.null(cuts0) && is.null(cuts1)){
    intervals = NULL
    lam0new = lam0
    lam1new = lam1
  }
  if (is.null(cuts0) && !is.null(cuts1)){
    intervals = rep(NA, length(cuts1)); intervals[1] = cuts1[1]
    if (length(cuts1) > 1){
    for (i in 2:length(cuts1)){
      intervals[i] = cuts1[i] - cuts1[i-1]
    }}
    lam0new = rep(lam0[1], length(cuts1))
    lam1new = lam1
  }
  if (!is.null(cuts0) && is.null(cuts1)){
    intervals = rep(NA, length(cuts0)); intervals[1] = cuts0[1]
    if (length(cuts0) > 1){
    for (i in 2:length(cuts0)){
      intervals[i] = cuts0[i] - cuts0[i-1]
    }}
    lam0new = lam0
    lam1new = rep(lam1[1], length(cuts0))
  }
  if (!is.null(cuts0) && !is.null(cuts1)){
    cuts = sort(unique(c(cuts0, cuts1)))
    intervals = rep(NA, length(cuts)); intervals[1] = cuts[1]
    if (length(cuts) > 1){
    for (i in 2:length(cuts)){
      intervals[i] = cuts[i] - cuts[i-1]
    }}
  
    lam0new = lam1new = rep(NA, (length(cuts)+1))
    for (i in 1:length(cuts)){
      for (j in 1:length(cuts0)){
        if (cuts0[j] == cuts[i]) {lam0new[i] = lam0[j]}
      }
      for (j in 1:length(cuts1)){
        if (cuts1[j] == cuts[i]) {lam1new[i] = lam1[j]}
      }
    }
    for (i in length(lam0new):1){ if (!is.na(lam0new[i])){ix0 = i; break;}}
    for (i in ix0:1){
      if (i > 1 && !is.na(lam0new[i]) && is.na(lam0new[i-1])){lam0new[i-1] = lam0new[i]}
    }
    if (ix0 < length(lam0new)) {
      for (i in (ix0+1):length(lam0new)){lam0new[i] = lam0[length(lam0)]}
    }
  for (i in length(lam1new):1){ if (!is.na(lam1new[i])){ix1 = i; break;}}
  for (i in ix1:1){
    if (i > 1 && !is.na(lam1new[i]) && is.na(lam1new[i-1])){lam1new[i-1] = lam1new[i]}
  }
  if (ix1 < length(lam1new)) {
    for (i in (ix1+1):length(lam1new)){lam1new[i] = lam1[length(lam1)]}
  }
  }
  o = nphsim::nphsim(nsim=nSim,lambdaC=lam0new,lambdaE=lam1new, ssC=N/(r+1), ssE=N*r/(r+1),
             intervals=intervals, gamma=gamma, R=rep(1, A),eta=eta0, etaE=eta1,fixEnrollTime = FALSE)
  dat = o$simd
  data.out = NULL
  #number of analyses
  if(is.null(DCO)){
    L = length(targetEvents)
  } else{
    L = length(DCO)
  }
  
  for (j in 1:nSim) {
    
    dataj = dat[dat$sim == j,]
    
    #rename variables
    dataj <- dplyr::rename(dataj, enterTime = enterT, calendarTime = ct, survTime= survival)
    dataj$group = ifelse(dataj$treatment == "control", 0, 1)
    
    #cut data according to the specified target events    
    dataj.cut <- lapply(1:L, function(k) {
      f.dataCut(data=dataj, targetEvents[k], DCO = DCO[k])
    })
    
    data.out=lapply(1:L, function(k) {
      rbind(data.out[[k]], dataj.cut[[k]])
    })
    
  }
  return(data.out)
}


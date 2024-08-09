#' Trial Data Simulations and Analysis Using Weighted Log-rank Test with non-uniform accrual pattern and lost to follow-up
#'
#' Simulate Randomized two-arm trial data with the following characteristics:  
#' (1) randomization time (entry time) is generated according to the specified non-uniform accrual pattern, 
#' i.e. the cumulative recruitment at calendar time t is (t/A)^w with weight w and enrollment complete in A months.
#' w = 1 means uniform enrollment, which is usually not realistic due to graduate sites activation process. It also allows customized enrollment pattern by specifying the cumulative enrollment distribution with domain (0, infinity).
#' (2) Survival time follows commonly used distributions and also allow customized survival functions. 
#' (3) Allow different drop off rates for both arms.
#' (4) Data cutoff dates are determined by specified vector of target events for all analyses.
#' (5) A dataset is generated for each analysis according to the specified number of target events. 
#'     Multiple analyses can be specified according to the vector of targetEvents, eg, targetEvents = c(100, 200, 300) 
#'     defines 3 analyses at 100, 200, and 300 events separately.
#' (6) Weighted log-rank test is then performed for each simulated group sequential dataset.
#'
#'
#' @param  nSim Number of trials
#' @param  n Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param Lambda Cumulative distribution function (CDF) for enrollment on (0, infinity). 
#' For example, uniform enrollment of 20 patients / month for 24 months has 
#' Lambda = function(t){t/24*as.numeric(t<= 24) + as.numeric(t>24)}. When Lambda is specified, A and w are ignored.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param dist0 Type of distribution for control arm. The following options are available. 
#' (1) "exponential": lam0 required;
#' (2) "weibull": shape0 and scale0 required; 
#' (3) "piecewise exponential": lamb0 and cuts0 required;
#' (4) "mixture cure rate of exponential": p10 and lam0 required; 
#' (5) "mixture cure rate of weibull": p10, shape0, and scale0 required
#' (6) "customized": S0 as a survival function is required.
#' @param  lam0 Hazard rates for control arm for intervals; for exponential distribution,
#'         lam0 = log(2) / median; For piecewise exponential distribution, lam0 is a vector of length 1 greater than length of cuts0.
#' @param shape0 shape parameter for weibull distribution for control arm. Refer to rweibull() for details.
#' @param scale0 scale parameter for weibull distribution for control arm. Refer to rweibull() for details. 
#' @param p10 cure rate parameter for mixture cure rate distribution for control arm
#' @param S0 Survival function for customized distribution for control arm.
#' @param cuts0 Cut points for piecewise exponential distribution for control arm
#' @param dist1 Type of distribution for experimental arm. See dist0. 
#' @param lam1 hazard rate for piecewise exponential distribution for experimental arm
#' @param shape1 shape parameter for weibull distribution for experimental arm. Refer to rweibull() for details.
#' @param scale1 scale parameter for weibull distribution for experimental arm. Refer to rweibull() for details. 
#' @param p11 cure rate parameter for mixture cure rate distribution for experimental arm
#' @param S1 Survival function for customized distribution for control arm.
#' @param cuts1 Cut points for piecewise exponential distribution for experimental arm
#' @param  lam1 Hazard rates for experimental arm for intervals; for exponential distribution,
#'         lam1 = log(2) / median; For piecewise exponential distribution, lam1 is a vector of length 1 greater than length of cuts1.
#' @param drop0 Drop Off rate per month, eg, 1%, for control arm
#' @param drop1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              100 target events are used to determine IA DCO; and 150 events are used 
#'              to determine the FA DCO. Must be integers.
#' @param  overall.alpha  Allocated overall alpha (one-sided) for group sequential design
#' @param  sf Spending function. LanDeMets O'Brien Fleming: "LDOF", LanDeMets Pocock: "LDPK", "HSD": Hwang-Shih-DeCani spending function with parameter param.
#' @param param parameter for Hwang-Shih-DeCani spending function
#' @param p1 A fixed p value boundary for IAs (one-sided), which is applicable to Haybittle-Peto alpha spending only.
#' @param cum.alpha Cumulative alpha spending by analysis, which is applicable to Bespoke method only. Cum.alpha must have the same length as timing.
#' @param alpha Allocated one-sided alpha levels. sum(alpha) is the total type I error.
#'           If alpha spending function a(t) is used for information time c(t1, ..., tK),
#'           then alpha1 = a(t1), alpha2 = a(t2)-a(t1), ..., alphaK = a(tK)-a(t_{K-1}),
#'           and the total alpha for all analyses is a(tK). When alpha is provided, sf is ignored.
#' @param logrank Indicator whether log-rank test is requested besides the weighted logrank tests. "Y" or "N". Default "Y".
#'                 If "Y", the traditional log-rank test will be used based on survdiff() function.
#'                 If "N", the weighted log-rank test with weighting function specified in fws will be used.
#' @param fws.options   Weighting strategies in the following format as examples. If fws.options is provided,
#' then the nphDesign object's weighting strategy will be ignored. fws can contain multiple weighting strategies.
#'        For example, fws = list(fws1, fws2) means 2 weighting strategies are evaluated, where
#'        fws1 = list(IA = list(lr), FA=list(lr, fh01)); fws2 = list(IA = list(lr), FA=list(fh01)). Each
#'        weighting strategy is specified as following examples for illustration.
#'        \itemize{
#'        \item (1) Single-time analysis using log-rank test: fws1 = list(FA = list(lr));
#'        \item (2) Two interim analyses and 1 final analysis with logrank at IA1, 
#'        max(logrank, fleming-harrington(0,1) at IA2, and 
#'        max(logrank, fleming-harrington(0,1), fleming-harrington(1,1)) at final): 
#'        fws2 = list(IA1 = list(lr), IA2=list(lr, fh01), FA=list(lr,fh01,fh11)).
#'        \item (3) One IA and one FA: stabilized Fleming-Harrington (0,1) at IA, 
#'        and max(logrank, stabilized Fleming-Harrington (0, 1)) at FA. 
#'        fw3 = list(IA = list(sfh01), FA=list(lr, sfh01)).
#'        \item General format of weighting strategy specification is: (a) the weight
#'        functions for each analysis must be provided in list() object even there is only
#'        1 weight function for that an analysis. (b) The commonly used functions 
#'        are directly available including lr: log-rank test; fh01: Fleming-Harrington (0,1) test;
#'        fh11: Fleming-Harrington(1,1) test; fh55: Fleming-Harrington(0.5, 0.5) test.
#'        stabilized versions at median survival time: sfh01, sfh11, sfh55. Modestly 
#'        log-rank test: mlr.
#'        (c) User-defined weight function can also be handled, but the weight function
#'        must be defined as a function of survival rate, i.e., fws = function(s){...}, 
#'        where s is the survival rate S(t-). 
#'        \item Options of weighting strategies for exploration: fws.options=list(fws1, fws2, fws3).
#'        }
#' @param H0 "Y" or "N" to indicate whether the simulation is for type I error
#' @param parallel True/False indicator. If true, use parallel computing weighted log-rank test for each analysis in strategy m
#' @param n.cores This will be used if parallelization is TRUE. Default is 8.
#' @param seed seed for generating samples. Default 2022
#' @param out.z Output test statistics, TRUE/FALSE
#'                  
#' @return An object with a dataframe for each analysis including the following variables:
#' \describe{
#' \item{power}{Power for each analysis}
#' \item{overall.power}{Overall power of the group sequential design}
#' \item{wlr.simulations}{Simulation results for each simulated study data. 
#' An array with dimensions: nSim simulations X M testing strategies 
#' (fws.options) X K analyses X variables:
#'      \itemize{ 
#'      \item z value
#'      \item p value
#'      \item analysis
#'      \item rejection boundary
#'      \item testing result (1 = positive; 0 = negative)
#'      }}
#' \item{lr.power}{Power for each analysis using log-rank test. Available if logrank ="Y"}
#' \item{lr.overall.power}{Overall power of the group sequential design using logrank test}
#' \item{lr.simulations}{Simulation results for each simulated study data using logrank test}
#' }
#' @examples
#' 
#' #Utility functions
#' lr = nphRshiny:::lr; fh01 = nphRshiny:::fh01; fh11 = nphRshiny:::fh11
#' fws1 = list(IA1 = list(lr), FA = list(lr))
#' fws2 = list(IA1 = list(lr), FA = list(fh01))
#' fws3 = list(IA1 = list(fh01), FA = list(fh01))
#' fws4 = list(IA1 = list(lr), FA = list(lr, fh01))
#' fws5 = list(IA1 = list(lr), FA = list(lr, fh01, fh11))
#' fws6 = list(IA1 = list(lr, fh01), FA = list(lr, fh01))
#' 
#' #Weighted logrank tests options
#' fws = list(fws1, fws2, fws3, fws4, fws5, fws6)
#' 
#' #Example (1): Simulate 10 samples from proportional HR with HR = 0.7
#' 
#' #Hazard and survival distributions
#' 
#' #Control Arm
#' m0 = 10; lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' 
#' #Experimental Arm
#' HR = 0.7; h1 = function(t){lambda0*HR}; 
#' S1 = function(t){exp(-lambda0 *HR* t)}
#' 
#' #Enrollment
#' F.entry = function(t){(t/21)^1.5*as.numeric(t <= 21) + as.numeric(t > 21)}
#' 
#' #Drop-off
#' drop0 = drop1 = 0.03/12
#' 
#' e = rep(NA, 2); DCO=c(24, 36)
#' eta0 = -log(1-drop0); G0=function(t){1-exp(-eta0*t)}
#' eta1 = -log(1-drop1); G1=function(t){1-exp(-eta1*t)}
#' for (i in 1:length(DCO)){e[i] = fe(DCO = DCO[i], r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
#' Lambda = F.entry, n = 100, G0=G0, G1=G1)$e}
#' 
#' #(a) Study design using weighted logrank test option 1
#' wlr.power.maxcombo(DCO = c(24, 36), overall.alpha=0.025, sf = "LDOF", 
#'   r = 1, n = 100, h0 = h0, S0=S0, h1 = h1, S1= S1, 
#'   f.ws = fws1, Lambda=F.entry, G0=G0, G1=G1, 
#'   mu.method = "Schoenfeld", cov.method = "H0")
#'   
#' #(b) Study design using weighting strategy 2
#' wlr.power.maxcombo(DCO = c(24, 36), overall.alpha=0.025, sf = "LDOF", 
#'   r = 1, n = 100, h0 = h0, S0=S0, h1 = h1, S1= S1, 
#'   f.ws = fws2, Lambda=F.entry, G0=G0, G1=G1, 
#'   mu.method = "Schoenfeld", cov.method = "H0")
#'   
#' #(c) Simulations for exploring weighted logrank tests option 1 and 2
#' 
#' #(c1) Type I error; also output logrank test simulation results
#' H0 = simulation.nphDesign.pwexp(nSim=5, N = 100, r = 1, 
#' A = 21, w=1.5,
#' lam0=lambda0, lam1=lambda0*0.7,
#' targetEvents = e, drop0 = 0.03/12, drop1=0.03/12,
#' overall.alpha = 0.025, sf = "LDOF",
#' H0 = "Y", logrank="Y", fws.options=list(fws1))
#' 
#' #same as above; using F.entry function to replac A and w specifications.
#' H0 = simulation.nphDesign.pwexp(nSim=5, N = 100, r=1, 
#' Lambda=F.entry, 
#' lam0=lambda0, lam1=lambda0*0.7,
#' targetEvents = e, drop0 = 0.03/12, drop1=0.03/12,
#' overall.alpha = 0.025, sf = "LDOF",
#' H0 = "Y", logrank="Y", fws.options=list(fws1))
#' 
#' #same as above using the general function
#' simulation.nphDesign(nSim=5, n = 100, r=1, Lambda=F.entry, drop0=0.03/12, drop1=0.03/12, 
#' dist0 = "exponential", lam0=lambda0, shape0 = NULL, scale0 = NULL, p10 = NULL, S0 = NULL, cuts0 = NULL,
#' dist1 = "exponential", lam1=lambda0, shape1 = NULL, scale1 = NULL, p11 = NULL, S1 = NULL, cuts1=NULL, 
#' targetEvents = e, sf = "LDOF", param = NULL, overall.alpha = 0.025, p1=NULL, cum.alpha=NULL,
#' logrank="N", fws.options=list(fws1), 
#' parallel=FALSE, n.cores=8, seed=2022, out.z = FALSE)
#' 
#' 
#' @export 
simulation.nphDesign = function(nSim=3, n = 100, r=1, A = 21, w=1.5, 
                                Lambda=NULL, drop0=0, drop1=0, 
                                dist0 = "exponential", lam0=log(2)/12,
                                shape0 = NULL, scale0 = NULL,
                                p10 = NULL, S0 = NULL, cuts0 = NULL,
                                dist1 = "exponential", lam1=log(2)/12*0.7,
                                shape1 = NULL, scale1 = NULL,
                                p11 = NULL, S1 = NULL, cuts1=NULL, 
                                targetEvents = c(30, 60), 
                                sf = "LDOF", param = NULL, overall.alpha = 0.025, p1=NULL, cum.alpha=NULL,
                                logrank="N", fws.options=list(fws5), 
                                parallel=FALSE, n.cores=8, seed=2022, out.z = FALSE) {

    set.seed(seed)
    side = 1 #always one-sided
    
    #Entry data
    nEachMonth = f.nEachMonth(N=n, A=A, w=w, r=r, Lambda=Lambda)
    n0 = sum(nEachMonth$n0); n1 = sum(nEachMonth$n1)
    
  ##############################
  #M Options of test strategies
  M = length(fws.options)
  
  #K analyses
  K=length(fws.options[[1]])
  targetEvents = ceiling(targetEvents)
  
  timing = targetEvents / targetEvents[K]
  
  #incremental alpha
  alpha = f.alpha(overall.alpha=overall.alpha, sf=sf, timing=timing, p1=p1, cum.alpha=cum.alpha, param=param)

  #weighted logrank test statistics
  wlr.sim = array(NA, dim=c(nSim, M, K, 5))
  if (logrank == "Y") {lr.sim = array(NA, dim=c(nSim, K, 5))}
  
  #output datasets
  dati.out = list(NULL)
  if (K > 1){for (k in 2:K){dati.out = c(dati.out, list(NULL))}}
  
  if(parallel){
    ##############
    #Run wlr test in parallel
    ##############
    library(parallel)   #For PSOCK parallel backend
    library(foreach)    #For 'foreach' loop
    library(doParallel) #For using PSOCK cluster with operator '%dopar%'
    n.cores <- n.cores
    my.cluster <- makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    clusterEvalQ(my.cluster, .libPaths(.libPaths())
    )
    registerDoParallel(cl = my.cluster)
    clusterEvalQ(my.cluster, library(nphDesign))
    sim.data=NULL
    
    for (i in 1:nSim) {
      #(1). Generate data
      dati = simulation.pwexp(nSim=1, N = n, A = A, w=w, Lambda=Lambda, r=r, lam0=lam0, lam1=lam1, 
                       cuts=cuts, drop0=drop0, drop1=drop1, targetEvents = targetEvents)
 
      sim.data[[i]]<-dati
    }
    
    #(3). Testing strategy m
    for (m in 1:M){    
      #Perform weighted log-rank test for each analysis in strategy m
      ## This next line takes a long-time to compute, making simulation difficult. Use parallel computing instead
      wlr <- foreach(i = 1:nSim) %dopar% {
        wlr.inference(data=sim.data[[i]], alpha = alpha, f.ws=fws.options[[m]])$test.results
      }
      stopCluster(cl = my.cluster) ## Stop cluster
      for (i in 1:nSim) {
        wlri = wlr[[i]][!duplicated(wlr[[i]]$analysis),]
        wlri$result = as.numeric(wlri$inference=="Positive")
        wlr.sim[i, m, , ] = as.matrix(wlri[,c(2,3,5,6,8)])
      }
    }    
  }else{
    for (i in 1:nSim) {
      #(1). Generate survival data for each arm
      T0 = genSurv(dist = dist0, n = n0, lam=lam0, shape=shape0, scale=scale0,
                   p1=p10, S=S0, cuts=cuts0)
      T1 = genSurv(dist = dist1, n = n1, lam=lam1, shape=shape1, scale=scale1,
                   p1=p11, S=S1, cuts=cuts1)
      
      #Permutation of the original ordered samples
      T0 = sample(T0); T1 = sample(T1)
      
      #(2). Drop Off data for each arm
      ############################
      if (drop0 > 0) {W0 = rexp(n0, rate=-log(1-drop0))} else {W0 = rep(Inf, n0)}
      if (drop1 > 0) {W1 = rexp(n1, rate=-log(1-drop1))} else {W1 = rep(Inf, n1)}
      
      #(3). Censor data from Drop-off
      ############################
      Y0 = apply(cbind(T0, W0), 1, min)
      Y1 = apply(cbind(T1, W1), 1, min)
      
      event0 = as.numeric(T0 < Inf)
      event0[W0 < T0] = 0
      event1 = as.numeric(T1 < Inf)
      event1[W1 < T1] = 0
      
      #(4). EnterTime, CalendarTime
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
      
      #(5). Assemble data before cut
      ############################
      sim = rep(i, n)
      treatment = c(rep("control", n0), rep("experimental", n1))
      enterTime = c(enterTime0, enterTime1)
      survTime = as.numeric(c(Y0, Y1))
      
      #trick infinity
      survTime[survTime > 1e6] = 1e6
      calendarTime = as.numeric(enterTime) + as.numeric(survTime)
      cnsr = c(1-event0, 1-event1)
      
      dati = data.frame(cbind(sim, enterTime, calendarTime, survTime, cnsr))
      dati$treatment = treatment
      dati$group = ifelse(dati$treatment == "control", 0, 1)
      
      #(6). Cut data
      ############################
      dati.cut = NULL
      for (ii in 1:K){
        dati.cut[[ii]] = f.dataCut(data=dati, targetEvents=targetEvents[ii])
        dati.out[[ii]] = rbind(dati.out[[ii]], dati.cut[[ii]])
      }
      
      #dati = simulation.pwexp(nSim=1, N = n, A = A, w=w, Lambda=Lambda, r=r, 
      #                        lam0=lam0, lam1=lam1, cuts0=cuts0, cuts1=cuts1,
      #                        drop0=drop0, drop1=drop1, targetEvents = targetEvents)
      
      #(7). Testing strategy m
      for (m in 1:M){    
        #Perform weighted log-rank test for each analysis in strategy m
        wlri = wlr.inference(data=dati.out, alpha = alpha, f.ws=fws.options[[m]])$test.results
        wlri = wlri[!duplicated(wlri$analysis),]
        wlri$result = as.numeric(wlri$inference=="Positive")
        wlr.sim[i, m, , ] = as.matrix(wlri[,c(2,3,5,6,8)])
      }
      
      #(4). Standard log-rank test for all analyses if requested
      if (logrank=="Y"){
        if(K > 1){
          #GSD boundary for each analysis
          if(side == 1) {
            z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha,timing=timing[1:(K-1)], 
                                       sfu=gsDesign::sfLDOF)$upper$bound 
          } else{
            z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha/2,timing=timing[1:(K-1)], 
                                       sfu=gsDesign::sfLDOF)$upper$bound 
          }
        } else {
          if(side == 1) {z.bd = qnorm(1-overall.alpha)} else{z.bd = qnorm(1-overall.alpha/2)}
        }
        for (j in 1:K){
          lr.test = survival::survdiff(survival::Surv(survTimeCut, 1-cnsrCut) ~ group, data = dati[[j]])
          
          #convert to z value in correct direction: z>0 means better experimental arm.
          better = as.numeric(lr.test$obs[1] > lr.test$obs[2])
          sign = 2*better - 1
          z = sqrt(lr.test$chisq) * sign
          
          #count power
          lr.sim[i, j, 1] = z
          if (side == 1){ p = 1 - pnorm(z)} else {p = 2*(1 - pnorm(z))}
          
          lr.sim[i, j, 2] = p
          lr.sim[i, j, 3] = j
          lr.sim[i, j, 4] = z.bd[j]
          lr.sim[i, j, 5] = as.numeric(z > z.bd[j])
        }
      }
      
    }
  }
  
  pow = matrix(NA, nrow=M, ncol=K)
  overall.pow = rep(0, M)
  
  for(m in 1:M){for(j in 1:K){pow[m, j] = sum(wlr.sim[,m,j,5])/nSim}}
  for(m in 1:M){for(i in 1:nSim){
    overall.pow[m] =  overall.pow[m] + as.numeric(sum(wlr.sim[i,m,,5])>0)
  }}
  overall.pow = overall.pow / nSim
  
  o=list()
  o$power = pow; o$overall.power = overall.pow
  if (out.z) {o$wlr.simulations = wlr.sim}
  
  if(logrank=="Y"){
    lr.pow = rep(NA, K)
    for (j in 1:K) {lr.pow[j] = sum(lr.sim[,j,5])/nSim}

    lr.overall.pow = 0
    for (i in 1:nSim){lr.overall.pow = lr.overall.pow + as.numeric(sum(lr.sim[i,,5])>0)}
    lr.overall.pow = lr.overall.pow / nSim
    
    o$lr.overall.power = lr.overall.pow
    o$lr.power = lr.pow
    if (out.z) {o$lr.simulations = lr.sim}
  }
  return(o)
}



#' Trial Data Simulations and Analysis Using Weighted Log-rank Test For Piecewise Exponential Distribution
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
#' (7) Weighted log-rank test is then performed for each simulated group sequential dataset.
#'
#'
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param Lambda Cumulative distribution function (CDF) for enrollment on (0, infinity). 
#' For example, uniform enrollment of 20 patients / month for 24 months has 
#' Lambda = function(t){t/24*as.numeric(t<= 24) + as.numeric(t>24)}. When Lambda is specified, A and w are ignored.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lam0 Hazard rates for control arm of intervals defined by cuts; for exponential(lambda0) distribution,
#'         lam0 = log(2) / median;
#' @param  lam1 Hazard rates for experimental arm for intervals; for exponential(lambda1) distribution,
#'         lam1 = log(2) / median; For delayed effect under H1, lambda1 is a vector (below).
#' @param  cuts Timepoints to form intervals for piecewise exponential distribution. For example,
#   \itemize{
#   \item Proportional hazards with hr = 0.65. Then lambda0 = log(2)/m0, lambda1 = log(2)/m0*hr, cuts = NULL. 
#   \item Delayed effect at month 6, and control arm has constant hazard (median m0) and 
#'       experimental arm has hr = 0.6 after delay, then cuts = 6, and 
#'       lam0 = log(2) / m0 or lam0 = rep(log(2) / m0, 2), 
#'       lam1 = c(log(2)/m0, log(2)/m0*hr). 
#   \item Delayed effect at month 6, and control arm has crossover to subsequent IO 
#'       treatment after 24 mo, so its hazard decreases 20%. Then, 
#'       lam0 = c(log(2)/m0, log(2)/m0, log(2)/m0*0.8), 
#'       lam1 = c(log(2)/m0, log(2)/m0*hr, log(2)/m0*hr), and
#'       cuts = c(6, 24), which forms 3 intervals (0, 6), (6, 24), (24, infinity)
#       }
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
#' #Example (1): Simulate 10 samples from proportional HR
#' 
#' #Hazard and survival distributions
#' 
#' #Control Arm
#' m0 = 10; lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' 
#' #Experimental Arm
#' h1 = function(t){lambda0*0.7}; 
#' S1 = function(t){exp(-lambda0 *0.7* t)}
#' 
#' #Enrollment
#' F.entry = function(t){(t/21)^1.5*as.numeric(t <= 21) + as.numeric(t > 21)}
#' 
#' #Drop-off
#' eta0 = -log(1-0.03/12) #control arm monthly drop off rate 0.03/12.
#' eta1 = -log(1-0.03/12) #control arm monthly drop off rate 0.03/12.
#' G0 = function(t){1-exp(-eta0*t)}; G1 = function(t){1-exp(-eta1*t)}; 
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
#' e = rep(NA, 2); DCO=c(24, 36)
#' for (i in 1:length(DCO)){e[i] = fe(DCO = DCO[i], r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
#' Lambda = F.entry, n = 100, G0=G0, G1=G1)$e}
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
#' @export 
simulation.nphDesign.pwexp = function(nSim=10000, N = 100, A = 21, w=1.5, Lambda=NULL, r=1, lam0=log(2)/12, lam1=log(2)/12*0.7, 
    cuts0=NULL, cuts1=NULL, drop0=0, drop1=0, targetEvents = c(30, 60), 
    sf = "LDOF", param = NULL, overall.alpha = 0.025, p1=NULL, cum.alpha=NULL, alpha = NULL,
    logrank="N", fws.options=list(fws5), H0 = "N", parallel=FALSE, n.cores=8, seed=2022, out.z = FALSE) {

    set.seed(seed)
    side = 1 #always one-sided
  #Simulation for checking type I error
  if (H0 == "Y"){lam1 = lam0} 
  
  ##############################
  #M Options of test strategies
  M = length(fws.options)
  
  #K analyses
  K=length(fws.options[[1]])
  targetEvents = ceiling(targetEvents)
  
  timing = targetEvents / targetEvents[K]
  
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(alpha) && !is.null(overall.alpha)){
    alpha = f.alpha(overall.alpha=overall.alpha, sf=sf, timing=timing, p1=p1, cum.alpha=cum.alpha, param=param)
  }
  
  wlr.sim = array(NA, dim=c(nSim, M, K, 5))
  if (logrank == "Y") {lr.sim = array(NA, dim=c(nSim, K, 5))}
  
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
      dati = simulation.pwexp(nSim=1, N = N, A = A, w=w, Lambda=Lambda, r=r, lam0=lam0, lam1=lam1, 
                       cuts0=cuts0, cuts1=cuts1, drop0=drop0, drop1=drop1, targetEvents = targetEvents)
 
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
      #(1). Generate data
      dati = simulation.pwexp(nSim=1, N = N, A = A, w=w, Lambda=Lambda, r=r, 
                              lam0=lam0, lam1=lam1, 
                       cuts0=cuts0, cuts1=cuts1, drop0=drop0, drop1=drop1, targetEvents = targetEvents)
      
      #(3). Testing strategy m
      for (m in 1:M){    
        #Perform weighted log-rank test for each analysis in strategy m
        wlri = wlr.inference(data=dati, alpha = alpha, f.ws=fws.options[[m]])$test.results
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
          better = as.numeric(lr.test$obs[2] < lr.test$exp[2])
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
  
  pow = cum.pow = matrix(0, nrow=M, ncol=K)
  overall.pow = rep(0, M)
  
  for(m in 1:M){for(j in 1:K){
    pow[m, j] = sum(wlr.sim[,m,j,5])/nSim
    for (i in 1:nSim){cum.pow[m,j]=cum.pow[m,j]+as.numeric(sum(wlr.sim[i,m,1:j,5])>=1)}
    cum.pow[m,j]=cum.pow[m,j]/nSim
    
    }}
  for(m in 1:M){for(i in 1:nSim){
    overall.pow[m] =  overall.pow[m] + as.numeric(sum(wlr.sim[i,m,,5])>0)
  }}
  overall.pow = overall.pow / nSim

  o=list()
  o$power = pow; o$overall.power = overall.pow
  o$cum.pow = cum.pow
  
  if (out.z) {o$wlr.simulations = wlr.sim}
  
  if(logrank=="Y"){
    lr.pow = lr.cum.pow = rep(0, K)
    for (j in 1:K) {
      lr.pow[j] = sum(lr.sim[,j,5])/nSim
      for (i in 1:nSim){lr.cum.pow[j]=lr.cum.pow[j]+as.numeric(sum(lr.sim[i,1:j,5])>=1)}
      lr.cum.pow[j]=lr.cum.pow[j]/nSim
    }

    lr.overall.pow = 0
    for (i in 1:nSim){lr.overall.pow = lr.overall.pow + as.numeric(sum(lr.sim[i,,5])>0)}
    lr.overall.pow = lr.overall.pow / nSim
    
    o$lr.overall.power = lr.overall.pow
    o$lr.power = lr.pow
    o$lr.cum.pow = lr.cum.pow
    if (out.z) {o$lr.simulations = lr.sim}
  }
  return(o)
}



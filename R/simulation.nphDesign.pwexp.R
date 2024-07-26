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
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lambda0 Hazard rates for control arm of intervals defined by cuts; for exponential(lambda0) distribution,
#'         lambda0 = log(2) / median;
#' @param  lambda1 Hazard rates for experimental arm for intervals; for exponential(lambda1) distribution,
#'         lambda1 = log(2) / median; For delayed effect under H1, lambda1 is a vector (below).
#' @param  cuts Timepoints to form intervals for piecewise exponential distribution. For example,
#   \itemize{
#   \item Proportional hazards with hr = 0.65. Then lambda0 = log(2)/m0, lambda1 = log(2)/m0*hr, cuts = NULL. 
#   \item Delayed effect at month 6, and control arm has constant hazard (median m0) and 
#'       experimental arm has hr = 0.6 after delay, then cuts = 6, and 
#'       lamda0 = log(2) / m0 or lambda0 = rep(log(2) / m0, 2), 
#'       lamda1 = c(log(2)/m0, log(2)/m0*hr). 
#   \item Delayed effect at month 6, and control arm has crossover to subsequent IO 
#'       treatment after 24 mo, so its hazard decreases 20%. Then, 
#'       lambda0 = c(log(2)/m0, log(2)/m0, log(2)/m0*0.8), 
#'       lambda1 = c(log(2)/m0, log(2)/m0*hr, log(2)/m0*hr), and
#'       cuts = c(6, 24), which forms 3 intervals (0, 6), (6, 24), (24, infinity)
#       }
#' @param dropOff0 Drop Off rate per month, eg, 1%, for control arm
#' @param dropOff1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @param logrank Indicator whether log-rank test is requested besides the weighted logrank tests. "Y" or "N". Default "Y".
#'                 If "Y", the traditional log-rank test will be used based on survdiff() function.
#'                 If "N", the weighted log-rank test with weighting function specified in fws will be used.
#' @fws.options   Weighting strategies in the following format as examples. If fws.options is provided,
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
#' #Example (1): Simulate 10 samples from proportional hazards scenario. 
#' fws1 = list(IA1 = list(lr), FA = list(lr))
#' fws2 = list(IA1 = list(lr), FA = list(fh01))
#' fws3 = list(IA1 = list(fh01), FA = list(fh01))
#' fws4 = list(IA1 = list(lr), FA = list(lr, fh01))
#' fws5 = list(IA1 = list(lr), FA = list(lr, fh01, fh11))
#' fws6 = list(IA1 = list(lr, fh01), FA = list(lr, fh01))
#' fws7 = list(IA1 = list(sfh01), FA = list(lr, fh01))
#' 
#' #7 weighting strategies for exploration 
#' fws = list(fws1, fws2, fws3, fws4, fws5, fws6, fws7)
#' 
#' alpha = f.alpha(overall.alpha=0.025, side = 1, sf="LDOF", targetEvents=c(60, 80))
#' m0 = 11.7
#' lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' h1 = function(t){lambda0*0.745}; 
#' S1 = function(t){exp(-lambda0 *0.745* t)}
#' 
#' F.entry = function(t){(t/21)^1.5*as.numeric(t <= 21) + as.numeric(t > 21)}
#' G.ltfu = function(t){0}
#' f.logHR = function(t){log(0.745)}
#' 
#' #study design using weighting strategy 1
#' wlr.power.maxcombo(T = NULL, events = c(60, 80), alpha=alpha, 
#'   power = NULL, side = 1, r = 1, n = 100, 
#'   h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
#'   f.ws = fws1, F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' #study design using weighting strategy 2
#' wlr.power.maxcombo(T = NULL, events = c(60, 80), alpha=alpha, 
#'   power = NULL, side = 1, r = 1, n = 100, 
#'   h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
#'   f.ws = fws2, F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' #Simulations for exploring 7 weighting strategies
#' H0 = simulate.nphDesign(nSim=5, N = 100, A = 21, w=1.5, r=1, lambda0=log(2)/11.7, lambda1=log(2)/11.7*0.745,
#' cuts=NULL, targetEvents = c(60, 80), 
#' sf = "LDOF", overall.alpha = 0.025, side = 1, alpha = NULL,
#' logrank="Y", fws=fws5)
#' 
#' @export 
simulation.nphDesign.pwexp = function(nSim=10000, N = 672, A = 21, w=1.5, r=1, lambda0=log(2)/11.7, lambda1=log(2)/11.7*0.745, 
    cuts=NULL, dropOff0=0, dropOff1=0, targetEvents = c(290, 397, 496), 
    sf = "LDOF", overall.alpha = 0.025, alpha = NULL,
    logrank="N", fws.options=NULL, H0 = "N", parallel=TRUE, n.cores=8, seed=2022) {

    set.seed(seed)
    side = 1 #always one-sided
  #Simulation for checking type I error
  if (H0 == "Y"){lambda1 = lambda0} 
  
  ##############################
  #M Options of test strategies
  M = length(fws.options)
  
  #K analyses
  K=length(fws.options[[1]])
  
  timing = targetEvents / targetEvents[K]
  
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(alpha) && !is.null(overall.alpha)){
    ld.obf = function(s){
      a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))
      return(a)
    }
    ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
    
    
    if (sf == "LDOF"){
      gs.alpha = ld.obf(s = timing)
    }
    if (sf == "LDPK") {
      gs.alpha = ld.pk(s = timing)
    }
    if (K == 1){alpha = overall.alpha} else{
      alpha[1] = gs.alpha[1]
      for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
    }
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
      dati = simulation.pwexp(nSim=1, N = N, A = A, w=w, r=r, lambda0=lambda0, lambda1=lambda1, 
                       cuts=cuts, dropOff0=dropOff0, dropOff1=dropOff1, targetEvents = targetEvents)
 
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
      dati = simulation.pwexp(nSim=1, N = N, A = A, w=w, r=r, lambda0=lambda0, lambda1=lambda1, 
                       cuts=cuts, dropOff0=dropOff0, dropOff1=dropOff1, targetEvents = targetEvents)
      
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
  o$wlr.simulations = wlr.sim
  if(logrank=="Y"){
    lr.pow = rep(NA, K)
    for (j in 1:K) {lr.pow[j] = sum(lr.sim[,j,5])/nSim}

    lr.overall.pow = 0
    for (i in 1:nSim){lr.overall.pow = lr.overall.pow + as.numeric(sum(lr.sim[i,,5])>0)}
    lr.overall.pow = lr.overall.pow / nSim
    
    o$lr.overall.power = lr.overall.pow
    o$lr.power = lr.pow
    o$lr.simulations = lr.sim
  }
  return(o)
}


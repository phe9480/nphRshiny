#' Covariance Matrix of Two Weighted Log-rank Tests At Two Analyses Times
#'
#' This function calculates the covariance matrix of two weighted log-rank tests at different analysis times.
#' The two weight functions are specified by stabilized Fleming-Harrington class
#' with parameters (rho, gamma, tau, s.tau), where tau and s.tau are thresholds for
#' survival time and survival rates, respectively. Either tau or s.tau can be specified.
#' tau = Inf or s.tau = 0 reduces to the Fleming-Harrington test (rho, gamma).
#' User-defined weight functions f.ws1 and f.ws2 can be used as well. For example,
#' f.ws1 = function(s){s^rho*(1-s)^gamma} is equivalent to Fleming-Harrington (rho, gamma)
#' test with parameters specified for rho and gamma. The first weighted log-rank statistic has
#' weight function w1 at one analysis and the second weighted log-rank test at a later analysis has weight function w2.
#' 
#' @param  time  Survival time
#' @param  event Event indicator; 1 = event, 0 = censor
#' @param  group Treatment group; 1 = experimental group, 0 = control
#' @param  rho1 Parameter for Fleming-Harrington (rho1, gamma1) weighted log-rank test.
#' @param  gamma1 Parameter for Fleming-Harrington (rho1, gamma1) weighted log-rank test.
#'         For log-rank test, set rho1 = gamma1 = 0.
#' @param  tau1  Cut point for stabilized FH test, sFH(rho1, gamma1, tau1); with weight
#'       function defined as w1(t) = s_tilda1^rho1*(1-s_tilda1)^gamma1, where
#'       s_tilda1 = max(s(t), s.tau1) or max(s(t), s(tau1)) if s.tau1 = NULL
#'       tau1 = Inf reduces to regular Fleming-Harrington test(rho1, gamma1)
#' @param  s.tau1  Survival rate cut S(tau1) at t = tau1; default 0.
#'       s.tau1 = 0 reduces to regular Fleming-Harrington test(rho1, gamma1)
#' @param  f.ws1  Self-defined weight function of survival rate. 
#'         For example, f.ws1 = function(s){1/max(s, 0.25)}
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#' @param  rho2 Parameter for Fleming-Harrington (rho2, gamma2) weighted log-rank test.
#' @param  gamma2 Parameter for Fleming-Harrington (rho2, gamma2) weighted log-rank test.
#'         For log-rank test, set rho2 = gamma2 = 0.
#' @param  tau2  Cut point for stabilized FH test, sFH(rho2, gamma2, tau2); with weight
#'       function defined as w2(t) = s_tilda2^rho2*(1-s_tilda2)^gamma2, where
#'       s_tilda2 = max(s(t), s.tau2) or max(s(t), s(tau2)) if s.tau2 = NULL
#'       tau2 = Inf reduces to regular Fleming-Harrington test(rho2, gamma2)
#' @param  s.tau2  Survival rate cut S(tau2) at t = tau2; default 0.
#'       s.tau2 = 0 reduces to regular Fleming-Harrington test(rho2, gamma2)
#' @param  f.ws2  Self-defined weight function of survival rate. 
#'         For example, f.ws2 = function(s){1/max(s, 0.25)}. 
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#' @param  strata1 Stratification variable 1
#' @param  strata2 Stratification variable 2
#' @param  strata3 Stratification variable 3
#'        
#' @return An object with dataframes below.
#' \describe{
#' \item{data}{dataframe with variables: time1, event1, group, time2, event2.}
#' \item{uni.event.time1}{dataframe at analysis time 1 with variables}
#'          \itemize{
#'          \item u.Ne:    Number of unique event times;
#'          \item u.eTime: Unique event times;
#'          \item Y0:      Risk set of control arm at each of unique event times;
#'          \item Y1:      Risk set of experimental arm at each of unique event times;
#'          \item Y:       Risk set of pooled data at each of unique event times;
#'          \item dN0:     Event set of control arm at each of unique event times;
#'          \item dN1:     Event set of experimental arm at each of unique event times;
#'          \item dN:      Event set of pooled data at each of unique event times;
#'          \item s1:      Survival time of pooled data by KM method;
#'          \item w1:      Weight function w1(t) at each of unique event times;
#'          \item w2t1:    Weight function w2(t) evaluated at each of unique event times at analysis time 1;
#'          \item V11:     Variance statistic at each of unique event times for w1;
#'          \item V12:     Covariance statistic at each of unique event times for weighted log-rank test 1 and 2;
#'          }
#' \item{uni.event.time2}{dataframe at analysis time 2 with variables}
#'          \itemize{
#'          \item u.Ne:    Number of unique event times;
#'          \item u.eTime: Unique event times;
#'          \item Y0:      Risk set of control arm at each of unique event times;
#'          \item Y1:      Risk set of experimental arm at each of unique event times;
#'          \item Y:       Risk set of pooled data at each of unique event times;
#'          \item dN0:     Event set of control arm at each of unique event times;
#'          \item dN1:     Event set of experimental arm at each of unique event times;
#'          \item dN:      Event set of pooled data at each of unique event times;
#'          \item s2:      Survival time of pooled data by KM method;
#'          \item w2:      Weight function w1(t) at each of unique event times;
#'          \item V22:     Variance statistic at each of unique event times for weighted log-rank test 2;
#'          }
#' \item{corr}{Correlation between two weigthed log-rank score statistics U1 and U2 evaluated at two analysis times, equivalent to the correlation between two normalized weigthed log-rank statistics Z1 and Z2, where Zi = Ui/sqrt(var(Ui))}
#' \item{cov}{Covariance between two weighted log-rank score statistics, U1 and U2  evaluated at two analysis times.}
#' \item{summary.events}{dataframe of events summary}
#'          \itemize{
#'          \item total.events1:    Total number of events at analysis time 1
#'          \item total.events2:    Total number of events at analysis time 2
#'          \item event.ratio:      Ratio of events for analysis time 1 vs 2
#'          \item sqrt.event.ratio: Square root of event.ratio. With large samples, 
#'          when the same weighted log-rank test is used at two analysis times, the
#'          square root of event.ratio is approximately the correlation, estimated 
#'          based on its consistent estimator, between the two weighted log-rank 
#'          test statistics evaluated at two different times.
#'          Note: sqrt.event.ratio is not an approximation to the correlation when two
#'          different weighted log-rank tests are used at two different analysis times.
#'          In addition, sqrt.event.ratio doesn't consider stratified analysis.
#'          }
#' }
#' @examples
#' 
#' #Example 1. Simulate a trial of 1:1 randomization with 600 patients. Enrollment 
#' #follows a distribution of (t/24)^1.5 for t<=24 months. Control arm has exponential dist. 
#' #with a median of 12 months, and experimental arm has delayed effect at 6 months, and after 6 months
#' #HR = 0.65. No dropoff. IA and FA will be performed at 400 and 500 events repectively.
#' 
#' data = simulation.pwexp(nSim=1, N = 600, A = 24, w=1.5, r=1, lam0=log(2)/12, 
#'       lam1=c(log(2)/12, log(2)/12*0.65), cuts=6, drop0= 0, drop1= 0, 
#'       targetEvents = c(400, 500))
#' data.IA = data[[1]]; data.FA = data[[2]]
#' data.IA$strata1 = data.FA$strata1 = sample(c(1,2), 600, replace = TRUE);
#' data.IA$strata2 = data.FA$strata2 = sample(c(1,2), 600, replace = TRUE);
#' data.IA$strata3 = data.FA$strata3 = sample(c(1,2), 600, replace = TRUE);
#' 
#' group = as.numeric(data.IA$treatment == "experimental")
#' 
#' #(1) Correlation between 2 unstratified logrank tests at IA and FA. The consistent estimator is very close to sqrt(event ratio).
#' v1 = wlr.cov2t(time1=data.IA$survTimeCut, event1=1-data.IA$cnsrCut, time2=data.FA$survTimeCut, event2=1-data.FA$cnsrCut, group=group, rho1=0, gamma1=0, tau1 = NULL, s.tau1=0,rho2=0, gamma2=0, tau2 = NULL, s.tau2=0,f.ws1=NULL, f.ws2=NULL, strata1=NULL, strata2=NULL, strata3=NULL)
#' 
#' #(2) Correlation between 2 stratified logrank tests at IA and FA
#' v2 = wlr.cov2t(time1=data.IA$survTimeCut, event1=1-data.IA$cnsrCut, time2=data.FA$survTimeCut, event2=1-data.FA$cnsrCut, group=group, rho1=0, gamma1=0, tau1 = NULL, s.tau1=0,rho2=0, gamma2=0, tau2 = NULL, s.tau2=0,f.ws1=NULL, f.ws2=NULL, strata1=data.IA$strata1, strata2=data.IA$strata2, strata3=data.IA$strata3)
#' 
#' #(3) Correlation between logrank test at IA and FH(0, 1) test at FA
#' v3 = wlr.cov2t(time1=data.IA$survTimeCut, event1=1-data.IA$cnsrCut, time2=data.FA$survTimeCut, event2=1-data.FA$cnsrCut, group=group, rho1=0, gamma1=0, tau1 = NULL, s.tau1=0,rho2=0, gamma2=1, tau2 = NULL, s.tau2=0,f.ws1=NULL, f.ws2=NULL, strata1=NULL, strata2=NULL, strata3=NULL)
#' 
#' #(4) Correlation between FH(0, 1) test at IA and FH(0, 1) test at FA
#' v4 = wlr.cov2t(time1=data.IA$survTimeCut, event1=1-data.IA$cnsrCut, time2=data.FA$survTimeCut, event2=1-data.FA$cnsrCut, group=group, rho1=0, gamma1=1, tau1 = NULL, s.tau1=0,rho2=0, gamma2=1, tau2 = NULL, s.tau2=0,f.ws1=NULL, f.ws2=NULL, strata1=NULL, strata2=NULL, strata3=NULL)
#' 
#' @keywords internal 
wlr.cov2t = function(time1=c(5,7,10,12,12,15,20,20), event1=c(1,0,0,1,1,0,1,1),
                     time2=c(5,10,13,12,14,15,20,20), event2=c(1,0,1,1,1,1,1,1),
                     group=c(0,1,0,1,0,1,0,1), strata1=NULL, strata2=NULL, strata3=NULL, 
                     rho1=0, gamma1=0, tau1 = NULL, s.tau1=0,
                     rho2=0, gamma2=1, tau2 = NULL, s.tau2=0,
                     f.ws1=NULL, f.ws2=NULL) {
  #Unstratified version
  wlr.cov2t0 = function(time1=c(5,7,10,12,12,15,20,20), event1=c(1,0,0,1,1,0,1,1),
                        time2=c(5,10,13,12,14,15,20,20), event2=c(1,0,1,1,1,1,1,1),
                        group=c(0,1,0,1,0,1,0,1), rho1=0, gamma1=0, tau1 = NULL, s.tau1=0,
                        rho2=0, gamma2=1, tau2 = NULL, s.tau2=0,
                        f.ws1=NULL, f.ws2=NULL) {
    
    f.s = function(time = time1, event = event1){
      u.eTime = unique(time[event==1]) #vector of unique event times
      u.Ne = length(u.eTime) #number of unique event times
      if (u.Ne == 0) {return(NULL); stop("No events")}
      
      u.eTime = sort(u.eTime) #sort by increasing order of unique event times
      
      Y0 = Y1 = Y = dN = dN0 = dN1 = rep(NA, u.Ne) #risk-set, death-set
      
      for (i in 1:u.Ne) {
        Y0[i] = sum(time>=u.eTime[i] & group == 0)
        Y1[i] = sum(time>=u.eTime[i] & group == 1)
        dN0[i] = sum(time == u.eTime[i] & group == 0)
        dN1[i] = sum(time == u.eTime[i] & group == 1)
      }
      
      Y=Y0+Y1
      dN = dN0 + dN1
      
      s = rep(NA, u.Ne)
      
      s[1] = 1 - dN[1]/Y[1]
      if (u.Ne > 1) {
        for (i in 2:u.Ne) {
          s[i] = s[i-1] * (1 - dN[i]/Y[i])
        }
      }
      o = list()
      o.out = data.frame(cbind(u.Ne, u.eTime, Y0, Y1, Y, dN0, dN1, dN))
      o$s = s; o$out=o.out
      return(o)
    }
    
    f.w = function(s=s, u.eTime=u.eTime, f.ws=f.ws1, tau=tau1, s.tau=s.tau1, rho=rho1, gamma=gamma1) {
      u.Ne = length(s)
      w = rep(NA, u.Ne)
      #If f.ws() function is provided, then use directly. 
      if(!is.null(f.ws)){
        w[1] = f.ws(1)
        for(i in 2:u.Ne){
          w[i] = f.ws(s[i-1])
          s.til=NULL
        }
      } else {
        #Find S(tau): = S(max(ti)|ti <= tau), where ti is unique event time.
        if (is.null(s.tau)) {
          if (is.null(tau)){stop("tau or s.tau or f.ws is required for weight function.")}else{
            s.tau = 1
            for (i in 1:u.Ne) {
              if (u.eTime[i] <= tau) {s.tau = s[i]} else {break}
            }
          }
        } #if s.tau is missing, find s.tau by tau.
        s.til = apply(cbind(s, s.tau),MARGIN=1,FUN=max)
        
        if(gamma == 0){w[1] = 1} else{w[1] = 0}; 
        for (i in 2:u.Ne){ w[i] = s.til[i-1]^rho*(1-s.til[i-1])^gamma }
        
        #w = s.til^rho*(1-s.til)^gamma
      }
      ot = list()
      ot$w = w; ot$s.til = s.til
      return(ot)
    }
    
    #w1: weight function based on (time1, event1; rho1, gamma1, tau1, s.tau1, f.ws1)
    #w2: weight function based on (time2, event2; rho2, gamma2, tau2, s.tau2, f.ws2)
    o.fs1 = f.s(time = time1, event = event1)
    o.fs2 = f.s(time = time2, event = event2)
    o.fw1 = f.w(s=o.fs1$s, u.eTime=o.fs1$out$u.eTime, f.ws=f.ws1, tau=tau1, s.tau=s.tau1, rho=rho1, gamma=gamma1)
    #w2 function evaluated at 1st analysis
    o.fw2t1 = f.w(s=o.fs1$s, u.eTime=o.fs1$out$u.eTime, f.ws=f.ws2, tau=tau2, s.tau=s.tau2, rho=rho2, gamma=gamma2)
    o.fw2 = f.w(s=o.fs2$s, u.eTime=o.fs2$out$u.eTime, f.ws=f.ws2, tau=tau2, s.tau=s.tau2, rho=rho2, gamma=gamma2)
    
    w1 = o.fw1$w; w2t1 = o.fw2t1$w; w2 = o.fw2$w; 
    
    V = matrix(NA, nrow=2, ncol=2)  
    s1.Y0 = o.fs1$out$Y0; s1.Y1 = o.fs1$out$Y1; s1.Y = o.fs1$out$Y; s1.dN = o.fs1$out$dN;
    s2.Y0 = o.fs2$out$Y0; s2.Y1 = o.fs2$out$Y1; s2.Y = o.fs2$out$Y; s2.dN = o.fs2$out$dN;
    v11 = w1*w1*s1.Y0*s1.Y1/(s1.Y^2)*(s1.Y-s1.dN)/(s1.Y-1)*s1.dN
    v12 = w1*w2t1*s1.Y0*s1.Y1/(s1.Y^2)*(s1.Y-s1.dN)/(s1.Y-1)*s1.dN
    v22 = w2*w2*s2.Y0*s2.Y1/(s2.Y^2)*(s2.Y-s2.dN)/(s2.Y-1)*s2.dN
    
    V[1,1]= sum(v11[!is.nan(v11)])
    V[1,2]=V[2,1]=sum(v12[!is.nan(v12)])
    V[2,2]=sum(v22[!is.nan(v22)])
    
    corr = V[1,2]/(sqrt(V[1,1]*V[2,2]))
    
    #create a dataframe to output the list of unique event times and statistics
    s1 = o.fs1$s; s2 = o.fs2$s
    uni.event.time1 = data.frame(cbind(o.fs1$out, s1, w1, w2t1, v11, v12))
    uni.event.time2 = data.frame(cbind(o.fs2$out, s2, w2, v22))
    
    #create a dataframe to output the original data
    data = data.frame(cbind(time1, event1, time2, event2, group))
    
    #create a dataframe to output the other statistics
    total.events1 = sum(event1==1)
    total.events2 = sum(event2==1)
    event.ratio = total.events1/total.events2
    sqrt.event.ratio = sqrt(event.ratio)
    
    summary.events = data.frame(cbind(total.events1, total.events2, event.ratio, sqrt.event.ratio))
    
    o=list()
    o$uni.event.time1 = uni.event.time1
    o$uni.event.time2 = uni.event.time2
    
    o$corr = corr
    o$cov = V
    o$summary.events = summary.events
    
    if(!is.null(f.ws1)){wt1 = f.ws1} else{
      wt1 = data.frame(cbind(rho1, gamma1, tau1, s.tau1))
    }
    if(!is.null(f.ws2)){wt2 = f.ws2} else{
      wt2 = data.frame(cbind(rho2, gamma2, tau2, s.tau2))
    }
    o$wt1 = wt1
    o$wt2 = wt2   
    return(o)
  }
  
  #Initialize strata1-3
  n = length(time1)
  if(is.null(strata1)){strata1 = rep(1,n)}
  if(is.null(strata2)){strata2 = rep(1,n)}
  if(is.null(strata3)){strata3 = rep(1,n)}
  u.strata1 = unique(strata1)
  u.strata2 = unique(strata2)
  u.strata3 = unique(strata3)
  
  n.strata1 = length(u.strata1)
  n.strata2 = length(u.strata2)
  n.strata3 = length(u.strata3)
  n.strata = n.strata1*n.strata2*n.strata3
  
  #unstratified analysis
  if (n.strata == 1) {
    o = wlr.cov2t0(time1=time1, event1=event1, time2=time2, event2=event2, group=group, 
                   rho1=rho1, gamma1=gamma1, tau1 = tau1, 
                   s.tau1=s.tau1, rho2=rho2, gamma2=gamma2, tau2 = tau2, 
                   s.tau2=s.tau2,f.ws1=f.ws1,f.ws2=f.ws2)
  }else{
    uni.event.time1=uni.event.time2=data=NULL; cov=0
    for (i in 1:n.strata1) {
      for (j in 1:n.strata2){
        for (k in 1:n.strata3){
          #Identifier of the unique stratum in combination of strata1, strata2, and strata3
          ix = (strata1==u.strata1[i] & strata2 == u.strata2[j] & strata3 == u.strata3[k])
            #Calculation for non-empty strata
            if(sum(ix) != 0){
            #survival data at 1st analysis within each stratum
            time1.s = time1[ix]
            event1.s = event1[ix]
            #survival data at 2nd analysis within each stratum
            time2.s = time2[ix]
            event2.s = event2[ix]
            group.s = group[ix]
            #Find the covariance within each stratum
            wlr.s = wlr.cov2t0(time1=time1.s, event1=event1.s, time2=time2.s, event2=event2.s, group=group.s, 
                             rho1=rho1, gamma1=gamma1, tau1 = tau1, 
                             s.tau1=s.tau1, rho2=rho2, gamma2=gamma2, tau2 = tau2, 
                             s.tau2=s.tau2,f.ws1=f.ws1,f.ws2=f.ws2)
            #Record the strata variables
            wlr.s$uni.event.time1$strata1 = u.strata1[i]
            wlr.s$uni.event.time1$strata2 = u.strata2[j]
            wlr.s$uni.event.time1$strata3 = u.strata3[k]
          
            wlr.s$uni.event.time2$strata1 = u.strata1[i]
            wlr.s$uni.event.time2$strata2 = u.strata2[j]
            wlr.s$uni.event.time2$strata3 = u.strata3[k]
          
            uni.event.time1 = rbind(uni.event.time1, wlr.s$uni.event.time1)
            uni.event.time2 = rbind(uni.event.time2, wlr.s$uni.event.time2)
          
            wlr.s$data$strata1 = u.strata1[i];
            wlr.s$data$strata2 = u.strata2[j]
            wlr.s$data$strata3 = u.strata3[k]
            data = rbind(data, wlr.s$data);
          
            cov = cov + wlr.s$cov
          }
        }  
      }       
    }  
    
    corr = cov[1,2]/(sqrt(cov[1,1]*cov[2,2]))    
    #create a dataframe to output the other statistics
    total.events1 = sum(event1==1)
    total.events2 = sum(event2==1)
    event.ratio = total.events1/total.events2
    sqrt.event.ratio = sqrt(event.ratio)
    
    summary.events = data.frame(cbind(total.events1, total.events2, event.ratio, sqrt.event.ratio))
    
    o=list()
    o$uni.event.time1 = uni.event.time1
    o$uni.event.time2 = uni.event.time2
    o$cov = cov
    o$corr = corr
    o$summary.events = summary.events
  }
  return(o)
}


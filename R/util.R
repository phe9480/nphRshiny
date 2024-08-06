#' Utility Functions
#'
#' This is a collection of utility functions.
#' 
#' 
#' @noRd
#' @keywords internal
lr = function(s){1}
fh01 = function(s){1-s}
fh11 = function(s){s*(1-s)}
sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 

fw <- function(s, rho, gamma, s.tau){
  s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max)
  w = s.til^rho*(1-s.til)^gamma
  return(w)
}

###########################
#Recruitment data utility
###########################
f.nEachMonth = function (N=600, A=24, w=2, r=2, Lambda=NULL) {
  
  N1 = N * (r/(r+1))
  N0 = N - N1
  
  #When r > 1, the control arm has smaller number of pts. 
  #Just need to determine enrollment for control arm per month, 
  #then to obtain enrollment for experimental arm by n1i = n0i * r.
  
  n1 = n0 = rep(NA, A) #enrollment by month
  randdt1 = rep(NA, N1) #randomization date
  randdt0 = rep(NA, N0)
  
  #Determine number of pts per month for control arm
  #(i-1)th month cumulative enrolled pts
  cLastN0 = 0
  for (i in 1:A) {
    #ith month: cumulative #pts
    if (is.null(Lambda)){
      cN0i = max(round((i/A)^w * N0), 1)
    } else {
      cN0i = max(round(Lambda(i/A) * N0), 1)
    }
    
    n0[i] = max(cN0i - cLastN0, 1)
    if (i == A) {n0[i] = N0 - sum(n0[1:(A-1)]) }
    cLastN0 = cN0i  
  }
  n1 = n0 * r
  
  #Patch for extreme rare scenarios that 0 enrollment in the last month
  if(n0[A] == 0 && n0[A-1] > 1){n0[A-1] = n0[A-1]-1; n0[A]=1}
  if(n1[A] == 0 && n1[A-1] > 1){n1[A-1] = n1[A-1]-1; n1[A]=1}
  
  o = list()
  o$n0 = n0
  o$n1 = n1
  return(o)
}
###########################
#Data cut utility
###########################
f.dataCut = function(data, targetEvents = 397, DCO = NULL) {
  data0 = data
  data0.order <- data0[order(as.numeric(data0$calendarTime)), ] #order by calendar time
  data.event <- data0.order[data0.order$cnsr == 0,] #Events Only
  
  data.event$event.seq <- seq.int(nrow(data.event)) #event sequence number
  if(is.null(DCO)){
    #Data cutoff in calendar time added to the original dataframe as a variable
    data0$calendarCutoff = as.numeric(data.event$calendarTime[data.event$event.seq == targetEvents])
  } else {
    data0$calendarCutoff = DCO
  }
  data0$survTimeCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), as.numeric(data0$survTime), as.numeric(data0$calendarCutoff) - as.numeric(data0$enterTime))
  data0$cnsrCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), data0$cnsr, 1)
  
  return(data0)
}


#Note:lr(s) = fw(s, rho=0, gamma=0, s.tau=0)
#   fh01(s) = fw(s, rho=0, gamma=1, s.tau=0)
#   fh11(s) = fw(s, rho=1, gamma=1, s.tau=0)
###################################################
#To use these internal functions, use nphRshiny:::lr, nphRshiny:::fh01, nphRshiny:::fh11
#Alternatively, directly load these functions when start the R Shiny app.


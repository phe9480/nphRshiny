#' Recruitment Utility Function to Determine the Number of Patients Enrolled per Month
#'
#' Determine the Number of Patients Enrolled per Month according to specified recruitment pattern.
#' 
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param Lambda Cumulative distribution function (CDF) for enrollment on (0, infinity). For example, uniform enrollment of 20 patients / month for 24 months has Lambda = function(t){t/24*as.numeric(t<= 24) + as.numeric(t>24)}.   
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' 
#' @return An object with elements:
#' \describe{
#' \item{n0}{number of subjects per month for control arm}
#' \item{n1}{number of subjects per month for experimental arm}
#' }
#' 
#' @export 
#' 
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

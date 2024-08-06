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


#Note:lr(s) = fw(s, rho=0, gamma=0, s.tau=0)
#   fh01(s) = fw(s, rho=0, gamma=1, s.tau=0)
#   fh11(s) = fw(s, rho=1, gamma=1, s.tau=0)
###################################################
#To use these internal functions, use nphRshiny:::lr, nphRshiny:::fh01, nphRshiny:::fh11
#Alternatively, directly load these functions when start the R Shiny app.


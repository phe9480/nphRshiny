#' Incremental Alpha for Using Haybittle-Peto Method in Group Sequential Design
#' 
#' The incremental alpha is the type I error allocated to a specific analysis while the sequential tests are negative for all previous analyses. Only one-sided type I error (alpha) is considered.
#' Specifically, incremental alpha_k = P(reject H0 at kth analysis while not rejecting H0 at all previous analyses)
#' Haybittle-Peto method has a fixed small p value rejection boundary at all IAs, which can preserve a large type I error for final analysis.
#' 
#' @param  overall.alpha  Allocated overall alpha (one-sided) for group sequential design
#' @param  timing Information fractions of analyses, for example timing = c(0.6, 0.75, 1.0).
#' @param  p1 A fixed p value boundary for IAs (one-sided)
#' 
#' @examples 
#' 
#' #(1) Group sequential tests at information times 0.5, 0.7, 1.0 with overall alpha 0.025.
#' #IA1 and IA2 have a rejection boundary p = 0.0003 (one-sided).
#'  
#' HP(p1 = 3e-04, overall.alpha = 0.025, timing = c(0.5, 0.7, 1))
#' 
#' (2) HP method at 0.7 information fraction after 0.0001/2 administrative alpha spent
#' HP(p1 = 0.002/2, overall.alpha = 0.025-0.0001/2, timing = c(0.7, 1))
#' 
#' HP(p1 = 0.002/2, overall.alpha = 0.025, timing = c(0.7, 1)) #without admin alpha spent
#' 
#' @return An object with dataframe below.
#'  \itemize{
#'  \item  p: Rejection boundary in p value (one-sided)
#'  \item  z: Rejection boundary in z value
#'  \item  alpha Incremental alpha(one-sided). 
#'  }
#'  
#' @export
#' 

HP <- function (p1 = 3e-04, overall.alpha = 0.025, timing = c(0.5, 0.7, 1)) {
  M = length(timing)
  corr = matrix(1, nrow = M, ncol = M)
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      corr[i, j] = corr[j, i] = sqrt(timing[i]/timing[j])
    }
  }
  if (M == 1) {
    pf = p1 = overall.alpha
    zf = z1 = qnorm(1 - overall.alpha)
  }
  else {
    z1 = qnorm(1 - p1)
    a = rep(NA, M)
    if (M == 2) {
      a[1] = p1
      a[2] = overall.alpha - a[1]
    }
    if (M > 2) {
      a[1] = p1
      for (i in 2:(M - 1)) {
        a[i] = mvtnorm::pmvnorm(lower = c(rep(-Inf, i - 1), z1), 
                                upper = c(rep(z1, i - 1), Inf), 
                                corr = corr[1:i, 1:i], 
                                abseps = 1e-08, maxpts = 1e+05)[1]
      }
      a[M] = overall.alpha - sum(a[1:(M - 1)])
    }
    f.x = function(x) {
      I = mvtnorm::pmvnorm(lower = c(rep(-Inf, M - 1), x), 
                           upper = c(rep(z1, M - 1), Inf), 
                           corr = corr, 
                           abseps = 1e-08, maxpts = 1e+05)[1]
      return(I - a[M])
    }
    zf = uniroot(f = f.x, interval = c(1, 10), tol = 1e-08)$root
    pf = 1 - pnorm(zf)
  }
  p = c(rep(p1, M - 1), pf)
  z = c(rep(z1, M - 1), zf)
  alpha = a
  
  o = data.frame(cbind(p, z, alpha))
  return(o)
}



  

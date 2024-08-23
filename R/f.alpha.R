#' Incremental Alpha for Group Sequential Design
#' 
#' The incremental alpha is the type I error allocated to a specific analysis while the sequential tests are negative for all previous analyses. Only one-sided type I error (alpha) is considered.
#' Specifically, incremental alpha_k = P(reject H0 at kth analysis while not rejecting H0 at all previous analyses)
#' 
#' @param  overall.alpha  Allocated overall alpha (one-sided) for group sequential design
#' @param  sf Spending function. Acceptable types include: (1) LanDeMets O'Brien Fleming: "LDOF", (2) LanDeMets Pocock: "LDPK", (3) Hwang-Shih-DeCani: "HSD" with parameter param. (4) Haybittle-Peto: "Haybittle-Peto". (5) Bespoke: "Bespoke"
#' @param  timing Information fractions of analyses, for example timing = c(0.6, 0.75, 1.0).
#' @param param parameter for Hwang-Shih-DeCani spending function
#' @param p1 A fixed p value boundary for IAs (one-sided), which is applicable to Haybittle-Peto alpha spending only.
#' @param cum.alpha Cumulative alpha spending by analysis, which is applicable to Bespoke method only. Cum.alpha must have the same length as timing.
#' 
#' @examples 
#' #Group sequential tests at information times 0.5, 0.7, 1.0 with overall alpha 0.025.
#' 
#' #(1) LD OBF spending
#' f.alpha(overall.alpha=0.025, sf="LDOF", timing=c(0.5, 0.7, 1.0))
#' 
#' #(2) LD Pocock spending
#' f.alpha(overall.alpha=0.025, sf="LDPK", timing=c(0.5, 0.7, 1.0))
#' 
#' #(3) Hwang-Shih-DeCani(-3)
#' f.alpha(overall.alpha=0.025, sf="HSD", timing=c(0.5, 0.7, 1.0), param=-3)
#' 
#' #(4) Haybittle-Peto with 0.0003 rejection boundaries at IAs
#' f.alpha(overall.alpha = 0.025, sf="Haybittle-Peto", timing = c(0.5, 0.7, 1), p1 = 0.0003)
#' 
#' #(5) Bespoke method with cumulative alpha spending as c(0.0003, 0.01, 0.025). In this method, the last element of cum.alpha must be equal to overall.alpha.
#' f.alpha(overall.alpha = 0.025, sf="Bespoke", timing = c(0.5, 0.7, 1), cum.alpha = c(0.0003, 0.01, 0.025))
#' 
#' #(6) Alpha spending for IA vs information fraction
#' t = seq(0.01, 1, by=0.01)
#' a = rep(NA, length(t))
#' for (i in 1:length(t)){
#'   a[i] = f.alpha(overall.alpha = 0.025, sf="LDOF", timing = c(t[i], 1), p1 = NULL)[1]
#' }
#' plot(t, a, type="l", xlab="Information Time", ylab = "IA Alpha Spending")
#' 
#' @export
#' 
f.alpha = function(overall.alpha=0.025, sf="LDOF", timing=c(0.75, 1), p1=NULL, cum.alpha=NULL, param=-3){
  K = length(timing)
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  
  ld.obf = function(s, overall.alpha=overall.alpha){
    2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))
  }
  ld.pk = function(s, overall.alpha=overall.alpha){overall.alpha * log(1 + (exp(1)-1)*s)}
  hsd <- function(s, overall.alpha=overall.alpha, param=-2){ if (param == 0) s * overall.alpha else overall.alpha * (1-exp(-s*param))/(1-exp(-param))}
  
  if (sf == "LDOF"){
    gs.alpha = ld.obf(s = timing, overall.alpha=overall.alpha)
  }
  if (sf == "LDPK") {
    gs.alpha = ld.pk(s = timing, overall.alpha=overall.alpha)
  }
  if (sf == "HSD") {
    gs.alpha = hsd(s = timing, overall.alpha=overall.alpha, param=param)
  }
  if (sf == "Bespoke"){
    gs.alpha = cum.alpha
  }
  
  if (sf == "Haybittle-Peto") {
    alpha = HP(p1 = p1, overall.alpha = overall.alpha, timing = timing)$alpha
  } else {
    if (K == 1){alpha = overall.alpha} else{
    alpha = rep(NA, K)
    alpha[1] = gs.alpha[1]
    for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
  }}
  
  return(alpha)
}


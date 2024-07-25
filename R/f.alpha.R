#' Incremental Alpha for Group Sequential Design
#' 
#' The incremental alpha is the type I error allocated to a specific analysis while the sequential tests are negative for all previous analyses. Only one-sided type I error (alpha) is considered.
#' Specifically, incremental alpha_k = P(reject H0 at kth analysis while not rejecting H0 at all previous analyses)
#' 
#' @param  overall.alpha  Allocated overall alpha (one-sided) for group sequential design
#' @param  sf Spending function. LanDeMets O'Brien Fleming: "LDOF", LanDeMets Pocock: "LDPK", "HSD": Hwang-Shih-DeCani spending function with parameter param.
#' @param  timing Information fractions of analyses, for example timing = c(0.6, 0.75, 1.0).
#' @param param parameter for Hwang-Shih-DeCani spending function
#' 
#' @examples 
#' (1) Group sequential tests at information times 0.5, 0.7, 1.0 with overall alpha 0.025.
#' 
#' #LD OBF spending
#' f.alpha(overall.alpha=0.025, sf="LDOF", timing=c(0.5, 0.7, 1.0))
#' 
#' #LD Pocock spending
#' f.alpha(overall.alpha=0.025, sf="LDPK", timing=c(0.5, 0.7, 1.0))
#' 
#' #HSD(-3)
#' f.alpha(overall.alpha=0.025, sf="HSD", timing=c(0.5, 0.7, 1.0), param=-3)
#' 
#' 
#' @export
#' 
f.alpha = function(overall.alpha=0.025, sf="LDOF", timing=c(0.75, 1), param=-3){
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
  
  if (K == 1){alpha = overall.alpha} else{
    alpha = rep(NA, K)
    alpha[1] = gs.alpha[1]
    for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
  }
  
  return(alpha)
}


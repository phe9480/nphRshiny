#' Display of Survival Curves per Study Design
#' 
#' This function plots the survival curves per study design
#' 
#' @param t  A sequence of survival time for the plot (x-axis)
#' @param S  An object containly a list of survival functions to draw
#' @param ... Other graphic parameters passed to the plot
#'
#' @return Display of the graph
#'  
#' @examples 
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' #Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' 
#' plot_S(Tmax = 50, S = list(S0, S1))
#' 
#' @export
#' 
#' 
plot_S = function(S = list(S0, S1), Tmax = 50,
                  leg=list(x=30, y=1, txt=c("Control Arm", "Experimental Arm")),
                  param = list(xlab="Time Since First Subject Randomized (mo)", 
                             ylab="Survival",
                             main="Survival Curve Per Study Design")){
  
  #number of survival curves
  g = length(S)
  t = seq(0, Tmax, by = 0.1)

  s0 = apply(cbind(t), MARGIN=1,FUN=S[[1]])
  s1 = apply(cbind(t), MARGIN=1,FUN=S[[2]])
  
  plot(t, seq(min(s0, s1), max(s0, s1), length.out=length(t)), type = "n", 
       bty = "l", xlab = param$xlab, 
       ylab=param$ylab, main = param$main)
  
  col.seq = c("seagreen3","blue3","turquoise4","deeppink3","orange")
  
  for (i in 1:g){
    si = apply(cbind(t), MARGIN=1,FUN=S[[i]])
    if (g <= 5){
      lines(t, si, lwd = 4, col = col.seq[i])
    }else{
      lines(t, si, lwd = 4, col = i)
    }
  }
  abline(h = seq(0, 1, 0.1), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
 
  if (g<= 5){
    legend(leg$x, leg$y, leg$txt, col=col.seq[1:g], lwd=4, bty="n", cex=0.8)
  } else {
    legend(leg$x, leg$y, leg$txt, col=1:g, lwd=4, bty="n", cex=0.8)
  }
}

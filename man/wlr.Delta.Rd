% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wlr.Delta.R
\name{wlr.Delta}
\alias{wlr.Delta}
\title{The Expectation of U(t)/n in Weighted Logrank Test Z = U(t) / sqrt(V(t))}
\usage{
wlr.Delta(
  DCO = 24,
  r = 1,
  h0 = function(t) {
     log(2)/12
 },
  S0 = function(t) {
     exp(-log(2)/12 * t)
 },
  h1 = function(t) {
     log(2)/12 * 0.7
 },
  S1 = function(t) {
     exp(-log(2)/12 * 0.7 * t)
 },
  cuts = NULL,
  rho = 0,
  gamma = 0,
  tau = NULL,
  s.tau = 0,
  f.ws = NULL,
  Lambda = function(t) {
     (t/18) * as.numeric(t <= 18) + as.numeric(t > 18)
 },
  G0 = function(t) {
     0
 },
  G1 = function(t) {
     0
 }
)
}
\arguments{
\item{DCO}{Analysis time, calculated from first subject in.}

\item{r}{Randomization ratio of experimental arm : control arm as r:1.
When r = 1, it is equal allocation. Default r = 1.}

\item{h0}{Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.}

\item{S0}{Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
but providing S0(t) can improves computational efficiency and
usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).}

\item{h1}{Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.}

\item{S1}{Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
but providing S1(t) can improves computational efficiency and
usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).}

\item{cuts}{A vector of cut points to define piecewise distributions.
If cuts is not specified or incorrectly specified, it might occasionally have numerical integration issue.}

\item{rho}{Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.}

\item{gamma}{Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
For log-rank test, set rho = gamma = 0.}

\item{tau}{Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)}

\item{s.tau}{Survival rate cut S(tau) at t = tau1; default 0.5, ie. cut at median.
s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)}

\item{f.ws}{Self-defined weight function of survival rate, eg, f.ws = function(s){1/max(s, 0.25)}
When f.ws is specified, sFH parameter will be ignored.}

\item{Lambda}{Cumulative distribution function of enrollment.}

\item{G0}{Cumulative distribution function of drop-off for control arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
is the distribution function for 3 percent drop-off in 12 months of followup.}

\item{G1}{Cumulative distribution function of drop-off for experimental arm, eg, G.ltfu=function(t){1-exp(-0.03/12*t)}
is the distribution function for 3 percent drop-off in 12 months of followup.}
}
\description{
This function calculates the asymptotic mean of U(t)/n,
where the weighted logrank test Z = U(t) / sqrt(V(t)).
Delta = E\link{U(t)/n}
Note: V(t)/n --> sigma2, Z = n^(-1/2)*U(t)/sqrt(V(t)/n).
E(Z) = mu = n^(1/2)*Delta/sigma
}
\examples{
#1:1 randomization, control~exp(median=12);
#experiment~delayed effect 6 mo and HR 0.6 after delay
#Fleming-Harrington (1,1) test; Enrollment: 18 mo with weight 1.5
#drop-off rate: 3\% in 1-year of follow-up in control arm; no drop-off in exp arm
HR = 0.6; delay = 6; lam0 = log(2) / 12; 
h0 = function(t){lam0}; S0 = function(t){exp(-lam0 * t)}
h1.D6 = function(t){lam0*as.numeric(t < delay)+HR*lam0*as.numeric(t >= delay)}
c = exp(-delay*lam0*(1-HR)); 
S1.D6 = function(t){exp(-lam0*t)*as.numeric(t<delay) + c*exp(-HR*lam0*t)*as.numeric(t>=delay)}
Lambda = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
drop0 = 0.03/12; drop1 = 0

wlr.Delta(DCO = 24, r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0})

wlr.sigma2(DCO = 24, r = 1, h0 = h0, S0 = S0,h1 = h1.D6,S1 = S1.D6, 
cuts=c(6), rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
Lambda = Lambda, G0 = function(t){1-exp(-drop0 * t)}, G1 = function(t){0})

}
\keyword{internal}

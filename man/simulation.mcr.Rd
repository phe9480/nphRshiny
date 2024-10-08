% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.mcr.R
\name{simulation.mcr}
\alias{simulation.mcr}
\title{Simulation of Two-Arm Trial Data Following Mixture Cure Rate Distribution With Delayed Effect and Non-uniform Accrual Pattern and Lost to Follow-up}
\usage{
simulation.mcr(
  nSim = 100,
  N = 600,
  A = 18,
  w = 1.5,
  Lambda = NULL,
  r = 1,
  p = c(0.1, 0.1),
  alpha = c(log(2)/12, log(2)/12),
  beta = c(1, 1),
  gamma = c(1, 1),
  lambda = c(0, 0),
  tau = c(0, 0),
  psi = c(1, 1),
  drop = c(0, 0),
  targetEvents = c(300, 420),
  DCO = NULL
)
}
\arguments{
\item{nSim}{Number of trials}

\item{N}{Total number patients in two arms.}

\item{A}{Total accrual period in months.}

\item{w}{Weight parameter in cumulative enrollment pattern. The
cumulative enrollment at month t is (t / A)^w, eg,
at month 6, the enrollment is N*(6/24)^2 = N/16 for
24 months planned accrual period.}

\item{Lambda}{Cumulative distribution function (CDF) for enrollment on (0, infinity). For example, uniform enrollment of 20 patients / month for 24 months has Lambda = function(t){t/24*as.numeric(t<= 24) + as.numeric(t>24)}.}

\item{r}{Randomization ratio \code{r:1}, where r refers to the experimental arm,
eg, \code{r=2} in 2:1 ratio}

\item{p}{Cure rate parameter. When p = 0, it reduces to \code{GMW(alpha, beta, gamma, lambda)}
distribution. Each of the following parameters is a vector of 2 components for
two treatment groups (control, experimental arm):
\code{p}, \code{alpha}, \code{beta}, \code{gamma}, \code{lambda}, \code{tau}, \code{psi}, \code{drop}.}

\item{alpha}{Generalized modified Weibull (GMW) distribution parameters. \code{alpha > 0}}

\item{beta}{Generalized modified Weibull (GMW) distribution parameters. \code{beta > 0}}

\item{gamma}{Generalized modified Weibull (GMW) distribution parameters. \code{gamma >= 0} but \code{gamma} and \code{lambda} cannot be both 0.}

\item{lambda}{Generalized modified Weibull (GMW) distribution parameters. \code{lambda >= 0} but \code{gamma} and `lambda`` cannot be both 0.}

\item{tau}{Threshold for delayed effect period. \code{tau = 0} reduces no delayed effect.}

\item{psi}{Hazard ratio after delayed effect. \code{psi = 1} reduces to the survival function without proportional hazards after delayed period (\code{tau}).}

\item{drop}{Drop-off rate per unit time. For example, if 3\% drop off in 1 year of followup, then \code{drop = 0.03/12}.}
}
\value{
A dataframe for each analysis including the following variables:
\describe{
\item{sim}{sequence number of simulated dataset;}
\item{treatment}{treatment group with values of "control" and "experimental"}
\item{enterTime}{Time of randomization in calendar time}
\item{calendarTime}{the time when event/censoring occurred in calendar time}
\item{survTime}{Survival time for analysis, = calendarTime - enterTime}
\item{cnsr}{censor status (0=event; 1=censor) before administrative censoring due to data cut}
\item{calendarCutOff}{Data CutOff Time (DCO);}
\item{survTimeCut}{Survival time after cut}
\item{cnsrCut}{Censor status after cut}
}
}
\description{
Consider the survival function of mixture cure rate (MCR) distribution:
\verb{S(t) = p + (1-p)S0(t)}, where \code{S0(t)} is a survival distribution
for susceptible subject, i.e., S0(0)=1 and S0(t) -> 0 as t -> Infinity.
S0(t) can be any proper survival function. For generality of S0(t), consider
the generalized modified Weibull (GMW) distribution with parameters
(alpha, beta, gamma, lambda) Martinez et al(2013).
}
\details{
alpha: scale parameter
beta and gamma: shape parameters
lambda: acceleration parameter

S0(t) = 1-(1-exp(-alpha\emph{t^gamma}exp(lambda*t)))^beta

Special cases:
(1) Weibull dist: lambda = 0, beta = 1. Beware of the parameterization difference.
(2) Exponential dist: lambda = 0, beta = 1, gamma = 1. The shape parameter (hazard rate) is alpha.
(3) Rayleigh dist: lambda = 0, beta = 1, gamma = 2.
(4) Exponentiated Weibull dist (EW): lambda = 0
(5) Exponentiated exponential dist (EE): lambda = 0 and gamma = 1
(6) Generalized Rayleigh dist (GR): lambda = 0, gamma = 2
(7) Modified Weibull dist (MW): beta = 1

Let T be the survival time according to survival function S1(t) below.
Denote T's distribution as MCR(p, alpha, beta, gamma, lambda, tau, psi),
where tau is the delayed effect and psi is the hazard ratio after delayed effect,
ie. proportional hazards to S(t) after delay tau.
\code{S(t) = p + (1-p)*S0(t)}
\verb{S1(t) = S(t)I(t<tau) + S(tau)^(1-psi)*S(t)^psi}
In brief, \code{S0(t)} is the proper GMW distribution \verb{(alpha, beta, gamma, lambda)};
\code{S(t)} is MCR with additional cure rate parameter \code{p};
In reference to \code{S(t)}, \code{S1(t)} is a delayed effect distribution and proportional hazards after delay.
}
\examples{
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
set.seed(2022)

#(6) Delayed effect 6 months, proportional hazards HR = 0.6
set.seed(2022)
data = simulation.mcr(nSim=10, N = 600, A = 18, w=1.5, r=1, p=c(0.1,0.1), 
alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
targetEvents = c(300, 420), DCO = NULL) 

#same as above  
data = simulation.mcr(nSim=10, N = 600, Lambda= function(t){(t/18)^1.5*as.numeric(t<= 18) + as.numeric(t>18)}, r=1, p=c(0.1,0.1), 
alpha = c(log(2)/12,log(2)/12), beta=c(1,1), gamma=c(1,1), 
lambda=c(0,0), tau=c(0,6), psi=c(1,0.6), drop=c(0,0),
targetEvents = c(300, 420), DCO = NULL) 

data.IA = data[[1]][data[[1]]$sim==1,]
data.FA = data[[2]][data[[2]]$sim==1,]
m0 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=0, psi=1)
m1 = qmcr(0.5, p=0.1, alpha = log(2)/12, beta=1, gamma=1, lambda=0, tau=6, psi=0.6)

km.IA<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data.IA)
plot(km.IA, ylab="Survival")
abline(h = c(0.1, 0.5), col="gray80", lty=2)
abline(v=c(6, m1, m0), col="gray80", lty=2)
 
km.FA<-survival::survfit(survival::Surv(survTimeCut,1-as.numeric(cnsrCut))~treatment, data=data.FA)
plot(km.FA, ylab="Survival") 
abline(h = c(0.1, 0.5), col="gray80", lty=2)
abline(v=c(6, m1, m0), col="gray80", lty=2)

}

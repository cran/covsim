% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vita.R
\name{vita}
\alias{vita}
\title{Calibrate a regular vine}
\usage{
vita(
  margins,
  sigma.target,
  vc = NULL,
  family_set = c("clayton", "gauss", "joe", "gumbel", "frank"),
  Nmax = 10^6,
  numrootpoints = 10,
  conflevel = 0.995,
  numpoints = 4,
  verbose = TRUE,
  cores = parallel::detectCores()
)
}
\arguments{
\item{margins}{A list where each element corresponds to a margin. Each
margin element is a list containing the distribution family ("distr") and
additional parameters. Must be a distribution available in the stats package.}

\item{sigma.target}{The target covariance matrix that is to be matched. The diagonal
elements must contain the variances of marginal distributions.}

\item{vc}{A vine dist object as specified by the rvinecopulib package. This object
specifies the  vine that is to be calibrated. If not provided, a D-vine is assumed.}

\item{family_set}{A vector of one-parameter pair-copula families that is to
be calibrated at each node in the vine. Possible entries are "gauss", "clayton", "joe", "gumbel" and "frank".
Calibration of pair-copula families is attempted in the order provided.}

\item{Nmax}{The sample size used for calibration. Reduce for faster calibration,
at the cost of precision.}

\item{numrootpoints}{The number of estimated roots at the initial calibration stage, which
determines a search interval where Nmax samples are drawn}

\item{conflevel}{Confidence level for determining search interval}

\item{numpoints}{The number of samples drawn with size Nmax, to determine the root within search interval
To increase precision increase this number. To calibrate faster (but less precisely), may be reduced to a number no lower than 2}

\item{verbose}{If TRUE, outputs details of calibration of each bicopula}

\item{cores}{Number of cores to use. If larger than 1, computations are done in parallel. May be determined with parallel:detectCores()}
}
\value{
If a feasible solution was found, a vine to be used for simulation
}
\description{
\code{vita} implements the VITA (VIne-To-Anything) algorithm.
Covariance matrix and margins are specified, and \code{vita} calibrates the
pair-copulas in each node of the tree to match the target covariance.
}
\examples{
set.seed(1)# define a target covariance. 3 dimensions.
sigma.target <- cov(MASS::mvrnorm(10, mu=rep(0,3), Sigma=diag(1, 3)))

#normal margins that match the covariances:
marginsnorm <- lapply(X=sqrt(diag(sigma.target)),function(X) list(distr="norm", sd=X) )

#calibrate with a default D-vine, with rather low precision (default Nmax is 10^6)
# if cores=1 is removed, all cores are used, with a speed gain
calibrated.vine <- vita(marginsnorm, sigma.target =sigma.target, Nmax=10^5, cores=1)
#check
#round(cov(rvinecopulib::rvine(10^5, calibrated.vine))-sigma.target, 3)

#margins are normal but dependence structure is not
#pairs(rvinecopulib::rvine(500, calibrated.vine))



}
\references{
Grønneberg, S., Foldnes, N., & Marcoulides, K. M. (2021). covsim: An r package for simulating non-normal data for structural equation models using copulas. Journal of Statistical Software. doi:10.18637/jss.v102.i03
}

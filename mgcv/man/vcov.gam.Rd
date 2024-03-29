\name{vcov.gam}
\alias{vcov.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract parameter (estimator) covariance matrix from GAM fit}
\description{ Extracts the Bayesian posterior covariance matrix of the
parameters or frequentist covariance matrix of the parameter estimators 
from a fitted \code{gam} object.
}
\usage{
\method{vcov}{gam}(object, sandwich=FALSE, freq = FALSE, dispersion = NULL,unconditional=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{fitted model object of class \code{gam} as produced by \code{gam()}.}
\item{sandwich}{compute sandwich estimate of covariance matrix. Currently expensive for discrete bam fits.}
\item{freq}{ \code{TRUE} to return the frequentist covariance matrix of the
parameter estimators, \code{FALSE} to return the Bayesian posterior covariance
matrix of the parameters. The latter option includes the expected squared bias
according to the Bayesian smoothing prior.} 
\item{dispersion}{ a value for the dispersion parameter: not normally used.}
\item{unconditional}{ if \code{TRUE} (and \code{freq==FALSE}) then the Bayesian 
smoothing parameter 
uncertainty corrected covariance matrix is returned, if available. }
\item{...}{ other arguments, currently ignored.}
}
\details{ Basically, just extracts \code{object$Ve}, \code{object$Vp} or \code{object$Vc} (if available) from a \code{\link{gamObject}}, unless \code{sandwich==TRUE} in which case the sandwich estimate is computed (with or without the squared bias component). 
}

\value{ A matrix corresponding to the estimated frequentist covariance matrix
of the model parameter estimators/coefficients, or the estimated posterior
covariance matrix of the parameters, depending on the argument \code{freq}.
}

\author{ Henric Nilsson. 
Maintained by Simon N. Wood \email{simon.wood@r-project.org}
}

\references{
Wood, S.N. (2017) Generalized Additive Models: An Introductio with R (2nd ed) CRC Press
}


\seealso{ \code{\link{gam}}}


\examples{ 
require(mgcv)
n <- 100
x <- runif(n)
y <- sin(x*2*pi) + rnorm(n)*.2
mod <- gam(y~s(x,bs="cc",k=10),knots=list(x=seq(0,1,length=10)))
diag(vcov(mod))
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...







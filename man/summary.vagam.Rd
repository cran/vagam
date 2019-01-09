\name{summary.vagam}
\alias{summary.vagam}
\alias{print.vagam}
\title{
Summary of generalized additive model (GAM) fitted using variational approximations (VA).
}
\description{
A summary of the results from applying \code{vagam}.
}
\usage{
\method{summary}{vagam}(object, ...)

\method{print}{vagam}(x,...)
}
\arguments{
  \item{object}{An object of class "vagam".}
  \item{x}{An object of class "vagam".}
  \item{...}{Not used.}
  }
\value{
A list (some of which is printed) containing the following elements:
\itemize{
\item{call:}{The matched call.}
\item{para.coeff:}{The estimated regression coefficients corresponding to the covariates in parametric component of the GAM. This includes the intercept term.}
\item{smooth.coeff:}{The estimated smoothing coefficients corresponding to the (P-spline bases set up for) covariates in the nonparametric component of the GAM. This corresponds to the mean vector of the variational distribution.}
\item{smooth.param:}{The estimated smoothing coefficients, or the fixed smoothing parameters if \code{lambda} was supplied.}
\item{phi:}{The estimated residual variance when the Gaussian distribution is assumed for the response.}
\item{logLik:}{The maximized value of the variational log-likelihood.}
\item{family:}{The assumed distribution for the response.}
\item{smooth.stat:}{A small table of summary statistics for the nonparametric component of the GAM, including an approximate Wald-type hypothesis test for the significance of each nonparametric covariate.}
\item{para.stat:}{If \code{para.se=TRUE}, then a small table containing summary statistics for the estimated parametric component of the GAM, including an approximate Wald-type hypothesis test for the significance of each parameteric covariate.}
} }
\author{
\packageAuthor{vagam}
}
\seealso{
\code{\link{vagam} for the main fitting function}
}
\examples{
## Please see examples in the help file for the vagam function.
}

\name{gamsim}
\alias{gamsim}
\title{
Simulate example datasets from a generalized additive models (GAM).
}
\description{
This function is a modification from example 7 of the gamSim function available in the mgcv package (Wood, 2017), which is turn is Gu and Wahba 4 univariate example with correlated predictors. Please see the source code for exactly what is simulated. The function is primarily used as the basis for conducting the simulation studies in Hui et al., (2018). 
}
\usage{
gamsim(n = 400, extra.X = NULL, beta = NULL, dist = "normal", scale = 1, offset = NULL)
}
\arguments{
  \item{n}{Sample size.}
  \item{extra.X}{Extra covariates, including critically an intercept if is to be included in the linear predictor for the GAM.}
  \item{beta}{Regression coefficient estimates.}
  \item{dist}{Currently only the "normal", "poisson" or "binomial" corresponding to the binomial distributions are available.}
  \item{scale}{Scale parameter in the Normal distribution.}
  \item{offset}{This can be used to specify an a-priori known component to be included in the linear predictor during fitting.  This should be \code{NULL} or a numeric vector of length equal to \code{n}.}
}
\value{
A data frame containing information such as the simulated responses, covariates, each of the 4 "truth" smooths, and the overall linear predictor. 
}
\references{
\itemize{
\item Hui, F. K. C., You, C., Shang, H. L., and Mueller, S. (2018). Semiparametric regression using variational approximations, \emph{Journal of the American Statistical Association}, \bold{forthcoming}.

\item Wood, S. N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman and Hall/CRC.
} }
\author{
\packageAuthor{vagam}
}
\seealso{
\code{\link{vagam} for the main fitting function}
}
\keyword{datagen}
\examples{
normal_dat = gamsim(n = 40, dist = "normal", 
                    extra.X = data.frame(int = rep(1,40), trt = rep(c(0,1), each = 20)),
                    beta = c(-1, 0.5))

pois_dat = gamsim(n = 40, dist = "poisson", 
                  extra.X = data.frame(int = rep(1, 40), trt = rep(c(0,1), each = 20)),
                  beta = c(-1, 0.5))

binom_dat = gamsim(n = 40, dist = "binomial",
                   extra.X = data.frame(int = rep(1, 40), trt = rep(c(0,1), each = 20)),
                   beta = c(0, 0.5))

## Please see examples in the help file for the vagam function.
}

\name{predict.vagam}
\alias{predict.vagam}
\title{
Predictions from a fitted generalized additive model (GAM). 
}
\description{
Takes a fitted vagam object produced by the main vagam function and produces predictions given a new set of values for the model covariates or the original values used for the model fit.
}
\usage{
\method{predict}{vagam}(object, new.smoothX, new.paraX = NULL, terms = NULL, 
alpha = 0.05, type = "link", ...)
}
\arguments{
  \item{object}{An object of class "vagam".}
  \item{new.smoothX}{A new matrix of covariates, each of which are were entered as additive smooth terms in the fitted GAM.}
  \item{new.paraX}{A new matrix of covariates, each of which were be entered as parametric terms in the fitted GAM. Note the predictions will account for the intercept ONLY if new.paraX is supplied.}
  \item{terms}{If \code{terms = NULL}, prediction is made across all smoothing covariates. Else, prediction is made to the specified smoothing covariate, with no intercept added.}
  \item{alpha}{Level of significance for the pointwise confidence bands for predictions.}
  \item{type}{When \code{type = "link"} (default) the linear predictor (with associated standard errors) is returned. When \code{type = "response"} predictions on the scale of the response are returned (with associated standard errors).}
  \item{...}{This is currently ignored.}
}
\details{
Current implemented a basic method of constructing predictions either for a single smoothing covariate, or across all the smoothing (and parametric if supplied) covariates, based on a GAM fitted using the main vagam function. By default, standard errors and this pointwise confidence bands are also produced based on fitted vagam object. Under the variational approximations framework, the smooths and confidence bands are constructed based on the variational approximation to the posterior distribution of the smoothing coefficients (which are treated as random effects with a normal prior under the mixed model framework). Please see Hui et al., (2018) for more information. 
}
\value{
A data frame containing information such as the predicted response, standard errors, and lower and upper bounds of the pointwise confidence bands.
}
\references{
\itemize{
\item Hui, F. K. C., You, C., Shang, H. L., and Mueller, S. (2018). Semiparametric regression using variational approximations, \emph{Journal of the American Statistical Association}, \bold{forthcoming}.
} }
\author{
\packageAuthor{vagam}
}
\seealso{
\code{\link{vagam} for the main fitting function}
}
\keyword{methods}
\examples{
## Please see examples in the help file for the vagam function.
}

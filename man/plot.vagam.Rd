\name{plot.vagam}
\alias{plot.vagam}
\title{
Basic plots for a fitted generalized additive model (GAMs). 
}
\description{
Takes a fitted vagam object produced by the main vagam function and plots the component smooth functions that make it up, on the scale of the linear predictor.
}
\usage{
\method{plot}{vagam}(x, n = 100, alpha = 0.05, rug = TRUE, se = TRUE, 
xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL, 
select = NULL, ...)
}
\arguments{
  \item{x}{An object of class "vagam".}
  \item{n}{Number of observations used in constructing predictions for plotting.}
  \item{alpha}{Level of significance for the pointwise confidence bands for predictions.}
  \item{rug}{If \code{rug = TRUE}, adds a rug representation (1-d plot) of the data to the plot.}
  \item{se}{If \code{se = TRUE}, adds lower and upper bounds of 95\% pointwise confidence bands.}
  \item{xlim}{Range of plotting for x variable.}
  \item{ylim}{Range of plotting for y variable.}
  \item{xlab}{Label for x variable.}
  \item{ylab}{Label for y variable.}
  \item{main}{Title of plotting.}
  \item{select}{Select which observations to plot.}
  \item{\dots}{Other plotting arguments.}
}
\details{
Currently implements a basic plot for each of the fitted smoothers from a GAM fitted using the main vagam function. This is done by making use of the predict function to construct the fitted smooths. There is also the option of adding pointwise confidence bands based on fitted vagam object. Under the variational approximations framework, the smooths and confidence bands are constructed based on the variational approximation to the posterior distribution of the smoothing coefficients (which are treated as random effects with a normal prior under the mixed model framework). Please see Hui et al., (2018) for more information. 
}
\value{
The functions main purpose is its side effect of generating a set of plots.
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
\keyword{hplot}
\examples{
## Please see examples in the help file for the vagam function.
}

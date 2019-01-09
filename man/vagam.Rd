\name{vagam}
\alias{vagam}
\title{
Fitting generalized additive models (GAMs) using variational approximations (VA). 
}
\usage{
vagam(y, smooth.X, para.X = NULL, lambda = NULL, int.knots, family = gaussian(), 
A.struct = c("unstructured", "block"), offset = NULL, save.data = FALSE, 
para.se = FALSE, doIC = FALSE, 
control = list(eps = 0.001, maxit = 1000, trace = TRUE, seed.number = 123, 
mc.samps = 4000, pois.step.size = 0.01))
}
\description{
Follows the variational approximation approach of Hui et al. (2018) for fitted generalized additive models. In this package, the term GAM is taken to be generalized linear mixed model, specifically, the nonparametric component is modeled using a P-splines i.e., cubic B-splines with a first order difference penalty. Because the penalty can be written as a quadratic form in terms of the smoothing coefficients, then it is treated a (degenerate) multivariate normal random effects distribution and a marginal log-likleihood for the resulting mixed model can be constructed. 

The VA framework is then utilized to provide a fully or at least closed to fully tractable lower bound approximation to the marginal likelihood of a GAM. In doing so, the VA framework aims offers both the stability and natural inference tools available in the mixed model approach to GAMs, while achieving computation times comparable to that of using the penalized likelihood approach to GAMs. 
}
\arguments{
  \item{y}{A response vector.}
  \item{smooth.X}{A matrix of covariates, each of which are to be entered as additive smooth terms in the GAM.}
  \item{para.X}{An optional matrix of covariates, each of which are to be entered as parametric terms in the GAM. Please note that NO intercept term needs to be included as it is included by default.}
  \item{lambda}{An optional vector of length \code{ncol(smooth.X)}, where each element corresponds to the smoothing parameter to be applied to the respective covariate in \code{smooth.X}. If supplied, then the GAM is fitted with the smoothing parameters held fixed at this values. If \code{lambda=NULL}, then smoothing parameters for all covariates to be smoothed are updated automatically as part of the VA algorithm.}
  \item{int.knots}{Either a single number of a vector of length \code{ncol(smooth.X)}, corresponding to the number of interior knots to be use for the respective covariate iin \code{smooth.X}. This argument is passed to the function \code{smooth.construct} from the mgcv package (Wood, 2017) in order to construct the P-splines bases. Equally spaced knots based on quantiles are used.}
  \item{family}{Currently only the \code{gaussian(link = "identity"), poisson(link = "log")}, and \code{binomial(link = "logit")} corresponding to Bernoulli distributions are available.}
  \item{A.struct}{The assumed structure of the covariance matrix in the variational distribution of the smoothing coefficients. Currently, the two options are \code{A.struc = "unstructured"} corresponding to assuming an fully unstructured covariance matrix, and \code{A.struc = "block"} which assumes a block diagonal structure where the all covariances between different smoothing covariates are assumed to be zero (but the covariance submatrix remains unstructured within the spline basis functions for a selected smoothing covariate). The latter is sub-optimal in the sense that the most appropriate variational distribution should use a completely unstructured covariance matrix , but MAY (but is not guaranteed) save computation time especially when the number of smoothing covariates and/or the number of interior knots is very large.}
  \item{offset}{This can be used to specify an a-priori known component to be included in the linear predictor during fitting.  This should be \code{NULL} or a numeric vector of length equal to \code{length(y)}.}  
  \item{save.data}{If \code{save.date=TRUE}, then the returned vagam object will also include \code{y, smooth.X, para.X}, and the full matrix of P-spline basis functions.}
  \item{para.se}{If \code{para.se=TRUE}, the standard errors based on the VA approach are returned for any covariates in \code{para.X} that are included as parametric terms. Note that if \code{para.se=FALSE} then a standard error for the intercept term will not be returned even though an intercept term is included by default.}
  \item{doIC}{If \code{doIC=TRUE}, then the AIC and BIC are returned, where the AIC is calculated as \eqn{AIC = -2\times} variation log-likelihood \eqn{+ 2\times trace(H)} with \eqn{trace(H)} bring a measure of the degrees of freedom of the model as based on the hat-matrix arising from iterative reweighted least squares, and the BIC replaces the 2 with \code{log(length(y))} for the model complexity penalty; please see Wood (2017) for more details. Note however that this out is largely mute as the VA approach provides an automatic method of selecting the smoothing parameters, meaning an external approach such as information criteria is not required.}
  \item{control}{A list controlling the finer details of the VA approach for fitting GAMs. These include:
  \itemize{
  \item{mc.samps:}{This controls Monte Carlo samples for calculating variational observed information matrix using Louis' method}
  \item{seed:}{This controls seed for starting values of the fitting algorithm in general}
  \item{pois.step.size:}{This controls step size for penalized iterative reweighted least squares (P-IRLS) portion of the VA approach when \code{family=poisson()}. This may be tweaked to use smaller step sizes as the approach here can be a tad unstable especially if there is possible overdispersion.}
  }
}
}
\details{
Please note that the package is still in its early days, and only a very basic form of GAMs with purely additive
terms and P-splines is fitted. The function borrows heavily from the excellent software available in the mgcv
package (Wood, 2017), in the sense that it uses the \code{smooth.construct} function with \code{bs = "ps"} to
set up the matrix of P-splines bases (so cubic B-splines with a first order difference penalty matrix) along with
imposing the necessary centering constraints. With these ingredients, it then maximizes the variational 
log-likelihood by iteratively updating the model and variational parameters. The variational log-likelihood
is obtained by proposing a variational distribution for the smoothing coefficients (in this case, a multivariate
normal distribution between unknown mean vector and covariance matrix), and then minimizing the Kullback-Leibler
distance between this variational distribution and the true posterior distribution of the smoothing coefficients.
In turn, this is designed to be (closed to) fully tractable lower bound approximation to the true marginal
log-likelihood for the GAM, which for non-normal responses does not possess a tractable form. Note that in
contrast to the marginal log-likelihood or many approximations such the Laplace approximation and adaptive
quadrature, the variational approximation typically presents a tractable form that is relatively straightforward
to maximize. At the same time, because it takes views the GAM as a mixed model, then it also possesses nice
inference tools such as an approximate posterior distribution of the smoothing coefficients available immediately
from maximizing the VA log-likelihood, and automatic choice of the smoothing parameters. We refer to readers to
Wood (2017) and Ruppert et al. (2003) for detailed introductions to GAMs and how many of them can be set up as
mixed models; Eilers and Marx (1996) for the seminal text on P-splines, and Hui et al. (2018) for the text on 
which this package is based.
}
\value{
An object of vagam class containing one or more of the following elements:
\itemize{
\item{call:}{The matched call.}
\item{kappa:}{The estimated regression coefficients corresponding to the covariates in \code{para.X}. This includes the intercept term.}
\item{a:}{The estimated smoothing coefficients corresponding to the (P-spline bases set up for) covariates in \code{smooth.X}. This corresponds to the mean vector of the variational distribution.}
\item{A:}{The estimated posterior covariance of the smoothing coefficients corresponding to the (P-spline bases set up for) covariates in \code{smooth.X}. This corresponds to the covariance matrix of the variational distribution.}
\item{lambda:}{The estimated smoothing parameters, or the fixed smoothing parameters if \code{lambda} was supplied.}
\item{IC:}{A vector containing the calculated values of and AIC and BIC if \code{doIC=TRUE}. Note this is largely a mute output.}
\item{phi:}{The estimated residual variance when \code{family=gaussian()}.}
\item{linear.predictors:}{The estimated linear predictor i.e., the parametric plus nonparametric component.}
\item{logL:}{The maximized value of the variational log-likelihood.}
\item{no.knots:}{The number of interior knots used, as per \code{int.knots}.}
\item{index.cov:}{A vector indexing which covariate each column in the final full matrix P-spline bases belongs to.}
\item{basis.info:}{A list with length equal to \code{ncol(smooth.X)}, with each element being the output from an application of \code{smooth.construct} to construct the P-spline for a selected covariate in \code{smooth.X}.}
\item{y, para.X, smooth.X, Z:}{Returned in \code{save.data=TRUE}. Note critically that \code{Z} final full matrix P-spline basis functions.} 
\item{smooth.stat:}{A small table of summary statistics for the nonparametric component of the GAM, including an approximate Wald-type hypothesis test for the significance of each nonparametric covariate.}
\item{para.stat:}{If \code{para.se=TRUE}, then a small table containing summary statistics for the estimated parametric component of the GAM, including an approximate Wald-type hypothesis test for the significance of each parameteric covariate.}
\item{obs.info:}{If \code{para.se=TRUE}, then the estimated variational observed information matrix for the parameteric component of the GAM; please see Hui et al. (2018) for more information.}
}
}
\references{
\itemize{
\item Eilers, P. H. C., and Marx, B. D. (1996) Flexible Smoothing with B-splines and Penalties.  Statistical Science, 11: 89-121

\item Hui, F. K. C., You, C., Shang, H. L., and Mueller, S. (2018). Semiparametric regression using variational approximations, \emph{Journal of the American Statistical Association}, \bold{forthcoming}.

\item Ruppert, D., Wand M. P., and Carroll, R. J. (2003) Semiparametric Regression. Cambridge University Press.

\item Wood, S. N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman and Hall/CRC.
} }
\author{
\packageAuthor{vagam}
}
\seealso{
\code{\link{summary.vagam}} for a basic summary of the fitted model; \code{\link{plot.vagam}} for basic plotting the component smooths; \code{\link{predict.vagam} for basic prediction}
}
\keyword{methods}
\examples{
## Example 1: Application to wage data
data(wage_data)

south_code <- gender_code <- race_code <- union_code <- vector("numeric", nrow(wage_data))
union_code[wage_data$union == "member"] <- 1
south_code[wage_data$south == "yes"] <- 1
gender_code[wage_data$gender == "female"] <- 1
race_code[wage_data$race == "White"] <- 1
para.X <- data.frame(south = south_code, gender = gender_code, race = race_code)

fit_va <- vagam(y = union_code, smooth.X = wage_data[,c("education", "wage", "age")],
                       para.X = para.X,
                       int.knots = 8, save.data = TRUE, 
                       family = binomial(), 
                       para.se = TRUE)
summary(fit_va)

a <- 1
par(mfrow = c(1, 3), las = 1, cex = a, cex.lab = a+0.2, cex.main = a+0.5, mar = c(5,5,3,2))
plot(fit_va, ylim = c(-2.7, 2.7), select = 1, 
                        xlab = "Education", ylab = "Smooth of Education", lwd = 3)
plot(fit_va, ylim = c(-2.7, 2.7), select = 2, 
                        xlab = "Wage", ylab = "Smooth of Wage", main = "Plots from VA-GAM", lwd = 3)
plot(fit_va, ylim = c(-2.7, 2.7), select = 3, 
                        xlab = "Age", ylab = "Smooth of Age", lwd = 3)

                      
\dontrun{
## Example 2: Simulated data with size = 50 and compare how GAMs can be fitted 
## in VA and mgcv (which uses penalized quasi-likelihood)
choose_k <- 5 * ceiling(50^0.2)
true_beta <- c(-1, 0.5)

poisson_dat <- gamsim(n = 50, dist = "poisson", extra.X = data.frame(intercept = rep(1,50), 
                        treatment = rep(c(0,1), each = 50/2)), beta = true_beta)

## GAM using VA
fit_va <- vagam(y = poisson_dat$y, smooth.X = poisson_dat[,2:5], 
		         para.X = data.frame(treatment = poisson_dat$treatment), 
                          int.knots = choose_k, save.data = TRUE, family = poisson(), 
                          para.se = TRUE)
summary(fit_va)
                       
## GAM using mgcv with default options
fit_mgcv1 <- gam(y ~ treatment + s(x0) + s(x1) + s(x2) + s(x3), 
                             data = poisson_dat, family = poisson())


## GAM using mgcv with P-splines and preset knots; 
## this is equivalent to VA in terms of the splines bases functions 
fit_mgcv2 <- gam(y ~ treatment + s(x0, bs = "ps", k = round(choose_k/2) + 2, m  = c(2,1)) + 
                             s(x1, bs = "ps", k = round(choose_k/2) + 2, m  = c(2,1)) + 
                             s(x2, bs = "ps", k = round(choose_k/2) + 2, m  = c(2,1)) + 
                             s(x3, bs = "ps", k = round(choose_k/2) + 2, m  = c(2,1)), 
                             data = poisson_dat, family = poisson())                     
}
}

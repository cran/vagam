print.vagam <- function(x, ...)
{
	cat("Variational approximation for GAMs\n")
	cat(paste("\nCall:", deparse(x$call, 200), "\n"))
	cat(paste("\nEstimated regression coefficients for parametric component:"), x$kappa)
	cat(paste("\nEstimated smoothing coefficients for nonparametric component:"), x$a)
	cat(paste("\nEstimated smoothing parameters (or fixed if lambda was supplied):"), x$lambda)
	cat(paste("\nNumber of interior knots used:"), x$no.knots)
	cat(paste("\nMaximized value of the variational log-likelihood:"), x$logL)
	cat("\n")
}

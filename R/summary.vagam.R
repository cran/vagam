summary.vagam <- function(object, ...)
{
	print(object)
	cat(paste("\nSummary statistics for nonparametric component:\n"))
	print(object$smooth.stat)
	cat(paste("\nSummary statistics for parametric component (if para.se = TRUE):\n"))
	print(object$para.stat)
	#cat(paste("\nSummary statistics for linear predictors:\n"))
	#print(summary(object$linear.predictors))
	cat("\n")
	
	out <- list(call = object$call, para.coeff = object$kappa, smooth.coeff = object$a, smooth.param = object$lambda, phi = object$phi, logLik = object$logL, family = object$family, smooth.stat = object$smooth.stat, para.stat = object$para.stat)
    class(out) <- "summary.vagam"
    return(out)
}

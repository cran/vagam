predict.vagam <- function(object, new.smoothX, new.paraX = NULL, terms = NULL, alpha = 0.05, type = "link", ...)
{
	if(class(object) != "vagam")
	{
		stop("object should be an object of class vagam. Thanks!")
	}
    x <- object
	if(alpha > 1 || alpha < 0)
	{
		stop("alpha should be a number between 0 and 1, for (1-alpha)% pointwise confidence intervals. Thanks!")
	}
	type <- match.arg(type, c("link", "response"))	
	if(!is.null(terms)) 
	{
		if(length(terms) != 1) 
			stop("terms should either be NULL, in which case the prediction is made on all terms, or should be a single number indicating which smoothing term is to be used for prediction. Thanks!")
		if(!is.vector(new.smoothX))
			stop("If terms is a single number indicating which smoothing term is to be used for prediction, then new.smoothX should be a vector corresponding to the covariate values for which predictions are to be calculated for that term. Thanks!")
	}
	if(is.null(terms))
	{
        if(ncol(new.smoothX) != length(x$no.knots))
        stop("If terms is NULL, in which case the prediction is made on all terms, then new.smoothX should be a matrix with the name number of columns as the number of smoothing covariates fitted in x. Thanks!")
	}
	if(!is.null(new.paraX))
	{
		if(any(apply(new.paraX, 2, function(x) all(x == 1))))
			stop("No intercept terms should be included in new.paraX, as this is included by default. Thanks!")
			new.paraX <- cbind(1, new.paraX) ## Intercept included by default
		if(is.null(x$para.stat))
			stop("If new.paraX is supplied, in which case the prediction involves parametric terms, then x should include the relevant standard errors i.e., x$para.stat. Thanks!")
	}

	
	if(!is.null(terms)) ## Predictions made for a single smoothing covariates
	{
		# calculate total variance for a smoothing covariate
		newZ <- Predict.matrix(x$basis.info[[terms]], data = data.frame(x0 = (new.smoothX)))
		newZ <- newZ %*% x$basis.info[[terms]]$transform_mat
		new_eta <- newZ %*% x$a[x$index.cov == terms]
		sub.A <- x$A[x$index.cov == terms, x$index.cov == terms]
		variance <- rowSums((newZ %*% sub.A) * newZ)
		
		# construct pointwise prediction interval based on normality assumption
		lbound <- new_eta - qnorm(alpha/2, lower.tail = FALSE) * sqrt(variance)
		ubound <- new_eta + qnorm(alpha/2, lower.tail = FALSE) * sqrt(variance)
	}
	
	
	if(is.null(terms)) ## Predictions made across the smoothing (and parametric if supplied) covariates. Note the predictions will account for the intercept ONLY if new.paraX is supplied
	{
		# calculate total variance for smoothing covariates
		newZ <- NULL
		for(k2 in 1:length(x$no.knots)) 
		{
			tmpZ <- Predict.matrix(x$basis.info[[k2]], data = data.frame(x0 = (new.smoothX[,k2])))
			newZ <- cbind(newZ, tmpZ %*% x$basis.info[[k2]]$transform.mat)
		}
		new_eta <- newZ %*% x$a
		variance <- rowSums((newZ %*% x$A) * newZ)

		if(!is.null(new.paraX))
		{
            new_eta <- new_eta + as.matrix(new.paraX) %*% x$kappa
            variance <- variance + rowSums((new.paraX %*% solve(x$obs.info)[1:length(x$kappa),1:length(x$kappa)]) * new.paraX)
		}
        if(!is.null(x$offset))
            new_eta <- new_eta + x$offset
		
		# construct pointwise prediction interval based on normality assumption
		lbound <- new_eta - qnorm(alpha/2, lower.tail = FALSE) * sqrt(variance)
		ubound <- new_eta + qnorm(alpha/2, lower.tail = FALSE) * sqrt(variance)
	}
	if(type == "response")
	{
		new_eta <- x$family$linkinv(new_eta)
		lbound <- x$family$linkinv(lbound)
		ubound <- x$family$linkinv(ubound)
	}
	out <- data.frame(prediction = new_eta, se = sqrt(variance), lower.bound = lbound, upper.bound = ubound,
					new.smoothX = new.smoothX)
	if(!is.null(new.paraX))
		out$new.paraX <- new.paraX

	return(out)
}

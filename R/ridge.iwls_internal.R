ridge.iwls <- function(x, y, family, penalty = NULL, initial.beta = NULL, max.steps = 1000, conv.eps = 1e-3, offset = NULL, gamma = 1, extra = FALSE)
{
	y <- as.vector(y)
	x <- as.matrix(x)
	p <- ncol(x)
	nobs <- nrow(x)
	if(is.null(penalty)) 
		penalty <- matrix(0, p, p)
	if(is.null(offset)) 
		offset <- numeric(length(y))
	if(!is.null(penalty))
	{
		if(nrow(penalty) != p || ncol(penalty) != p) 
			stop("penalty could a square matrix with dimension ncol(x) by ncol(x). Thanks!")
	}
	converged <- FALSE
	stop_at <- max.steps

	beta_mat <- matrix(0, nrow = max.steps, ncol = p)
	if(is.null(initial.beta)) 
		initial.beta <- rep(0.01, p)
	eta.new <- x %*% initial.beta + offset

	for(i in 1:max.steps)
	{
		beta_mat[i,] <- initial.beta
		mu.new <- family$linkinv(eta.new)
		d.new <- family$mu.eta(eta.new)
		v.new <- family$variance(mu.new)
		weights <- c(d.new/sqrt(v.new))
		x.star <- weights * x
		y.tilde.star <- weights * (eta.new - offset + (y - mu.new)/d.new)
		p.imat.new <- crossprod(x.star) + penalty
		# print(p.imat.new)
		inv.pimat.new <- chol2inv(chol(p.imat.new))		
		beta_new <- gamma * inv.pimat.new %*% crossprod(x.star, y.tilde.star) + (1 - gamma) * beta_mat[i,]
		
		if((sum(abs(beta_new - initial.beta))/sum(abs(initial.beta)) <= conv.eps)) 
		{
			converged <- TRUE
			stop_at <- i
			if(i < max.steps) 
				break;
		}
		else
		{
			initial.beta <- beta_new
			eta.new <- x %*% beta_new + offset
		}
		# print(c(beta_new))
	}
	
	out <- list(coefficients = as.vector(beta_new), family = family, converged = converged,
                stop_at = stop_at, linear.predictors = eta.new)

     if(extra) {
          Infmat <- inv.pimat.new %*% crossprod(x.star)
          tr.Infmat <- sum(diag(Infmat))
          out$tr.Inf <- tr.Infmat
          }
     return(out)
}

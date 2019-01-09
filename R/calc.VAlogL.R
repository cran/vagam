calc.VAlogL <-
function(y, Z, para.X, family, lambda, kappa, a, A, phi, S, d, offset = NULL, index.cov, eta = NULL) 
{
	if(is.null(offset)) 
		offset <- numeric(length(y))

    if(is.null(eta)) 
          eta <- para.X %*% kappa + Z %*% a + offset
	if(family$family == "gaussian")
	{
		out <- sum(dnorm(y, mean = eta, sd = sqrt(phi), log = TRUE)) - (0.5/phi) * sum((Z %*% A) * Z) + 0.5 * determinant(A)$mod + 0.5 * sum(d * log(lambda[lambda > 0]))
	}
	if(family$family == "binomialp")
	{
        ZAZ <- (Z %*% A) * Z
		out <- sum(dbinom(y, 1, prob = binomial(link = "probit")$linkinv(eta + 0.5 * rowSums(ZAZ)), log = TRUE)) - (0.5/phi) * sum(ZAZ) + 0.5 * determinant(A)$mod + 0.5 * sum(d * log(lambda[lambda > 0]))
	}
	if(family$family == "binomial")
	{
        ZAZ <- rowSums((Z %*% A) * Z)
		out <- sum(dbinom(y, 1, prob = binomial()$linkinv(eta + 0.5 * ZAZ), log = TRUE)) - 0.5 * sum(y * ZAZ) + 0.5 * determinant(A)$mod + 0.5 * sum(d * log(lambda[lambda > 0]))
	}
	if(family$family == "poisson") 
	{
        ZAZ <- rowSums((Z %*% A) * Z)
		out <- sum(dpois(y, lambda = exp(eta + 0.5 * ZAZ), log = TRUE)) - 0.5 * sum(y * ZAZ) + 0.5 * determinant(A)$mod + 0.5 * sum(d * log(lambda[lambda > 0]))
	}

	for(k2 in 1:length(lambda)) 
	{
		out <- out - 0.5 * lambda[k2] * sum(diag(crossprod(S[[k2]], A[index.cov == k2, index.cov == k2]))) - 0.5 * lambda[k2] * t(a[index.cov == k2]) %*% S[[k2]] %*% a[index.cov == k2]
	}
	return(as.numeric(out))
}

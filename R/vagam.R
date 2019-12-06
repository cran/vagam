######################
## TODO
## 2) Redo calculation of information matrix to do sandwich formula
## 3) Allow more spline types, specifically more straight additive splines options and then possibly tensor products?
## 4) Port ridge.iwla to use C++
## 5) Allow subset, and na.action options

vagam <-
function(y, smooth.X, para.X = NULL, lambda = NULL, int.knots, family = gaussian(), A.struct = c("unstructured", "block"), offset = NULL, save.data = FALSE, para.se = FALSE, doIC = FALSE, control = list(eps = 1e-3, maxit = 1000, trace = TRUE, seed.number = 123, mc.samps = 4000, pois.step.size = 0.01))
{
	if(!(family$family %in% c("gaussian", "poisson", "binomial"))) #"binomialp"
		stop("Specified family currently no supported. Sorry!")
 	if(family$family == "gaussian")
 		family <- gaussian(link = "identity")
 	if(family$family == "poisson")
 		family <- poisson(link = "log")
 	if(family$family == "binomial")
 		family <- binomial(link = "logit")
 	if(family$family == "binomialp")
 		family <- binomial(link = "probit")


	A.struct <- match.arg(A.struct, c("unstructured", "block"))
	ncov_smooth_X <- ncol(smooth.X)


	if(is.null(offset))
		offset <- numeric(length(y))


    if(!(length(int.knots) %in% c(1,ncov_smooth_X)))
		stop("int.knots should either be a single number of a vector with length equal to ncol(smooth.X). Thanks!")
	if(length(int.knots) == 1)
		actual_int_knots <- rep(int.knots, ncov_smooth_X)
	if(length(int.knots) > 1)
		actual_int_knots <- int.knots


	if(is.null(colnames(smooth.X)))
        colnames(smooth.X) <- paste("smoothX", 1:ncov_smooth_X, sep = "")
	if(!is.null(para.X))
	{
		para.X <- as.matrix(para.X)
		if(is.null(colnames(para.X)))
            colnames(para.X) <- paste("paraX", 1:ncol(para.X), sep = "")
		if(any(apply(para.X, 2, function(x) all(x == 1))))
			stop("No intercept terms should be included in para.X, as this is included by default. Thanks!")
		para.X <- cbind(1, para.X) ## Intercept included by default
		colnames(para.X)[1] <- "Intercept"
	}
	if(is.null(para.X))
	{
		para.X <- matrix(1, nrow = length(y), ncol = 1)
		colnames(para.X) <- "Intercept"
	}
	if(is.null(lambda))
	{
		cw_lambda <- new_lambda <- rep(2, ncov_smooth_X)
		message("Lambda updated as part of VA estimation. Yeah baby!")
	}
	if(!is.null(lambda))
	{
		if(length(lambda) != ncov_smooth_X)
			stop("lambda should a vector with length equal to ncol(smooth.X). Thanks!")
		cw_lambda <- new_lambda <- lambda
		message("Lambda given. Thanks!")
	}


	## Construct B-splines basis, penalty matrix, and some other quantities for smooth.X, by extracting attributes from a mgcv gam fit
    ## no. of total knots = no. of interior knots[1] + 2*degree + 1, but then -1 to remove the intercept term
    ## no. of basis dimension = d_j = no. of interior knots[1]+ degree + 1, where degree of spline is typically set to 3 for cubic B-splines. But then we typically -1 to remove the intercept term. Then to ensure centering constraint is the satisfied, reduce by one to produce to d_j = no. of interior knots[1] + degree - 1
	S <- basis_info <- vector("list", ncov_smooth_X)
	d <- numeric(ncov_smooth_X) ## rank of Sj
	index_cov <- NULL ## Indexes which covariate each column in Z belongs to
	Z <- NULL

	for(k in 1:ncov_smooth_X)
	{
        x0 <- smooth.X[,k]
 		get_basis <- smooth.construct(s(x0, bs = "ps", k = actual_int_knots[k]+3, m  = c(2,1)), data = data.frame(x0), knots = NULL)
 		Z_k <- get_basis$X
 		transform_mat <- qr.Q(qr(as.matrix(colSums(Z_k))), complete = TRUE)[,-1] ## Centering constraint
		Z_k <- Z_k %*% transform_mat
		index_cov <- c(index_cov, rep(k, ncol(Z_k)))
		colnames(Z_k) <- paste("smooth.X", k, 1:sum(index_cov == k), sep = "")
		Z <- cbind(Z, Z_k)
 		S[[k]] <- crossprod(transform_mat, get_basis$S[[1]]) %*% transform_mat
 		d[k] <- get_basis$rank
 		get_basis$X <- NULL;
 		get_basis$transform_mat <- transform_mat;
 		basis_info[[k]] <- get_basis
 		rm(transform_mat, get_basis, x0)
	}

	cw_kappa <- new_kappa <- numeric(ncol(para.X)) ## coefficients for parametric part
	names(new_kappa) <- names(para.X)
	cw_a <- new_a <- numeric(ncol(Z)) ## coefficients for smooth part, also VA means
	if(family$family == "poisson")
	{
		old <- .Random.seed
		on.exit({ .Random.seed = old})
		set.seed(control$seed.number)
		cw_a <- new_a <- rnorm(ncol(Z)) ## Use random coefficients as it tends to work better for Poisson
	}
	names(new_a) <- colnames(Z)
	cw_A <- new_A <- matrix(0, nrow = ncol(Z), ncol = ncol(Z)) ## VA covariance
	rownames(new_A) <- colnames(new_A) <- colnames(Z)
	cw_phi <- new_phi <- 1 ## variance
	cw_logL <- -Inf
	new_logL <- 10
	diff_logL <- 10
	diff_lambda <- 0
	counter <- 1

	if(family$family %in% c("gaussian","binomialp"))
        ZZt <- crossprod(Z)
	tic <- proc.time()


	## VA...Let it rip!
	while(diff_logL > control$eps)
	{
		if(counter > control$maxit)
               break;

		## Update kappa and a and lambda if required
		## Calculation for Poisson somewhat unstable
		ZA <- Z %*% cw_A
		new_offset <- Z %*% cw_a + as.numeric(family$family %in% c("poisson","binomial")) * 0.5 * rowSums(ZA * Z) + offset
 		fit0 <- suppressWarnings(glm.fit(x = para.X, y = y, family = family, offset = new_offset, intercept = FALSE))
 		new_kappa <- fit0$coefficients

		Q <- matrix(0, ncol(Z), ncol(Z))
		for(k2 in 1:ncov_smooth_X)
            Q[index_cov == k2,index_cov == k2] <- S[[k2]]*cw_lambda[k2]
		new_offset <- para.X %*% new_kappa + as.numeric(family$family %in% c("poisson","binomial")) * 0.5 * rowSums(ZA * Z) + offset
		if(family$family == "poisson")
        {
 			fit0 <- ridge.iwls(x = Z, y = y, family = family, penalty = Q, initial.beta = cw_a, offset = new_offset, gamma = 0.1 + min(counter * control$pois.step.size, 0.9)) ## Gradually increase step size with counter
		}
		if(family$family != "poisson")
		{
			fit0 <- ridge.iwls(x = Z, y = y, family = family, penalty = cw_phi*Q, initial.beta = cw_a, offset = new_offset)
		}
 		new_a <- fit0$coefficients
		rm(new_offset)


		## Update phi for Gaussian
		new.eta <- para.X %*% new_kappa + Z %*% new_a + offset
		if(family$family == "gaussian")
		{
			new_phi <- (sum((y - new.eta)^2) + sum(ZA * Z))/length(y)
		}


		## Update A
		if(A.struct == "unstructured")
		{
			if(family$family %in% c("gaussian","binomialp"))
				new_A <- chol2inv(chol(Q + (1/new_phi)*ZZt))
			if(family$family == "poisson")
			{
				err2 <- 10;
				cw_A <- new_A
				while(err2 > 0.01)
				{
					denom <- Z * matrix(sqrt(exp(new.eta + 0.5 * rowSums(ZA * Z))), nrow = length(y), ncol = ncol(Z), byrow = FALSE)
					denom <- Q + crossprod(denom)
					new_A <- chol2inv(chol(denom))
					err2 <- sum((cw_A - new_A)^2)/2
					cw_A <- new_A
				}
			}
			if(family$family == "binomial")
			{
				err2 <- 10;
				cw_A <- new_A
				while(err2 > 0.01)
				{
					denom <- Z * matrix(sqrt(binomial()$linkinv(new.eta + 0.5 * rowSums(ZA * Z))), nrow = length(y), ncol = ncol(Z), byrow = FALSE)
					denom <- Q + crossprod(denom)
					new_A <- chol2inv(chol(denom))
					err2 <- sum((cw_A - new_A)^2)/2
					cw_A <- new_A
				}
			}
		}
		if(A.struct == "block") ## This might actually be slower than doing the full thing!
		{
			for(k2 in 1:ncov_smooth_X)
			{
				if(family$family %in% c("gaussian", "binomialp"))
				{
					new_A[index_cov == k2, index_cov == k2] <- chol2inv(chol(Q[index_cov == k2,index_cov == k2] + (1/new_phi) * crossprod(Z[,index_cov == k2])))
 				}
				if(family$family == "poisson")
				{
					err2 <- 10;
					while(err2 > 0.01)
						{
						denom <- Z[,index_cov == k2] * matrix(sqrt(exp(new.eta + 0.5 * rowSums(ZA * Z))), nrow = length(y), ncol = sum(index_cov==k2), byrow = FALSE)
						denom <- Q[index_cov == k2, index_cov == k2] + crossprod(denom)
						new_A[index_cov == k2, index_cov == k2] <- chol2inv(chol(denom))
						err2 <- sum((cw_A - new_A)^2)/2
						cw_A <- new_A
						}
				}
				if(family$family == "binomial")
				{
					err2 <- 10;
					while(err2 > 0.01)
						{
						denom <- Z[,index_cov == k2] * matrix(sqrt(binomial()$linkinv(new.eta + 0.5 * rowSums(ZA * Z))), nrow = length(y), ncol = sum(index_cov==k2), byrow = FALSE)
						denom <- Q[index_cov == k2, index_cov == k2] + crossprod(denom)
						new_A[index_cov == k2, index_cov == k2] <- chol2inv(chol(denom))
						err2 <- sum((cw_A - new_A)^2)/2
						cw_A <- new_A
						}
				}
			}
		}

		if(is.null(lambda))
		{
			for(k2 in 1:ncov_smooth_X)
			{
				new_lambda[k2] <- d[k2]/(sum(diag(crossprod(S[[k2]], new_A[index_cov == k2, index_cov == k2]))) + tcrossprod(new_a[index_cov == k2], S[[k2]]) %*% new_a[index_cov == k2])
			}
		}


		## Calculate new VA logL
		new_logL <- calc.VAlogL(y = y, Z = Z, para.X = para.X, family = family, lambda = new_lambda, kappa = new_kappa, a = new_a, A = new_A, phi = new_phi, S = S, d = d, index.cov = index_cov, eta = new.eta)
		diff_logL <- new_logL - cw_logL
		diff_lambda <- sum((new_lambda - cw_lambda)^2)
		if(control$trace)
			cat("Iteration:", counter, "\t Current VA logL:", cw_logL, " | New VA logL:", new_logL, " | Difference:", new_logL - cw_logL, "\n")
		counter <- counter + 1
		cw_logL <- new_logL
		cw_lambda <- new_lambda
		cw_kappa <- new_kappa
		cw_a <- new_a
		cw_A <- new_A
		cw_phi <- new_phi
		}
	toc <- proc.time()


	## Calculate stuff related to external smoothing parameter selection
    ic_out <- rep(NA,2)
    names(ic_out) <- c("AIC", "BIC")
	if(doIC)
	{
        ic_out <- -2 * new_logL - new_phi + c(2 * fit0$tr.Inf * new_phi, log(length(y)) * fit0$tr.Inf * new_phi)
    }

	out <- list(kappa = new_kappa, a = new_a, A = new_A, lambda = new_lambda, IC = ic_out, phi = new_phi, linear.predictors = new.eta, offset = offset, logL = new_logL, no.knots = actual_int_knots, index.cov = index_cov, basis.info = basis_info, family = family, time.taken = toc-tic)
	if(save.data)
	{
		out$y
		out$para.X <- para.X
		out$smooth.X <- smooth.X
		out$Z <- Z
	}


	## Wald test for each smooth curve
	smooth.X.waldstat <- smooth.X.waldpval <- numeric(ncov_smooth_X)
 	new_A.inv <- chol2inv(chol(new_A))
	for(k2 in 1:ncov_smooth_X)
	{
		smooth.X.waldstat[k2] <- crossprod(new_a[index_cov == k2], new_A.inv[index_cov == k2,index_cov == k2]) %*% new_a[index_cov == k2]
		smooth.X.waldpval[k2] <- pchisq(smooth.X.waldstat[k2], df = sum(index_cov==k2), lower.tail = FALSE)
	}
	smooth.X.wald <- round(rbind(smooth.X.waldstat,smooth.X.waldpval),5)
	rownames(smooth.X.wald) <- c("Wald Statistic", "p-value")
	colnames(smooth.X.wald) <- colnames(smooth.X)
	rm(smooth.X.waldstat,smooth.X.waldpval)
	out$smooth.stat <- smooth.X.wald


	## Variatonal Observed information matrix for model params and lambas
	if(para.se)
	{
		message("Calculating information matrix for model parameters...")
 		obs_info <- info.valouis(y = y, para.X = para.X, Z = Z, kappa = new_kappa, phi = new_phi, lambda = new_lambda, a = new_a, A = new_A, S = S, d = d, index.cov = index_cov, mc.samps = control$mc.samps, family = family)
		## If there the information matrix is not +ve definite, then do a one-step to produce an unstructured A and recalculate the info matrix; hopefully this helps!
		if(any(diag(obs_info) < 0) & family$family == "poisson")
		{
		message("Redoing...")
		err2 <- 10;
		cw_A <- new_A
		#new.eta <- para.X %*% new_kappa + Z %*% new_a
		while(err2 > 0.01)
			{
			denom <- Z * matrix(sqrt(exp(new.eta + 0.5 * rowSums((Z %*% cw_A) * Z))), nrow = length(y), ncol = ncol(Z), byrow = FALSE)
			denom <- Q + crossprod(denom)
			new_A <- chol2inv(chol(denom))
			err2 <- sum((cw_A - new_A)^2)/2
			cw_A <- new_A
			}
		obs_info <- info.valouis(y = y, para.X = para.X, Z = Z, kappa = new_kappa, phi = new_phi, lambda = new_lambda, a = new_a, A = cw_A, S = S, d = d, index.cov = index_cov, mc.samps = control$mc.samps, family = family)
		}

		obs_info <- info.valouis(y = y, para.X = para.X, Z = Z, kappa = new_kappa, phi = new_phi, lambda = new_lambda, a = new_a, A = new_A, S = S, d = d, index.cov = index_cov, mc.samps = control$mc.samps, family = family)
		para_sderr <- sqrt(diag(solve(obs_info)))[1:ncol(para.X)]
		para_X_stat <- round(rbind(new_kappa, para_sderr, new_kappa/para_sderr, 2 * pnorm(abs(new_kappa/para_sderr), lower.tail = FALSE)), 5)
		rownames(para_X_stat) <- c("Estimate", "Std. Error", "Wald Statistic", "p-value")
		colnames(para_X_stat) <- colnames(para.X)
		out$para.stat <- para_X_stat
		out$obs.info <- obs_info
    }
    if(!para.se)
    {
    out$para.stat <- NULL
    }

	class(out) <- "vagam"
	out$call <- match.call()
	return(out)
}

info.valouis <- function(y, para.X, Z, kappa, phi = 1, lambda, a, A, S, d, index.cov, mc.samps = 2000, family = gaussian(), offset = NULL, seed = NULL)
{
    if(!is.null(seed))
    {
        old <- .Random.seed
        on.exit( { .Random.seed <<- old } )
        set.seed(seed)
    }
    
    
	if(is.null(offset)) 
		offset <- numeric(length(y))

		
    diminfo <- ncol(para.X) + length(lambda) + as.numeric(family$family == "gaussian")
    n <- length(y)
    mc_betas <- rmvnorm(mc.samps, mean = a, sigma = A)
    mc_betasTSmc_betas <- matrix(0, nrow = mc.samps, ncol = length(lambda))
    for(k in 1:length(lambda)) ## beta^T_j%*%S_j%*%beta_j, across the MC samples
    {
        mc_betasTSmc_betas[,k] <- rowSums((mc_betas[,index.cov==k] %*% S[[k]]) * mc_betas[,index.cov == k])
    }
    score.lambda <- matrix(0.5 * d/lambda, nrow = mc.samps, ncol = length(lambda), byrow = TRUE) - 0.5 * mc_betasTSmc_betas
    rm(mc_betasTSmc_betas)
    
    
    scorescoreT_out <- matrix(0, nrow = diminfo, ncol = diminfo) ## \sum\limits_{i=1}^n E{\partial \ell_{com}(y_i, beta) \partial \ell_{com}(y_i, beta)^T}
    neghess_out <- matrix(0, nrow = diminfo, ncol = diminfo) ## \sum\limits_{i=1}^n E{\partial^2 \ell_{com}(y_i, beta)}
    eta.para <- para.X %*% kappa
    
    if(family$family == "gaussian")
    {
        for(i in 1:mc.samps)
        {
            res <- y - eta.para - as.vector(Z %*% mc_betas[i,]) - offset
            cw_score <- c(crossprod(para.X,res)/phi, -0.5 * n/phi + 0.5/phi^2 * sum(res^2))
            cw_score <- c(cw_score, score.lambda[i,])
            scorescoreT_out <- scorescoreT_out + tcrossprod(cw_score)
        }
        scorescoreT_out <- scorescoreT_out/mc.samps
        res <- y - eta.para - Z %*% a - offset
        neghess_out[1:ncol(para.X),1:ncol(para.X)] <- crossprod(para.X)/phi
        neghess_out[1:ncol(para.X),ncol(para.X)+1] <- crossprod(para.X, res)/phi^2
        neghess_out[ncol(para.X)+1,1:ncol(para.X)] <- t(neghess_out[1:ncol(para.X),ncol(para.X)+1])
        neghess_out[ncol(para.X)+1,ncol(para.X)+1] <- -0.5 * n/phi^2 + 1/phi^3 * (sum(res^2) + sum((Z %*% A) * Z))
        diag(neghess_out)[(ncol(para.X)+2):diminfo] <- 0.5 * d/lambda^2
    }
    if(family$family == "poisson")
    {
        for(i in 1:mc.samps)
        {
            res <- y - exp(eta.para + as.vector(Z %*% mc_betas[i,]) + offset)
            cw_score <- c(crossprod(para.X, res), score.lambda[i,])
            scorescoreT_out <- scorescoreT_out + tcrossprod(cw_score)
        }
        scorescoreT_out <- scorescoreT_out/mc.samps
        res <- exp(eta.para + Z %*% a + 0.5 * rowSums((Z %*% A) * Z) + offset)
        neghess_out[1:ncol(para.X),1:ncol(para.X)] <- t(para.X) %*% (para.X * as.vector(res))
        diag(neghess_out)[(ncol(para.X)+1):diminfo] <- 0.5 * d/lambda^2
    }
    if(family$family == "binomialp")
    {
        all.eta <- eta.para + Z %*% a + offset
        seqa <- -Inf * (y == 0)
        seqa[is.nan(seqa)] <- 0
        seqb <- Inf * (y == 1)
        seqb[is.nan(seqb)] <- 0
        
        for(i in 1:mc.samps)
        {
            mc.u <- rtruncnorm(n, a = seqa, b = seqb, mean = all.eta)
            res <- (mc.u - eta.para -  as.vector(Z %*% mc_betas[i,]) - offset)
            cw_score <- c(crossprod(para.X, res), score.lambda[i,])
            scorescoreT_out <- scorescoreT_out + tcrossprod(cw_score)
        }
        scorescoreT_out <- scorescoreT_out/mc.samps
        neghess_out[1:ncol(para.X),1:ncol(para.X)] <- crossprod(para.X)
        diag(neghess_out)[(ncol(para.X)+1):diminfo] <- 0.5 * d/lambda^2
    }
    if(family$family == "binomial")
    {
        meanw <- numeric(n)
        for(i in 1:mc.samps)
        {
            res <- binomial()$linkinv(eta.para + as.vector(Z %*% mc_betas[i,]) + offset)
            cw_score <- c(crossprod(para.X, y-res), score.lambda[i,])
            scorescoreT_out <- scorescoreT_out + tcrossprod(cw_score)
            meanw <- meanw + binomial()$variance(res)
        }
        scorescoreT_out <- scorescoreT_out/mc.samps
        meanw <- meanw/mc.samps
        neghess_out[1:ncol(para.X),1:ncol(para.X)] <- crossprod(para.X, para.X * as.vector(meanw))
        diag(neghess_out)[(ncol(para.X)+1):diminfo] <- 0.5 * d/lambda^2
    }
    
    rm(cw_score, score.lambda, res)
    get.info <- neghess_out - scorescoreT_out
    return(get.info)
}

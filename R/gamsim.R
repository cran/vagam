gamsim <- function(n = 400, extra.X = NULL, beta = NULL, dist = "normal", scale = 1, offset = NULL)
{
	if(is.null(offset)) 
		offset <- numeric(n)

    x0 <- runif(n, 0, 1)
  	x1 <- x0 * 0.7 + runif(n, 0, 0.3)
  	x2 <- runif(n, 0, 1)
  	x3 <- x2 * 0.9 + runif(n, 0, 0.1)
  	
  	f0 <- function(x) 
  	{
    		2 * sin(pi * x)
  	}
  
  	f1 <- function(x)
  	{
    		exp(2 * x)
  	}
  
  	f2 <- function(x)
  	{
  	  	0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
  	}
  
  	f3 <- function(x)
  	{
  		  0 * x
  	}
	
  	# nonparametric component (centered or not)
  
  	f_noncenter <- f0(x0) + f1(x1) + f2(x2)
  	f_smooth <- f_noncenter - mean(f_noncenter)
  	if(!is.null(extra.X))
  	{
		if(is.null(colnames(extra.X))) colnames(extra.X) <- paste("Para.X", 1:ncol(extra.X), sep = "")
  	    f_all <- f_smooth + as.matrix(extra.X) %*% beta + offset
  	}
    else
    {
        f_all <- f_smooth + offset
    }  	
  	if(dist == "normal")
  	{
        y <- rnorm(n, f_all, scale)
  	}
  	else if(dist == "poisson")
  	{
  	    y <- rpois(n, exp(f_all))
  	}
  	else if(dist == "binomialp")
  	{
        g <- binomial(link = "probit")$linkinv(f_all)
        y <- rbinom(n, 1, g)
  	}
  	else if(dist == "binomial")
  	{
        g <- binomial()$linkinv(f_all)
        y <- rbinom(n, 1, g)
  	}
  	else
  	{
        stop("dist not recognised")
  	}

  	if(!is.null(extra.X))
  	{
	## f is linear predictor just for the smoothed part; linear.predictor includes parametric components if required
  	    data <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2, x3 = x3, f = f_smooth, f0 = f0(x0), f1 = f1(x1), f2 = f2(x2), f3 = x3 * 0, linear.predictor = f_all, extra.X)
  	}
  	else
  	{
  	   data <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2, x3 = x3, f = f_smooth, f0 = f0(x0), f1 = f1(x1), f2 = f2(x2), f3 = x3 * 0, linear.predictor = f_all)	
  	}
  	
  	return(data)
}

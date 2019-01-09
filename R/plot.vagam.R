plot.vagam <- function(x, n = 100, alpha = 0.05, rug = TRUE, se = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL, select = NULL, ...)
{
	if(class(x) != "vagam")
	{
		stop("x should be an object of class vagam. Thanks!")
	}
    if(is.null(x$smooth.X))
    {
        stop("Fitted GAM must contain the data. Thanks!")
    }
    
	if(is.null(select)) 
        sel_plots <- 1:ncol(x$smooth.X)
	if(!is.null(select)) 
        sel_plots <- select
	for(k2 in sel_plots) 
	{
		calc_pred <- predict.vagam(object = x, new.smoothX = seq(min(x$smooth.X[,k2]), max(x$smooth.X[,k2]),
								length = n), terms = k2, alpha = alpha)
		
		plot(calc_pred$new.smoothX, calc_pred$prediction, type = "l", xlim = xlim, ylim = ylim, main = main,
                                xlab = ifelse(is.null(xlab), colnames(x$smooth.X)[k2], xlab), 
                                ylab = ifelse(is.null(ylab), paste("Smooth of", colnames(x$smooth.X)[k2]), ylab), ...)
		if(se)
		{
			lines(calc_pred$new.smoothX, calc_pred$lower, lty = 2, ...)
			lines(calc_pred$new.smoothX, calc_pred$upper, lty = 2, ...)			
		}
		if(rug)
		{
			rug(x$smooth.X[,k2], ...)
		}
	}
}

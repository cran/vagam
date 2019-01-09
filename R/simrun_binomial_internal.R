simrun_binomial <- function(t, band_const, choose.n = 100,
num.holdout.points = 10, true.beta = c(-1, 0.5),
CI = 0.95)
{
    set.seed(123 + t)
    choose.k <- band_const * ceiling(choose.n^0.2)
    
    all.results <- matrix(NA, nrow = 9, ncol = 5)
    rownames(all.results) = c("comptime", "mse", "bias-para", "mse-para", "confint-para", "width-para", "width-smooth", "intervalscore-smooth","mse-meanresp")
    colnames(all.results) = c("mgcv-default", "mgcv-psplines", "gamm4", "vagams-unstructuredA", "vagams-bdiagA")
    
    sel.holdout.points <- sample(1:choose.n, num.holdout.points)
    binomial_dat <- gamsim(n = choose.n, dist = "binomial", extra.X = data.frame(intercept = rep(1,choose.n), treatment = rep(c(0,1), each = choose.n/2)), beta = true.beta)
    
    #############################
    ## 1) GAM mgcv with default
    #############################
    
    tic <- proc.time()
    fit.mgcv1 <- gam(y~treatment + s(x0) + s(x1) + s(x2) + s(x3), data = binomial_dat, family = binomial(link = "logit"))
    all.results[1,1] <- (proc.time() - tic)[3]
    all.results[2,1] <-  mean((fit.mgcv1$linear.predictors - binomial_dat$linear.predictor)^2)
    all.results[3,1] <- fit.mgcv1$coefficients[2] - true.beta[2]
    all.results[4,1] <- (fit.mgcv1$coefficients[2] - true.beta[2])^2
    make.ci <- c(summary(fit.mgcv1)$p.coeff[2] - qnorm(0.5 + CI/2) *summary(fit.mgcv1)$se[2],
    summary(fit.mgcv1)$p.coeff[2] + qnorm(0.5 + CI/2) * summary(fit.mgcv1)$se[2])
    all.results[5,1] <- findInterval(true.beta[2], make.ci) == 1
    all.results[6,1] <- diff(make.ci)
    
    holdout.fit <- function(data, holdout.points)
    {
        new.dat <- data[-holdout.points,]
        fit1 <- gam(y ~ treatment + s(x0) + s(x1) + s(x2) + s(x3), family = binomial(link = "logit"), data = new.dat)
        get.pred <- predict.gam(fit1, newdata = data[holdout.points,], se.fit = TRUE)
        get.pred <- list(fit = c(get.pred$fit) - cbind(1,data[holdout.points,"treatment"])%*%fit1$coefficients[1:2], se.fit = c(get.pred$se.fit)) ## Only want smooth fit, so subtract the parametric component out (although se does not remove this!)
        make.ci <- cbind(get.pred$fit - qnorm(0.5 + CI/2) * get.pred$se.fit, get.pred$fit + qnorm(0.5 + CI/2) * get.pred$se.fit)
        
        all.widths <- apply(make.ci,1,diff)
        all.coverage <- sapply(1:length(holdout.points), function(x) findInterval(data$f[x], make.ci[x,]) != 1)
        all.interval.score <- all.widths + (2/(1 - CI))*as.numeric(all.coverage)
        return(cbind(all.widths, all.coverage, all.interval.score))
    }
    
    do.holdfits <- holdout.fit(data = binomial_dat, holdout.points = sel.holdout.points)
    all.results[7,1] <- colMeans(do.holdfits)[1]
    all.results[8,1] <- colMeans(do.holdfits)[3]
    all.results[9,1] <-  mean((binomial()$linkinv(fit.mgcv1$linear.predictors) - binomial()$linkinv(binomial_dat$linear.predictor))^2)
    rm(holdout.fit)
    
    ###############################################
    ## 2) GAM mgcv with P splines and preset knots
    ###############################################
    
    tic <- proc.time()
    fit.mgcv2 <- gam(y ~ treatment + s(x0, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x1, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x2, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x3, bs = "ps", k = choose.k + 2, m  = c(2,1)), data = binomial_dat, family = binomial(link = "logit"))
    all.results[1,2] <- (proc.time() - tic)[3]
    all.results[2,2] <- mean((fit.mgcv2$linear.predictors - binomial_dat$linear.predictor)^2) ## Averaged across choose.n points in dataset
    all.results[3,2] <- fit.mgcv2$coefficients[2] - true.beta[2]
    all.results[4,2] <- (fit.mgcv2$coefficients[2] - true.beta[2])^2
    make.ci <- c(summary(fit.mgcv2)$p.coeff[2] - qnorm(0.5 + CI/2) * summary(fit.mgcv2)$se[2],
    summary(fit.mgcv2)$p.coeff[2] + qnorm(0.5 + CI/2) * summary(fit.mgcv2)$se[2])
    all.results[5,2] <- findInterval(true.beta[2], make.ci) == 1
    all.results[6,2] <- diff(make.ci)
    
    holdout.fit <- function(data, holdout.points) ## Not sure if you are holding out each point at a time, or you are holding out all chosen points in one go!
    {
        new.dat <- data[-holdout.points,]
        fit1 <- gam(y ~ treatment + s(x0, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x1, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x2, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x3, bs = "ps", k = choose.k + 2, m  = c(2,1)), family = binomial(link = "logit"), data = new.dat)
        get.pred <- predict.gam(fit1, newdata = data[holdout.points,], se.fit = TRUE)
        get.pred <- list(fit = c(get.pred$fit) - cbind(1,data[holdout.points,"treatment"])%*%fit1$coefficients[1:2], se.fit = c(get.pred$se.fit)) ## Only want smooth fit, so subtract the parametric component out (although se does not remove this!)
        make.ci <- cbind(get.pred$fit - qnorm(0.5 + CI/2) * get.pred$se.fit, get.pred$fit + qnorm(0.5 + CI/2) * get.pred$se.fit)
        
        all.widths <- apply(make.ci,1,diff)
        all.coverage <- sapply(1:length(holdout.points), function(x) findInterval(data$f[x], make.ci[x,]) != 1)
        all.interval.score <- all.widths + (2/(1 - CI))*as.numeric(all.coverage)
        return(cbind(all.widths, all.coverage, all.interval.score))
    }
    
    do.holdfits <- holdout.fit(data = binomial_dat, holdout.points = sel.holdout.points)
    all.results[7,2] <- colMeans(do.holdfits)[1]
    all.results[8,2] <- colMeans(do.holdfits)[3]
    all.results[9,2] <- mean((binomial()$linkinv(fit.mgcv2$linear.predictors) - binomial()$linkinv(binomial_dat$linear.predictor))^2)
    rm(holdout.fit)
    
    ###############################################
    ## 3) GAMM4 using mixed model parameterization
    ###############################################
    
    tic <- proc.time()
    fit.gamm4 <- gamm4(y ~ treatment + s(x0, bs = "ps", k = choose.k+2, m  = c(2,1)) + s(x1, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x2, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x3, bs = "ps", k = choose.k + 2, m  = c(2,1)), data = binomial_dat, family = binomial(link = "logit"))
    all.results[1,3] <- (proc.time() - tic)[3]
    all.results[2,3] <- mean((fit.gamm4$gam$linear.predictors - binomial_dat$linear.predictor)^2) ## Averaged across choose.n points in dataset
    all.results[3,3] <- fit.gamm4$gam$coefficients[2] - true.beta[2]
    all.results[4,3] <- (fit.gamm4$gam$coefficients[2] - true.beta[2])^2
    make.ci <- c(summary(fit.gamm4$gam)$p.coeff[2] - qnorm(0.5 + CI/2) * summary(fit.gamm4$gam)$se[2],
    summary(fit.gamm4$gam)$p.coeff[2] + qnorm(0.5 + CI/2) * summary(fit.gamm4$gam)$se[2])
    all.results[5,3] <- findInterval(true.beta[2], make.ci) == 1
    all.results[6,3] <- diff(make.ci)
    
    holdout.fit <- function(data, holdout.points)
    {
        new.dat <- data[-holdout.points,]
        fit1 <- gamm4(y ~ treatment + s(x0, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x1, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x2, bs = "ps", k = choose.k + 2, m  = c(2,1)) + s(x3, bs = "ps", k = choose.k + 2, m  = c(2,1)), family = binomial(link = "logit"), data = new.dat)
        get.pred <- predict.gam(fit1$gam, newdata = data[holdout.points,], se.fit = TRUE)
        get.pred <- list(fit = c(get.pred$fit) - cbind(1,data[holdout.points,"treatment"])%*%fit1$gam$coefficients[1:2], se.fit = c(get.pred$se.fit)) ## Only want smooth fit, so subtract the parametric component out (although se does not remove this!)
        make.ci <- cbind(get.pred$fit - qnorm(0.5 + CI/2) * get.pred$se.fit, get.pred$fit + qnorm(0.5 + CI/2) * get.pred$se.fit)
        
        all.widths <- apply(make.ci,1,diff)
        all.coverage <- sapply(1:length(holdout.points), function(x) findInterval(data$f[x], make.ci[x,]) != 1)
        all.interval.score <- all.widths + (2/(1 - CI))*as.numeric(all.coverage)
        return(cbind(all.widths, all.coverage, all.interval.score))
    }
    
    do.holdfits <- try(holdout.fit(data = binomial_dat, holdout.points = sel.holdout.points), silent = TRUE)
    if(!inherits(do.holdfits, "try-error"))
    {
        all.results[7,3] <- colMeans(do.holdfits)[1]
        all.results[8,3] <- colMeans(do.holdfits)[3]
        rm(holdout.fit)
    }
    all.results[9,3] <-  mean((binomial()$linkinv(fit.gamm4$gam$linear.predictors) - binomial()$linkinv(binomial_dat$linear.predictor))^2)
    
    #######################################################
    ## 4) VA with unstructed variational covariance matrix
    #######################################################
    
    tic <- proc.time()
    Rprof()
    #     fit.va1 <- vagam(y = binomial_dat$y, smooth.X = binomial_dat[,c(2:5)], para.X = data.frame(treatment = binomial_dat$treatment), int.knots = choose.k, A.struct = "unstructured", save.data = TRUE, family = binomial(), para.se = TRUE, control = list(eps = 1e-3, maxit = 1000, trace = TRUE, seed_number = t, mc.samps = 4000, pois.step.size = 0.01))
    fit.va1 <- vagam(y = binomial_dat$y, smooth.X = binomial_dat[,c(2:5)], para.X = data.frame(treatment = binomial_dat$treatment), int.knots = choose.k, A.struct = "unstructured", save.data = TRUE, family = binomial(), para.se = TRUE, control = list(eps = 1e-3, maxit = 1000, trace = FALSE, seed_number = t, mc.samps = 4000, pois.step.size = 0.01))
    Rprof(NULL)
    all.results[1,4] <- (proc.time() - tic)[3]
    all.results[2,4] <- mean((fit.va1$linear.predictors - binomial_dat$linear.predictor)^2) ## Averaged across choose.n points in dataset
    all.results[3,4] <- fit.va1$kappa[2] - true.beta[2]
    all.results[4,4] <- (fit.va1$kappa[2] - true.beta[2])^2
    make.ci <- c(fit.va1$para.stat[1,2] - qnorm(0.5 + CI/2) * fit.va1$para.stat[2,2],
    fit.va1$para.stat[1,2] + qnorm(0.5 + CI/2) * fit.va1$para.stat[2,2])
    all.results[5,4] <- findInterval(true.beta[2], make.ci) == 1
    all.results[6,4] <- diff(make.ci)
    
    holdout.fit <- function(data, holdout.points)
    {
        new.dat <- data[-holdout.points,]
        fit1 <- vagam(y = new.dat$y, smooth.X = new.dat[,2:5], para.X = data.frame(treatment = new.dat$treatment), int.knots = choose.k, A.struct = "unstructured", save.data = TRUE, para.se = FALSE, family = binomial(), control = list(eps = 1e-3, maxit = 1000, trace = FALSE, seed_number = t, mc.samps = 4000, pois.step.size = 0.01))
        get.pred <- predict.vagam(fit1, new.smoothX = data[holdout.points,2:5])
        make.ci <- cbind(get.pred$lower.bound, get.pred$upper.bound)
        
        all.widths <- apply(make.ci,1,diff)
        all.coverage <- sapply(1:length(holdout.points), function(x) findInterval(data$f[x], make.ci[x,]) != 1)
        all.interval.score <- all.widths + (2/(1 - CI))*as.numeric(all.coverage)
        return(cbind(all.widths, all.coverage, all.interval.score))
    }
    
    do.holdfits <- holdout.fit(data = binomial_dat, holdout.points = sel.holdout.points)
    all.results[7,4] <- colMeans(do.holdfits)[1]
    all.results[8,4] <- colMeans(do.holdfits)[3]
    all.results[9,4] <-  mean((binomial()$linkinv(fit.va1$linear.predictors) - binomial()$linkinv(binomial_dat$linear.predictor))^2)
    rm(holdout.fit)
    
    ###########################################################
    ## 5) VA with block diagonal variational covariance matrix
    ###########################################################
    
    tic <- proc.time()
    fit.va2 <- vagam(y = binomial_dat$y, smooth.X = binomial_dat[,2:5], para.X = data.frame(treatment = binomial_dat$treatment), int.knots = choose.k, A.struct = "block", save.data = TRUE, family = binomial(), para.se = TRUE, control = list(eps = 1e-3, maxit = 1000, trace = FALSE, seed_number = t, mc.samps = 4000, pois.step.size = 0.01))
    all.results[1,5] <- (proc.time() - tic)[3]
    all.results[2,5] <- mean((fit.va2$linear.predictors - binomial_dat$linear.predictor)^2) ## Averaged across choose.n points in dataset
    all.results[3,5] <- fit.va2$kappa[2] - true.beta[2]
    all.results[4,5] <- (fit.va2$kappa[2] - true.beta[2])^2
    make.ci <- c(fit.va2$para.stat[1,2] - qnorm(0.5 + CI/2) * fit.va2$para.stat[2,2],
    fit.va2$para.stat[1,2] + qnorm(0.5 + CI/2) * fit.va2$para.stat[2,2])
    all.results[5,5] <- findInterval(true.beta[2], make.ci) == 1
    all.results[6,5] <- diff(make.ci)
    
    holdout.fit <- function(data, holdout.points)
    {
        new.dat <- data[-holdout.points,]
        fit1 <- vagam(y = new.dat$y, smooth.X = new.dat[,2:5], para.X = data.frame(treatment = new.dat$treatment), int.knots = choose.k,
        A.struct = "block", save.data = TRUE, para.se = FALSE, family = binomial(), control = list(eps = 1e-3, maxit = 1000, trace = FALSE, seed_number = t, mc.samps = 4000, pois.step.size = 0.01))
        get.pred <- predict.vagam(fit1, new.smoothX = data[holdout.points,2:5])
        make.ci <- cbind(get.pred$lower.bound, get.pred$upper.bound)
        
        all.widths <- apply(make.ci,1,diff)
        all.coverage <- sapply(1:length(holdout.points), function(x) findInterval(data$f[x], make.ci[x,]) != 1)
        all.interval.score <- all.widths + (2/(1 - CI))*as.numeric(all.coverage)
        return(cbind(all.widths, all.coverage, all.interval.score))
    }
    
    do.holdfits <- holdout.fit(data = binomial_dat, holdout.points = sel.holdout.points)
    all.results[7,5] <- colMeans(do.holdfits)[1]
    all.results[8,5] <- colMeans(do.holdfits)[3]
    all.results[9,5] <-  mean((binomial()$linkinv(fit.va2$linear.predictors) - binomial()$linkinv(binomial_dat$linear.predictor))^2)
    return(all.results)
}


'degreedistr' <-
function(web, plot.it=TRUE, pure.call=TRUE, silent=TRUE, level="both", ...){

    # calculates cumulative degree distributions and fits exponential, power law
    # and truncated power law functions to it
    # return fits
    web <- empty(web)
    web <- (web>0)*1 #turns it into a qualitative network
    k <- sum(web) # number of links in network
    S <- sum(dim(web)) # number of species in network
    ddlower <- rowSums(web)
    ddhigher <- colSums(web)

	if (level == "both"){
		lower <- TRUE; higher <- TRUE
	} else {
		lower <- FALSE; higher <- FALSE
		if (level =="lower") lower <- TRUE else higher <- TRUE
	}


    Plo <- sapply(sort(unique(ddlower)), function(x) sum(ddlower>=x))
    Plower <- cbind.data.frame(k=sort(unique(ddlower)), P=Plo/max(Plo))
    if (lower & nrow(Plower) < 5) warning("Too few data points (< 5) for lower trophic level! Fitting makes no sense! The truncated fit is the first to fail because it has one parameter more.")
    Phi <- sapply(sort(unique(ddhigher)), function(x) sum(ddhigher>=x))
    Phigher <- cbind.data.frame(k=sort(unique(ddhigher)), P=Phi/max(Phi))
    if (higher & max(Phigher) < 5) warning("Too few data levels of degrees (< 5) for higher trophic level! Fitting makes no sense! The truncated fit is the first to fail because it has one parameter more.")

    fitdd <- function(...){
    	 # Obviously, an intercept "b" has to be fitted, too.
    	 # This was wrongly omitted in the versions <1.09!
				
        start.trials.b <- c(0.01, .5, 1, 2, 4, 10, 100, 1000)
        start.trials.gamma <- c(0.01, 0.1, 1, 10)
        starts <- expand.grid(start.trials.b, start.trials.gamma)
        # exponential
        for (i in 1:nrow(starts)){
        	EXP <- try(nls(P ~ b*exp(-gamma*k), start=list(gamma=starts[i,2], b=starts[i,2]), ...), silent=silent)
        	if (!inherits(EXP, "try-error")) break;
        }
        # power law
        for (i in 1:nrow(starts)){
   		    PL <- try(nls(P ~ b*k^(-gamma), start=list(gamma=starts[i,2], b=starts[i,1]), ...), silent=silent)
        	if (!inherits(PL, "try-error")) break;
        }
        # truncated power law
        if (!inherits(PL, "try-error")){ # try only if PL converged!
        	for (kx.try in 10^c(-4:4)){
        		TPL <- try(nls(P ~ b*(k^(-gamma))*exp(-k/kx), start=list(gamma=coef(PL)[1], b=coef(PL)[2], kx=kx.try), ...), silent=silent)
        		if (!inherits(TPL, "try-error")) break;
         } 
        } else {TPL <- try(sqrt("w"), silent=TRUE)} #only to produce an error!
        
        list(EXP, PL, TPL)
    }

	# for 3 or fewer data points, fitting is silly. Thus, a first check is performed that returns try-errors if too few data points are found:
	if (lower & length(Plo) < 4){fitl <- list(try(sqrt("w"), silent=TRUE))} else {fitl <- fitdd(data=Plower, nls.control(maxiter=1000))}
	if (higher & length(Phi) < 4){fith <- list(try(sqrt("w"), silent=TRUE))} else {fith <- fitdd(data=Phigher, nls.control(maxiter=1000))}
#was:    fith <- fitdd(data=Phigher, nls.control(maxiter=1000))

	# set NAs right for plotting:
	indexl <- which(sapply(fitl, function(x) !inherits(x, "try-error"))!=0)
    fitnew <- fitl[indexl]
    
    indexh <- which(sapply(fith, function(x) !inherits(x, "try-error"))!=0)
    fithnew <- fith[indexh]


    if (plot.it){
      plotfit <- function(data, fit, ...){
          plot(data$P ~ data$k, log="xy", pch=16, xlab="number of links [k]", ylab="cumulative distribution", ...)
          abline(h=1, lty=2)
          for (i in 1:3){
              if (!inherits(fit[[i]], "try-error")) lines(seq(1, max(ddlower), by=0.1), predict(fit[[i]],
                newdata=data.frame(k=seq(1, max(ddlower), by=0.1))), col=paste("grey", i*20, sep=""), ...)
          }
          
      }
      if (pure.call & level == "both") par(mfrow=c(1,2), mar=c(5,5,4,1)) else par(mar=c(5,5,4,1))
      if (lower) plotfit(data=Plower, fit=fitl, lwd=2, cex=1.8, cex.lab=1.5, main="lower trophic level", ...)
      if (higher) plotfit(data=Phigher, fit=fith, lwd=2, cex=1.8, cex.lab=1.5, main="higher trophic level", ...)
    }

    res.out <- matrix(ncol=5, nrow=3)
    rownames(res.out) <- c("exponential", "power law", "truncated power law")
    colnames(res.out) <- c("Estimate", "Std. Error", "Pr(>|t|)", "R2", "AIC")

	# if all fits are "try-error":
	if (lower & all(sapply(fitl, function(x) inherits(x, "try-error")))) only.na.low <- TRUE else only.na.low <- FALSE
	if (higher & all(sapply(fith, function(x) inherits(x, "try-error")))) only.na.high <- TRUE else only.na.high <- FALSE
	
	# lower trophic level:
    resl.out <- res.out
	if (!only.na.low){ # catch the case that none of the regressions converges:
		resl2 <- t(sapply(fitl[indexl], function(mod) try(c(R2=cor(eval(mod$data)$P, predict(mod)), AIC=AIC(mod)), silent=TRUE )))
		resl.out[indexl, 4:5] <- resl2
		resl1 <- t(sapply(fitl[indexl], function(mod) coef(summary(mod))[1, c(1,2,4)]))
		resl.out[indexl, 1:3] <- resl1
	}

    # higher trophic level:
    resh.out <- res.out
    if (!only.na.high){ # catch the case that none of the regressions converges:
		resh2 <- t(sapply(fith[indexh], function(mod) try(c(R2=cor(eval(mod$data)$P, predict(mod)), AIC=AIC(mod)), silent=TRUE )))
		resh.out[indexh, 4:5] <- resh2
		resh1 <- t(sapply(fith[indexh], function(mod) coef(summary(mod))[1, c(1,2,4)]))
		resh.out[indexh, 1:3] <- resh1
	}

    if (level == "both"){
    	 out <- list("lower level dd fits"=resl.out, "higher level dd fits"=resh.out)
 	} else {
 		out <- if(lower) list("lower level dd fits"=resl.out) else "higher level dd fits"=resh.out
 	}
 	return(out)

}

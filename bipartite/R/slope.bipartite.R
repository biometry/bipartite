`slope.bipartite` <-
function(object, plot.it=TRUE, ...){
    # function to calculate the slope of the extinction curve
    # of a random extinction result and plot it (optionally)
    # object object of class bipartite from second.extinction analysis
    # plot.it   logical: should data and fit be plotted (defaults to TRUE)
    # ...       arguments passed on to plot (not to curve!)

    if (class(object)!="bipartite") stop("This function cannot be meaningfully applied to objects of this class.")
    if (is.list(object)) stop("You seem to have computed extinction slopes for 'both' trophic levels using 'random' extermination and requesting the details of the sequences. Since this leads to a list of different-length sequences, you cannot use this function. \n\nPlease use method 'random' for each level separately OR set details to FALSE. \n\nWe apologise for the inconvenience.")
    N <- colSums(object)

    if (all(object[-nrow(object), 2] == 1)) y <- -object[, 3] else y <- -object[, 2] #selects the correct column
    #y <- 100 - (sum(y)-cumsum(y))/sum(y) * 100 #ranged between 0 and 100

    y <- (sum(y)-cumsum(y))/sum(y)  #ranged between 0 and 1
    x <- (object[,"no"] / max(object[,"no"]))  #ranges x between 0 and 1

#    fit<- nls( y ~ 1 - b*x^a, start=list(a=2, b=1), lower=c(-1, 0.001), upper=c(500, 500), algorithm="port")
    fit <- try(nls(y ~ 1 - x^a, start=list(a=1)), silent=TRUE)
        if (class(fit)=="try-error") fit <- nls( (y+rnorm(length(y), s=0.01)) ~ 1 - x^a, start=list(a=1))

    if(plot.it)
    {
        par(mar=c(5, 5, 1, 1))
        plot(x, y, xlab="proportion of primary extinctions", ylab="proportion of species in other trophic level still alive",
            axes=TRUE, type="n", cex.lab=1.5)
        legend("bottomleft", legend=paste("killed: ", attr(object, "exterminated"), ""), cex=2, bty="n")
        abline(h=1)
        abline(v=1)
        points(x, y, ...)
        lines(seq(0, 1, 0.1), predict(fit, newdata=data.frame(x=seq(0, 1, 0.1))), col="red", lwd=2)
    }

    return(c("exponent"=as.numeric(coef(fit)[1])))#, "coefficient"=as.numeric(coef(fit)[2])))
}


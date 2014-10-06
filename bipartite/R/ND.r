# normalised degree
ND <- function(web, normalised=TRUE){
    # calculates the (normalised) degree of a species
    # by Carsten F. Dormann, 14 Dec 2010
    web <- (web > 0) * 1
    k <- sum(web)
    dlower <- rowSums(web)
    dhigher <- colSums(web)  
    Nlow <- Nhigh <- 2 # effectively unnormalised
    if (normalised){
      Nlow <- length(dhigher)
      Nhigh <- length(dlower)
    }
    low <- dlower/Nlow; names(low) <- rownames(web)
    high <- dhigher/Nhigh; names(high) <- colnames(web)
    list("lower"=low, "higher"=high)
}

CC <- function(web, cmode="suminvundir", rescale=TRUE, weighted=TRUE, ...){
    # closeness centrality as used in Gonzales et al. (2009)
    # by Carsten F. Dormann, 14 Dec 2010
    # ... options passed on to closeness in package sna
    # uses a version to calculate centrality that allows for disconnected graphs
    wh <- as.one.mode(web, project="higher", weighted=weighted)
    cch <- closeness(wh, cmode=cmode, rescale=rescale, ...)

    wl <- as.one.mode(web, project="lower", weighted=weighted)
    ccl <- closeness(wl, cmode=cmode, rescale=rescale, ...)

    list("lower"=ccl, "higher"=cch)
}

BC <- function(web, rescale=TRUE, cmode="undirected", weighted=TRUE, ...){
    # betweenness centrality as used in Gonzales et al. (2009)
    # by Carsten F. Dormann, 14 Dec 2010
    # ... options passed on to closeness in package sna; particularly: rescale=TRUE!
    wh <- as.one.mode(web, project="higher", weighted=weighted)
    bch <- betweenness(wh, rescale=FALSE, cmode=cmode, ...)
    if (rescale & sum(bch != 0)) bch <- bch/sum(bch, na.rm=TRUE)

    wl <- as.one.mode(web, project="lower", weighted=weighted)
    bcl <- betweenness(wl, rescale=FALSE, cmode=cmode, ...)
    if (rescale & sum(bcl != 0)) bcl <- bcl/sum(bcl, na.rm=TRUE)

    list("lower"=bcl, "higher"=bch)
}



## example:
#data(olesen2002flores)
#(ndi <- ND(olesen2002flores))
#(cci <- CC(olesen2002flores))
#(bci <- BC(olesen2002flores))
#
#cor.test(bci[[1]], ndi[[1]], method="spear") # 0.779
#cor.test(cci[[1]], ndi[[1]], method="spear") # 0.826
#
#cor.test(bci[[2]], ndi[[2]], method="spear") # 0.992
#cor.test(cci[[2]], ndi[[2]], method="spear") # 0.919
#
### PLANTS:
#bc <- bci[[1]]
#cc <- cci[[1]]
#nd <- ndi[[1]]
## CC:
#summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) # lower RSE
#summary(nls(cc ~ c*nd^d, start=list(c=0.02,d=2))) 
## BC:
#summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
#summary(nls(bc ~ c*nd^d, start=list(c=0.02,d=2))) # lower RSE
#
### ANIMALS:
#bc <- bci[[2]]
#cc <- cci[[2]]
#nd <- ndi[[2]]
## CC:
#summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) 
#summary(nls(cc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
## BC:
#summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
#summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
#
#
##see also, for whole web measures:
#centralization(web, "degree")
#centralization(web, "betweenness")
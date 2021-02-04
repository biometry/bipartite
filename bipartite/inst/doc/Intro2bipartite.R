### R code from vignette source 'Intro2bipartite.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(knitr)
opts_chunk$set(fig.path='figure/twocolumn-', fig.align='center', fig.show='hold', cache=TRUE, fig.width=5, fig.height=5, fig.show='hold', cache=TRUE, tidy=T, tidy.opts=list(width.cutoff=70))
render_listings()


###################################################
### code chunk number 2: Intro2bipartite.Rnw:153-154
###################################################
library(bipartite)


###################################################
### code chunk number 3: Intro2bipartite.Rnw:162-165 (eval = FALSE)
###################################################
## par(xpd=T)
## plotweb(motten1982)
## visweb(motten1982)


###################################################
### code chunk number 4: Intro2bipartite.Rnw:189-190 (eval = FALSE)
###################################################
## plotPAC(PAC(motten1982), outby=0.9)


###################################################
### code chunk number 5: Intro2bipartite.Rnw:205-207 (eval = FALSE)
###################################################
## mod <- computeModules(motten1982)
## plotModuleWeb(mod)


###################################################
### code chunk number 6: Intro2bipartite.Rnw:220-227 (eval = FALSE)
###################################################
## par(mfrow=c(1,2), xpd=T)
## gplot(as.one.mode(motten1982, project="higher"), 
##  label=colnames(motten1982), gmode="graph", 
## label.cex=0.6, vertex.cex=2)
## gplot(as.one.mode(motten1982, project="lower"), 
## 	label=rownames(motten1982), gmode="graph", 
## 	label.cex=0.6, vertex.cex=2, vertex.col="green")


###################################################
### code chunk number 7: Intro2bipartite.Rnw:290-292
###################################################
networklevel(bezerra2009, index=c("ISA", "weighted NODF", "Fisher alpha"), 
             SAmethod="log")


###################################################
### code chunk number 8: Intro2bipartite.Rnw:322-324
###################################################
grouplevel(bezerra2009, level="both", index=c("mean number of links", "weighted 
     cluster coefficient", "effective partners", "niche overlap"), dist="bray")


###################################################
### code chunk number 9: Intro2bipartite.Rnw:337-338
###################################################
str(linklevel(bezerra2009, index=c("dependence", "endpoint")))


###################################################
### code chunk number 10: Intro2bipartite.Rnw:369-371
###################################################
specieslevel(bezerra2009, level="lower", index=c("normalised degree", "PDI", 
      "effective partners"), PDI.normalise=F)


###################################################
### code chunk number 11: networkPCA (eval = FALSE)
###################################################
## web.names <- data(package="bipartite")$results[,3]
## data(list=web.names) #loads all webs
## # the next step takes around 10 minutes:
## netw.indic.webs <- t(sapply(web.names, function(x) networklevel(get(x), 
##               index="ALLBUTDD")))


###################################################
### code chunk number 12: load external data
###################################################
#load("/Volumes/Macintosh HD/Users/Carsten/Data/aktuell/bipartite/bipartite/vignettes/figures/netw.indic.webs.Rdata")
load("./figures/netw.indic.webs.Rdata") # loads the files stored above; avoids re-running this time-consuming analysis!


###################################################
### code chunk number 13: Intro2bipartite.Rnw:422-424 (eval = FALSE)
###################################################
## PCA.out <- prcomp(netw.indic.webs[,-5], scale.=T)
## biplot(PCA.out, xpd=T, las=1)


###################################################
### code chunk number 14: Intro2bipartite.Rnw:430-431
###################################################
summary(PCA.out)


###################################################
### code chunk number 15: Intro2bipartite.Rnw:435-436
###################################################
PCA.out <- prcomp(netw.indic.webs[,-5], scale.=T)


###################################################
### code chunk number 16: Intro2bipartite.Rnw:438-439
###################################################
round(PCA.out$rotation[, 1:4], 3)


###################################################
### code chunk number 17: varclus (eval = FALSE)
###################################################
## library(Hmisc)
## plot(varclus(netw.indic.webs), cex=0.8)
## abline(h=0.5, lty=2, col="grey")


###################################################
### code chunk number 18: Intro2bipartite.Rnw:537-543 (eval = FALSE)
###################################################
## data(Safariland)
## Iobs <- nestednodf(Safariland)$statistic[3]
## nulls <- nullmodel(web=Safariland, N=1000, method='r2d') # takes a while!
## Inulls <- sapply(nulls, function(x) nestednodf(x)$statistic[3])
## plot(density(Inulls), xlim=c(0, 100), lwd=2, main="NODF")
## abline(v=Iobs, col="red", lwd=2)


###################################################
### code chunk number 19: Intro2bipartite.Rnw:575-584
###################################################
weblist <- lapply(c("Safariland", "vazarr", "vazllao", "vazcer", "vazmasc", 
                       "vazmasnc", "vazquec", "vazquenc"), get)
# Write a function to compute the desired statistic, e.g. the difference 
# between grazed and ungrazed:
meandiff <- function(webs){
   obs <- sapply(webs, networklevel, index="linkage density")  
   mean(obs[1:4] - obs[5:8])
}
(observed <- meandiff(weblist))


###################################################
### code chunk number 20: Intro2bipartite.Rnw:587-589
###################################################
nulllist <- lapply(weblist, nullmodel, N=1, method="r2d")
meandiff(weblist)


###################################################
### code chunk number 21: Intro2bipartite.Rnw:592-597
###################################################
res <- 1:5000
for (i in 1:5000){ # takes a few minutes !!
   nulllist <- sapply(weblist, nullmodel, N=1, method="r2d")
   res[i] <- meandiff(nulllist)  
}


###################################################
### code chunk number 22: Intro2bipartite.Rnw:608-612 (eval = FALSE)
###################################################
## hist(res, xlim=c(-0.3, 0.3), border="white", col="grey")
## abline(v=observed, col="red", lwd=2)
## # compute p-value as proportion smaller or than observed
## sum(res < observed)/length(res) * 2 # *2 for two-tailed test


###################################################
### code chunk number 23: Intro2bipartite.Rnw:704-706 (eval = FALSE)
###################################################
## library(devtools)
## install_github(rep="pedroj/bipartite_plots")



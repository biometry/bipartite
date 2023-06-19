### R code from vignette source 'Intro2bipartite.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(knitr)
opts_chunk$set(fig.path='figures/twocolumn-', fig.align='center', fig.show='hold', cache=TRUE, fig.width=5, fig.height=5, fig.show='hold', cache=TRUE, tidy=F, tidy.opts=list(width.cutoff=70))
#render_listings()


###################################################
### code chunk number 2: Intro2bipartite.Rnw:155-156
###################################################
library(bipartite)


###################################################
### code chunk number 3: Intro2bipartite.Rnw:164-167 (eval = FALSE)
###################################################
## par(xpd=T)
## plotweb(motten1982)
## visweb(motten1982)


###################################################
### code chunk number 4: Intro2bipartite.Rnw:191-192 (eval = FALSE)
###################################################
## plotPAC(PAC(motten1982), outby=0.9)


###################################################
### code chunk number 5: Intro2bipartite.Rnw:207-209 (eval = FALSE)
###################################################
## mod <- computeModules(motten1982)
## plotModuleWeb(mod)


###################################################
### code chunk number 6: Intro2bipartite.Rnw:222-229 (eval = FALSE)
###################################################
## par(mfrow=c(1,2), xpd=T)
## gplot(as.one.mode(motten1982, project="higher"), 
##  label=colnames(motten1982), gmode="graph", 
## label.cex=0.6, vertex.cex=2)
## gplot(as.one.mode(motten1982, project="lower"), 
## 	label=rownames(motten1982), gmode="graph", 
## 	label.cex=0.6, vertex.cex=2, vertex.col="green")


###################################################
### code chunk number 7: Intro2bipartite.Rnw:292-294
###################################################
networklevel(bezerra2009, index=c("ISA", "weighted NODF", "Fisher alpha"), 
             SAmethod="log")


###################################################
### code chunk number 8: Intro2bipartite.Rnw:324-326
###################################################
grouplevel(bezerra2009, level="both", index=c("mean number of links", "weighted 
     cluster coefficient", "effective partners", "niche overlap"), dist="bray")


###################################################
### code chunk number 9: Intro2bipartite.Rnw:339-340
###################################################
str(linklevel(bezerra2009, index=c("dependence", "endpoint")))


###################################################
### code chunk number 10: Intro2bipartite.Rnw:371-373
###################################################
specieslevel(bezerra2009, level="lower", index=c("normalised degree", "PDI", 
      "effective partners"), PDI.normalise=F)


###################################################
### code chunk number 11: betweenPlot
###################################################
data(Safariland)

# plot the one-mode projection for the lower level:
set.seed(4) # don't ask me why gplot is stochastic ...
par(xpd=T, mar=c(0,6,0,6))
gplot(as.one.mode(Safariland, project="lower"), label=rownames(Safariland))


###################################################
### code chunk number 12: Intro2bipartite.Rnw:396-398
###################################################
# convert matrix into one-mode edgelist:
SafPlantsEL <- as.tnet(as.one.mode(Safariland, project="lower"))


###################################################
### code chunk number 13: Intro2bipartite.Rnw:400-408 (eval = FALSE)
###################################################
## # compute betweenness:
## tnet::betweenness_w(SafPlantsEL) # 6
## bipartite::BC(Safariland, rescale=F)$lower # 6
## DiagrammeR::get_betweenness(DiagrammeR::from_igraph(tnet_igraph(SafPlantsEL))) # 6
## igraph::betweenness(tnet_igraph(SafPlantsEL), cutoff=9) # 6
## influenceR::betweenness(tnet_igraph(SafPlantsEL)) # 12
## sna::betweenness(SafPlantsEL) # length 36!!
## sna::betweenness(as.matrix(SafPlantsEL)) # length 27!!


###################################################
### code chunk number 14: Intro2bipartite.Rnw:412-419 (eval = FALSE)
###################################################
## btws <- cbind("t:betw"=tnet::betweenness_w(SafPlantsEL)[,2],
##     "b:BC"=bipartite::BC(Safariland, rescale=F)$lower,
##     "DR:betw"=DiagrammeR::get_betweenness(DiagrammeR::from_igraph(tnet_igraph(SafPlantsEL)))[,2],
##     "i:estbetw"=igraph::betweenness(tnet_igraph(SafPlantsEL), cutoff=9),
##     "inf:betw"=influenceR::betweenness(tnet_igraph(SafPlantsEL))/2 )
## rownames(btws) <- rownames(Safariland)
## btws


###################################################
### code chunk number 15: networkPCA (eval = FALSE)
###################################################
## web.names <- data(package="bipartite")$results[,3]
## data(list=web.names) #loads all webs
## # the next step takes around 10 minutes:
## netw.indic.webs <- t(sapply(web.names, function(x) networklevel(get(x), 
##               index="ALLBUTDD")))


###################################################
### code chunk number 16: load external data
###################################################
#load("/Volumes/Macintosh HD/Users/Carsten/Data/aktuell/bipartite/bipartite/vignettes/figures/netw.indic.webs.Rdata")
load("./figures/netw.indic.webs.Rdata") # loads the files stored above; avoids re-running this time-consuming analysis!


###################################################
### code chunk number 17: Intro2bipartite.Rnw:496-498 (eval = FALSE)
###################################################
## PCA.out <- prcomp(netw.indic.webs[,-5], scale.=T)
## biplot(PCA.out, xpd=T, las=1)


###################################################
### code chunk number 18: Intro2bipartite.Rnw:504-505
###################################################
summary(PCA.out)


###################################################
### code chunk number 19: Intro2bipartite.Rnw:509-510
###################################################
PCA.out <- prcomp(netw.indic.webs[,-5], scale.=T)


###################################################
### code chunk number 20: Intro2bipartite.Rnw:512-513
###################################################
round(PCA.out$rotation[, 1:4], 3)


###################################################
### code chunk number 21: varclus (eval = FALSE)
###################################################
## library(Hmisc)
## plot(varclus(netw.indic.webs), cex=0.8)
## abline(h=0.5, lty=2, col="grey")


###################################################
### code chunk number 22: Intro2bipartite.Rnw:611-617 (eval = FALSE)
###################################################
## data(Safariland)
## Iobs <- nestednodf(Safariland)$statistic[3]
## nulls <- nullmodel(web=Safariland, N=1000, method='r2d') # takes a while!
## Inulls <- sapply(nulls, function(x) nestednodf(x)$statistic[3])
## plot(density(Inulls), xlim=c(0, 100), lwd=2, main="NODF")
## abline(v=Iobs, col="red", lwd=2)


###################################################
### code chunk number 23: Intro2bipartite.Rnw:649-658
###################################################
weblist <- lapply(c("Safariland", "vazarr", "vazllao", "vazcer", "vazmasc", 
                       "vazmasnc", "vazquec", "vazquenc"), get)
# Write a function to compute the desired statistic, e.g. the difference 
# between grazed and ungrazed:
meandiff <- function(webs){
   obs <- sapply(webs, networklevel, index="linkage density")  
   mean(obs[1:4]) - mean(obs[5:8])
}
(observed <- meandiff(weblist))


###################################################
### code chunk number 24: Intro2bipartite.Rnw:661-663
###################################################
nulllist <- lapply(weblist, nullmodel, N=1, method="r2d")
meandiff(weblist)


###################################################
### code chunk number 25: Intro2bipartite.Rnw:666-671
###################################################
res <- 1:5000
for (i in 1:5000){ # takes a few minutes !!
   nulllist <- sapply(weblist, nullmodel, N=1, method="r2d")
   res[i] <- meandiff(nulllist)  
}


###################################################
### code chunk number 26: Intro2bipartite.Rnw:682-686 (eval = FALSE)
###################################################
## hist(res, xlim=c(-0.3, 0.3), border="white", col="grey")
## abline(v=observed, col="red", lwd=2)
## # compute p-value as proportion smaller or than observed
## sum(res < observed)/length(res) * 2 # *2 for two-tailed test


###################################################
### code chunk number 27: Intro2bipartite.Rnw:780-782 (eval = FALSE)
###################################################
## library(devtools)
## install_github(rep="pedroj/bipartite_plots")



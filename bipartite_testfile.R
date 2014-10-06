## bipartite test file ##

## CHECK R-CONSOLE FOR RED TEXT INDICATING ERRORS!! ##



## run this file after every change in bipartite before submitting it to CRAN!!
## compile and install bipartite before running this! Do not run on R-functions directly (see bottom!)
# command line, in the respective folder:
# R CMD build bipartite
# R CMD install bipartite_2.03.tar.gz
require(bipartite)

# lazy load data does not require data to be loaded via "data(.)"!
#  as.one.mode
image(as.one.mode(Safariland))
visweb(as.one.mode(Safariland, project="lower", fill=NA), NA.col="green") #NA.col should have no effect and cause no problem!
# check fill=NA and visweb's NA.col:
visweb(as.one.mode(vazquenc, fill=NA), NA.col="green") #slow!

# C.score
C.score(Safariland, FUN=sd, normalise=T) # check that "any" functions works!

# compart
compart(Safariland) # a comparted network
compart(bezerra2009) # an uncomparted network

# computeModules
## a lot to test here! let's start with the problem of calling computeModules twice in a row:
comp1 <- computeModules(vazquenc)
comp2 <- computeModules(vazquenc)
# to be continued ...

# czvalues
czvalues(comp1)

# degreedistr
degreedistr(web=memmott1999)
degreedistr(web=memmott1999, pure.call=F)
degreedistr(web=memmott1999, pure.call=F, level="lower")
degreedistr(web=memmott1999, pure.call=T, level="lower")
degreedistr(web=memmott1999, pure.call=T, level="higher", las=1)

# dfun
set.seed(2)
dfun(vazquenc, abuns=runif(24, 1,9)) # checks whether external abundances are accepted
dfun(t(vazquenc), abuns=runif(7, 1,9)) # checks whether external abundances are accepted

# discrepancy
discrepancy(vazquenc)
#nested(vazquenc, method="discrepancy2")
## Not nice, this; it's a namespace issue; somehow permute::allPerms is not available to nesteddisc, forcing me to load all of vegan here! Don't know what's going wrong.
## OUTCOMMENTED because otherwise vegan would be loaded (through nested)! That would obviously affect the check of all following functions!

# empty
vazquenc[,3] <- 0
empty(vazquenc, count=TRUE)
rm(vazquenc)

# endpoint
endpoint(vazquenc)

# extinction: see second.extinct

# fc
fc(t(Safariland), dist="canberra", method="complete")
fc(t(Safariland), dist="euc", method="single")

# frame2webs ...

# genweb
genweb(10, 20, 4)

# grouplevel
grouplevel(vazquenc)
grouplevel(vazquenc, level="higher")
grouplevel(vazquenc, level="lower", weighted=FALSE)
grouplevel(vazquenc, level="lower", weighted=TRUE, index=c("niche overlap"))
grouplevel(vazquenc, level="higher", weighted=TRUE, index=c("niche overlap", "fd", "weighted cluster coefficient"))

# H2fun
set.seed(2)
hist(sapply(nullmodel(vazquenc, 50),  function(x) H2fun(x)[[1]])) # should be close to 0
H2fun(matrix(runif(100, 2, 100), 10, 10), H2_integer=FALSE)

# linklevel
str(linklevel(memmott1999, index=c("dependence", "endpoint")))

# listModuleInformation
if (!exists("comp1")) comp1 <- computeModules(vazquenc)
listModuleInformation(comp1)

# ND, BC, CC
ND(vazquenc)
BC(vazquenc)
BC(vazquenc, rescale=FALSE, weighted=FALSE)
CC(vazquenc)

# nested
#nested(vazquenc, method="ALL")
## OUTCOMMENTED because otherwise vegan would be loaded (through nested)! That would obviously affect the check of all following functions!

# nestedness
nestedness(Safariland, n.nulls=20)[c(4, 9:20)]
nestedness(Safariland, null.models=FALSE)$temperature

# nestedrank
nestedrank(Safariland, normalise=F)
nestedrank(Safariland, weighted=F, normalise=F)
nestedrank(Safariland, method="binmatnest")
nestedrank(Safariland, method="wine")
nestedrank(Safariland, method="wine", weighted=FALSE)
nestedrank(Safariland, method="sort")
nestedrank(Safariland, method="otto")


# networklevel
networklevel(vazquenc)
networklevel(vazquenc, index="networklevel")
networklevel(vazquenc, index="networklevel", weighted=FALSE)
networklevel(vazquenc, index="quantitative")
networklevel(vazquenc, index="H2")
networklevel(vazquenc, index=c("H2", "cluster coefficient"))
networklevel(vazquenc, legacy=TRUE)
# check that double calls to an index do work:
networklevel(matrix(rpois(16,4),nrow=4),c("H2","H2"))

# nodespec
nodespec(Safariland)
nodespec(Safariland, inf.replace=Inf)

# npartite ...

# null.distr
null.distr(2, vazquenc)
null.distr(2, vazquenc>0, distr="negative binomial") # tests whether it works with binary data

# null.t.test
null.t.test(vazquenc, N=30, index="connectance") # does it work with one index? 
null.t.test(vazquenc, N=30, index=c("connectance", "linkage density")) # does it work with one index?

# nullmodel
lapply(1:5, function(x) nullmodel(vazquenc, N=2, method=x))
lapply(c("r2d", "swap.web", "vaznull", "shuffle.web", "mgen"), function(x) nullmodel(vazquenc, N=2, method=x))

# PAC
PAC(vazquenc)
PAC(vazquenc>0)

# PDI
PDI(vazquenc, normalise=FALSE, log=TRUE)
PDI(vazquenc>0)

# plotModuleWeb
if (!exists("comp1")) comp1 <- computeModules(vazquenc)
plotModuleWeb(comp1)
plotModuleWeb(comp1, plotModules = TRUE, rank = TRUE, weighted = TRUE, displayAlabels = TRUE, displayBlabels = TRUE, labsize = 0.6, xlabel = "", ylabel = "", square.border = "lightgreen", fromDepth = 0, upToDepth = -1)

# plotPAC
plotPAC(PAC(vazquenc))
plotPAC(PAC(vazquenc), scaling = 2, plot.scale = 1.5, fill.col = rgb(0.2, 0.3, 0.4, 0.5), arrow.col = rgb(0.4, 0.3, 0.2, 0.5), outby = 0.5, text = TRUE, circles = TRUE, radius = 1.5)

# plotweb
plotweb(vazquenc)
# test of Jochen's workaround:
plotweb(vazquenc,abuns.type='independent')
plotweb(Safariland, abuns.type='independent')
plotweb(Safariland, abuns.type='independent',arrow="down.center")
plotweb(Safariland, abuns.type='additional',arrow="down.center")
# an example set with abundances
myabuns.low <- rowSums(Safariland)
myabuns.low[rownames(vazquenc)] <- myabuns.low[rownames(vazquenc)] + rowSums(vazquenc)          # "total abundances" (might be independent)
myabuns.low.unused <- myabuns.low - rowSums(Safariland)[names(myabuns.low)]                     # "unused abundances"
plotweb(Safariland, abuns.type='independent',arrow="down.center",low.abun=myabuns.low)          # that's correct
plotweb(Safariland, abuns.type='additional',arrow="down.center",low.abun=myabuns.low)           # that's wrong
plotweb(Safariland, abuns.type='additional',arrow="down.center",low.abun=myabuns.low.unused)    # that's how it should look like for the additional case
plotweb(Safariland, abuns.type='independent',arrow="down.center",low.abun=myabuns.low*0.2)      # still works
plotweb(Safariland, abuns.type='additional',arrow="down.center",low.abun=myabuns.low*0.2)       # always shows marginals as abundances
plotweb(Safariland, abuns.type='independent',arrow="no",low.abun=myabuns.low*0.2)               # currently gives a warning, but the ratio between upper width and lower width could actually tell you something about preferences!




# plotweb2 ...

# printoutModuleInformation
if (!exists("comp1")) comp1 <- computeModules(vazquenc)
printoutModuleInformation(comp1)

# robustness

# second.extinct
bs <- second.extinct(Safariland, method="random", participant="both", details=T) 
slope.bipartite(bs) # should return an error with an explanation
bs <- second.extinct(Safariland, method="random", participant="both", details=F) 
slope.bipartite(bs) # should work


# slope.bipartite

# sortweb

# specieslevel
web <- matrix(c(0,0,0,0,0,0,0,2,0,0,1,3,1,0,5,8),nrow=4) # check for warnings (see email by Jochen 7.3.2014)
specieslevel(bezerra2009)
specieslevel(bezerra2009, index="betweenness") # lower level is all NaN (I guess because bezerra is nearly complete as one-mode?
specieslevel(Safariland, index="betweenness")
specieslevel(Safariland, index="closeness")
specieslevel(bezerra2009, level="higher", index=c("proportional similarity", "proportional generality"))
specieslevel(Safariland, index="nestedrank", nested.weighted=T)
specieslevel(Safariland, index="nestedrank", nested.weighted=F)
specieslevel(Safariland, index="nestedrank", nested.weighted=T, nested.normalise=F)
specieslevel(Safariland, index="nestedrank", nested.weighted=F, nested.normalise=F)
specieslevel(Safariland, index="nestedrank", nested.weighted=F, nested.normalise=T, nested.method="wine")
specieslevel(Safariland, index="nestedrank", nested.weighted=F, nested.normalise=F, nested.method="sort")

# from Jochen's reported errors in specieslevel:
specieslevel(matrix(rpois(16,3),nrow=4),index="ALLBUTD")
specieslevel(matrix(c(0,0,5,0,0,5,0,0,5,0,0,5),nrow=3,byrow=T),index="closeness", level="higher")
specieslevel(matrix(c(0,0,5,0,0,5,0,0,5,0,0,5),nrow=3,byrow=T),index="closeness", level="lower")
specieslevel(matrix(c(0,0,5,0,0,5,0,0,5,0,0,5),nrow=3,byrow=T),index="betweenness", level="lower")
specieslevel(matrix(c(0,0,5,0,0,5,0,0,5,0,0,5),nrow=3,byrow=T),index="betweenness", level="higher")


specieslevel(vazquec, index="interaction push pull") # from Natacha Chacoff's error report
# strength

# swap.web

# togetherness

# V.ratio

# visweb

# web2edges

# webs2array

# wine



## ----------some not so bright ideas ---------##
###source function to make sure we use the most recent ones:
# NOT SO BRIGHT! because the namespace will import specified functions. When we now load the functions via source, these imported functions are NOT available, unless we explicitly load the package, too:
#setwd("/Users/cdormann/Data/aktuell/BESS/R_FoodWeb/")
#rfiles <- list.files("bipartite/R")
## source files, excluding zzz.r
#for (i in seq_along(rfiles)[-length(rfiles)]) source(paste0("bipartite/R/", rfiles[i]))
# C.score(Safariland) # error: could not find "designdist"

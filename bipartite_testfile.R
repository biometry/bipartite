## bipartite test file ##

## CHECK R-CONSOLE FOR RED TEXT INDICATING ERRORS!! ##


# run 
# tools::compactPDF("/Users/Carsten/Data/aktuell/bipartite/bipartite/inst/doc/Dormann2011NetworkBiology.pdf", gs_quality = "ebook") 
# on the large PDFs in inst/ !!
# check link to external functions, use syntax of Hadley's texinfo page





## run this file after every change in bipartite before submitting it to CRAN!!
## compile and install bipartite before running this! Do not run on R-functions directly (see bottom!)
# command line, in the respective folder:
# R CMD build bipartite
# R CMD install bipartite_2.10.tar.gz
library(bipartite)

# lazy load data does not require data to be loaded via "data(.)"!
#  as.one.mode
image(as.one.mode(Safariland))
visweb(as.one.mode(Safariland, project="lower", fill=NA), NA.col="green") #NA.col should have no effect and cause no problem!
# check fill=NA and visweb's NA.col:
visweb(as.one.mode(vazquenc, fill=NA), NA.col="green") #slow!

# C.score
C.score(Safariland, FUN=sd, normalise=T) # check that "any" functions works!
C.score(t(Safariland))
m <- matrix(c(1,0, 0,1, 1,0, 0,1, 1,0, 0,1, 1,0), 7,2,TRUE)
C.score(t(m), normalise=FALSE)
C.score(m, normalise=FALSE)
C.score(t(m), normalise=TRUE) # should be 0.57, since half of the species have non-checkerboard distribution
C.score(m, normalise=TRUE) # should be 1!
n <- matrix(c(1,0,0, 0,1,0, 0,0,1), 3,3)
C.score(t(n), normalise=T)
system.time(C.score(memmott1999)) # 3.6s

# CC (closeness centrality)
CC(Safariland)
CC(Safariland, rescale=F)
CC(Safariland, weighted=F) # weighed=F same as T!
CC(Safariland, rescale=F, weighted=F) # weighed=F same as T!
wf <- as.one.mode(Safariland, project="lower", weighted=F)
wt <- as.one.mode(Safariland, project="lower", weighted=T)
closeness(wf, cmode="suminvundir")
closeness(wt, cmode="suminvundir") # makes no difference!
specieslevel(Safariland, index="closeness", rescale=T) # same as CC(. , rescale=F)



# compart
compart(Safariland) # a comparted network
compart(bezerra2009) # an uncomparted network

# computeModules
## a lot to test here! let's start with the problem of calling computeModules twice in a row:
comp1 <- computeModules(vazquenc)
comp2 <- computeModules(vazquenc, forceLPA=TRUE)
comp3 <- computeModules(vazquenc, method="DormannStrauss")
comp3 <- computeModules(vazquenc, method="DormannStrauss")
plotModuleWeb(comp1)
plotModuleWeb(comp2)
plotModuleWeb(comp3, weighted=F)
# to be continued ...
web <- matrix(c(0,1), 3, 3)
web[1,] <- 1
computeModules(web) # test ability to remove all-1s
# check it works with fully connecte network (error fixed in 2.09):
web <- matrix(runif(150, 0.1, 10), 30, 50)
computeModules(web)

# czvalues
czvalues(comp1)
czvalues(comp3)

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

# DIRT_LPA_wb_plus
res <- DIRT_LPA_wb_plus(Safariland, mini=3, reps=20)
mod <- convert2moduleWeb(Safariland, res)
plotModuleWeb(mod)

# discrepancy
discrepancy(vazquenc)
#nested(vazquenc, method="discrepancy2")
## Not nice, this; it's a namespace issue; somehow permute::allPerms is not available to nesteddisc, forcing me to load all of vegan here! Don't know what's going wrong.
## OUTCOMMENTED because otherwise vegan would be loaded (through nested)! That would obviously affect the check of all following functions!

# Check that this works, because vegan is being restructured and now has a function "nullmodel"!
nulls <- simulate(vegan::nullmodel(Safariland, method="quasiswap"), nsim = 10)
apply(nulls, 3, discrepancy)

# empty
vazquenc[,3] <- 0
empty(vazquenc, count=TRUE)
rm(vazquenc)
M <- matrix(c(0, 0, 1, 0, 1), nrow=1)
empty(M) # should be matrix with 1, 1!

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

# mgen
obs.mat <- matrix(c(1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0), 5, 5)
rs <- rowSums(obs.mat)
cs <- colSums(obs.mat)
web <- rs %*% t(cs)
web <- web/sum(web)
n = sum(obs.mat)
mgen(web, n, keep.species=FALSE, rep.cell=FALSE) # Allowing zero marginal sums
mgen(web, n, keep.species=TRUE, rep.cell=FALSE) # Not allowing zero marginal sums

mgen(mosquin1967, keep.species=FALSE, rep.cell=TRUE) # Allowing zero marginal sums
mgen(mosquin1967, keep.species=TRUE, rep.cell=TRUE) # Not allowing zero marginal sums


# ND, BC, CC
ND(vazquenc)
BC(vazquenc)
BC(vazquenc, rescale=FALSE, weighted=FALSE)
CC(vazquenc)

# nest.smdm
nest.smdm(Safariland)
nest.smdm(Safariland, weighted=TRUE)
nest.smdm(Safariland, weighted=TRUE, decreasing="abund")
nest.smdm(Safariland, weighted=T, decreasing="abund", sort=F)
# identify modules using computeModules:
mod <- computeModules(Safariland)
const <- module2constraints(mod)
nest.smdm(Safariland, constraint=const)
nest.smdm(Safariland, constraint=const, weighted=T)



# nested
nested(vazquenc, method="ALL")

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
networklevel(Safariland, index="NODF")

# nodespec
nodespec(Safariland)
nodespec(Safariland, inf.replace=Inf)

# npartite ...

# nullmodel
nullmodel(motten1982, 2, "r2d")
for (i in 1:6) nullmodel(motten1982, 2, i) # i=6 should return an error message!
nullmodel(motten1982, 2, "mgen")

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
second.extinct(Safariland, participant="both", method="external", ext.row=1:9, ext.col=27:1) # should break!

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

# vaznull
set.seed(1)
m <- matrix(rpois(4, 2), 2, 2)
vaznull(2, m) # Error when m is full (i.e. no 0); works fine with a single 0 already (set seed to 4)  

#vaznullexternal
abun.lower <- c(15,5,2,7,4,8,6,0.01,6)
set.seed(2)
(abun.higher <- rpois(27, lambda=4))
abun.higher[1] <- 0.001
nulls <- vaznullexternal(2, Safariland, abun.higher=abun.higher, abun.lower=abun.lower)
cor(colSums(nulls[[1]]), abun.higher) # close to 1

nulls <- vaznullexternal(2, Safariland)
cor(colSums(nulls[[1]]), colSums(Safariland)) # very close to 1


# V.ratio

# visweb

# web2edges

# webs2array
data(Safariland, vazquenc, vazquec)
allin1 <- webs2array()
allin1 <- webs2array(Safariland)
allin1 <- webs2array(Safariland, vazquenc, vazquec)
str(allin1)
## now we can compute distance between two webs:
vegdist(t(cbind(as.vector(allin1[,,2]), as.vector(allin1[,,3]))), method="jacc")
webinput <- substitute(list(Safariland, vazquenc, vazquec))
as.character(webinput)


# wine



## ----------some not so bright ideas ---------##
###source function to make sure we use the most recent ones:
# NOT SO BRIGHT! because the namespace will import specified functions. When we now load the functions via source, these imported functions are NOT available, unless we explicitly load the package, too:
#setwd("/Users/cdormann/Data/aktuell/BESS/R_FoodWeb/")
#rfiles <- list.files("bipartite/R")
## source files, excluding zzz.r
#for (i in seq_along(rfiles)[-length(rfiles)]) source(paste0("bipartite/R/", rfiles[i]))
# C.score(Safariland) # error: could not find "designdist"

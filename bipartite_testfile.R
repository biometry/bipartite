# This now is the command line version; in particular the R CMD CHECK is better here than in the devtools version (which throws an error related to roxygen)

# workflow: 
## At the bottom of the testfile is the workflow on the old computer, where the compression did not work and thus several detours had to be taken.

Sys.setenv("R_CHECK_RD_VALIDATE_RD2HTML"=FALSE) # to switch off html-syntax check in R CMD check
## DOES NOT WORK! (I used "_" before and after the option, as indicate in the NEWS-section.)
#Sys.getenv()

R CMD build bipartite --compact-vignettes=gs+qpdf
R CMD check bipartite_2.23.tar.gz --as-cran
R CMD install bipartite_2.23.tar.gz
# now check in the testfile below anything that changed whether it actually works!
# upload to https://win-builder.r-project.org/upload.aspx and check on R-devel!


#### check reverse dependencies ####
# identify, which packages depend on bipartite:
devtools::revdep("bipartite")
# [1] "bioregion"     "bipartiteD3"   "bootPLS"       "cassandRa"    
# [5] "econetwork"    "econullnetr"   "HiveR"         "leiden"       
# [9] "NIMAA"         "nos"           "PhageCocktail" "plsRbeta"     
# [13] "plsRglm"       "primer"        "RPANDA"        "tapnet"  

# check whether changes break reverse dependencies:
pak::pkg_install("r-lib/revdepcheck")
library(revdepcheck)
setwd("bipartite/") # move into package-directory for this test
getwd()
revdep_check(num_workers = 4)

# ── CHECK ───────────────────────────────────────────────────────────────────────────────────── 16 packages 
# ✔ bipartiteD3 0.3.2                     ── E: 1     | W: 0     | N: 0                                          
# ✔ cassandRa 0.2.0                       ── E: 0     | W: 0     | N: 0                                          
# ✔ econetwork 0.7.0                      ── E: 1     | W: 0     | N: 0                                          
# ✔ bootPLS 1.0.1                         ── E: 1     | W: 0     | N: 0                                          
# ✖ econullnetr 0.2.1                     ── E: 0  +2 | W: 0     | N: 0                                          
# ✔ bioregion 1.2.0                       ── E: 0     | W: 0     | N: 0                                          
# ✔ HiveR 0.4.0                           ── E: 0     | W: 0     | N: 0                                          
# ✔ leiden 0.4.3.1                        ── E: 1     | W: 0     | N: 0                                          
# ✔ nos 2.0.0                             ── E: 0     | W: 0     | N: 0                                          
# ✔ PhageCocktail 1.0.3                   ── E: 0     | W: 0     | N: 1                                          
# ✔ NIMAA 0.2.1                           ── E: 0     | W: 0     | N: 0                                        
# ✔ plsRbeta 0.3.1                        ── E: 0     | W: 0     | N: 0                                          
# ✔ plsRglm 1.5.1                         ── E: 0     | W: 0     | N: 0                                          
# ✔ primer 1.2.0                          ── E: 0     | W: 0     | N: 2                                          
# ✔ tapnet 0.3                            ── E: 0     | W: 0     | N: 1                                          
# ✖ RPANDA 2.4                            ── E: 0     | W: 0  +1 | N: 0  +1                                      OK: 14  


# econullnetr
library(econullnetr)
sil.null <- generate_null_net(Silene[, 2:7], Silene.plants[, 2:6], sims = 2, c.samples = Silene[, 1], r.samples = Silene.plants[, 1])
plot_bipartite(sil.null)
# change in econullnetr::plot_bipartite: method -> sorting

# RPANDA
plot_network: possible error in plotweb(as.matrix(link), col.low =
                                          col.P, col.high = col.H, col.interaction = "lightgray", empty =
                                          TRUE): unused arguments (col.low = col.P, col.high = col.H,
                                                                   col.interaction = "lightgray")
# offensive code:
try(plotweb(as.matrix(link),col.low=col.P,col.high = col.H,col.interaction = "lightgray",empty=TRUE))
col.low --> lower_color
col.high --> higher_color
col.interaction --> link_color

#### bipartite test file ####

## CHECK R-CONSOLE FOR RED TEXT INDICATING ERRORS!! ##


## run this file after every change in bipartite before submitting it to CRAN!!
library(bipartite)
library(testthat) # if we want to formally test expressions and functions!

#source("/Users/Carsten/Data/aktuell/Networks/bipartite/bipartite/R/computeModules.R")
#source("/Users/Carsten/Data/aktuell/Networks/bipartite/bipartite/R/restrictednull.R")

# lazy load data does not require data to be loaded via "data(.)"!


#  as.one.mode
image(as.one.mode(Safariland))
visweb(as.one.mode(Safariland, project="lower", fill=NA), NA.col="green") #NA.col should have no effect and cause no problem!
# check fill=NA and visweb's NA.col:
visweb(as.one.mode(vazquenc, fill=NA), NA.col="green")

as.one.mode(Safariland, project="lower", legacy=F)
as.one.mode(Safariland, project="lower", legacy=T)

as.one.mode(Safariland, project="lower", weighted=F, legacy=T)
as.one.mode(Safariland, project="lower", legacy=T)

gplot(as.one.mode(Safariland, project="higher", weighted=T, legacy=F))


# array2linkmx
array2linkmx(webs2array(Safariland, vazquenc))

# betalinkr 
#--case1: testdata of frame2webs (a case with low overlap and no shared links) --
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
  lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
  webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
testarray <- frame2webs(testdata, type.out="array")
betalinkr(testarray, distofempty="na", partitioning="poisot")  # returns NA for OS and ST
# betalink(prepare_networks(list(testarray[,,1]))[[1]],  prepare_networks(list(testarray[,,2]))[[1]]) # betalink also returns NaN for this case (but requires package betalink)
betalinkr(testarray, partitioning="poisot") # fixed, now OS is zero
betalinkr(testarray, partitioning="commondenom", partition.st=TRUE) # Novotny-style: plant differences rule here
#--case2: Vazquez-data --
testarray <- webs2array(Safariland, vazarr)  
# various options
# some binary examples
betalinkr(testarray)
betalinkr(testarray, partitioning="poisot", binary=TRUE, index="bray")
betalinkr(testarray, index="jaccard")
betalinkr(testarray, partitioning="poisot", binary=TRUE, index="jaccard")
# some quantitative examples
betalinkr(testarray, binary=FALSE, index="sorensen")
betalinkr(testarray, binary=FALSE, partitioning="poisot")
betalinkr(testarray, partitioning="commondenom", binary=FALSE, proportions=FALSE)
betalinkr(testarray, partitioning="poisot", binary=FALSE, proportions=FALSE)
betalinkr(testarray, index="horn", binary=F, partitioning="poisot")
# partition.st:
betalinkr(testarray, partitioning="commondenom", binary=TRUE, index="sorensen",partition.st=TRUE)
betalinkr(testarray, partitioning="commondenom", binary=FALSE, index="sorensen",partition.st=TRUE)
#--case3: two fully connected  and completely shared webs --
testarray <- array(1:24, dim=c(2, 3, 2))
betalinkr(testarray)
# mostly shared webs
testarray <- array(1:24, dim=c(2, 3, 4))
  set.seed(23) # creates a special case where OS is smaller for "poisot" than for "commondenom"
  testarray[sample(1:24, 10)] <- 0 # setting some entries of above matrix to zero
  testarray <- testarray[, , sample(1:4, 2)] # selecting two sites at random (for now, just developing the function for a 2row-matrix)
betalinkr(testarray, partitioning="poisot", binary=F)
betalinkr(testarray, partitioning="commondenom", binary=F, partition.rr=T) # warning, but shows that OS_poisot < OS_commondenom can be explained by size difference of sharedsp subwebs (OS.rich) being removed by standardization to proportions (only in poisot this is done for the subweb)
#--case4: example from Figure 1 in Poisot et al. 2012 --
metaweb <- matrix(rep(0,25),nrow=5)
metaweb[as.matrix(data.frame(c(2,4,5,5),c(1,2,2,3)))] <- 1  # the adj. matrix for the metaweb
dimnames(metaweb) <- list(letters[1:5], letters[1:5])  # creating species names
web1 <- metaweb
web1[5,2] <- 0 # one link removed from metaweb
web2 <- web1[-1,-1] # top predator removed from web1
web3 <- metaweb[-1,-1] # top predator removed from metaweb
betalinkr_multi(webs2array(web1,web2,web3),partitioning="poisot")       # OS>ST (web1 vs web3)
betalinkr_multi(webs2array(web1,web2,web3), partitioning="commondenom") # OS=ST (web1 vs web3)
# partition.st:
betalinkr_multi(webs2array(web1,web2,web3), partition.st=T) 
# partition.rr:
betalinkr_multi(webs2array(web1,web2,web3), partition.rr=T) # WN.rich=0, but OS.rich>0 (replacement leaves sharedspweb)
# both secondary partitions in one call
betalinkr_multi(webs2array(web1,web2,web3), partition.st=T, partition.rr=T) 


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
system.time(C.score(memmott1999)) # 3s

# CC (closeness centrality)
CC(Safariland)
CC(Safariland, rescale=F)
wf <- as.one.mode(Safariland, project="lower", weighted=F)
wt <- as.one.mode(Safariland, project="lower", weighted=T)
closeness(wf, cmode="suminvundir")
closeness(wt, cmode="suminvundir") # makes no difference!
specieslevel(Safariland, index="closeness") # same as CC(. , rescale=T); no option for rescale=F 


# compart
compart(Safariland) # a comparted network
compart(bezerra2009) # an uncomparted network

# computeModules
## a lot to test here! let's start with the problem of calling computeModules twice in a row:
comp1 <- computeModules(vazquenc)
comp1 <- LPA_wb_plus(vazquenc)$modularity
comp1plus <- LPA_wb_plus(cbind(vazquenc, rep(0, nrow(vazquenc))))$modularity
comp1plus2 <- LPA_wb_plus(rbind(cbind(vazquenc, rep(0, nrow(vazquenc))), rep(0, ncol(vazquenc)+1)) )$modularity
comp2 <- computeModules(vazquenc, forceLPA=TRUE)
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

a<-matrix(c(1,1,1,1,1,1,
            1,1,1,1,1,1,
            1,1,1,1,1,1,
            0,1,0,0,0,0,
            0,1,0,0,0,0,
            0,1,0,0,0,0,
            0,1,0,0,0,0,
            0,1,0,0,0,0,
            1,1,1,1,1,1,
            1,1,0,0,0,0,
            1,1,0,0,0,0), ncol=11, nrow =6)
computeModules(a)
computeModules(a, empty.web=T) # should cause an error!
slot(computeModules(a), "likelihood")


# czvalues
comp1 <- computeModules(memmott1999)
czvalues(comp1, weighted=T)
czvalues(comp1, level="lower")
# plotModuleWeb(comp1)
czvalues(comp3)


# decimalr2dtable
nulls <- decimalr2dtable(100, Safariland)
g.dec <- sapply(nulls, networklevel, index="generality")
plot(density(g.dec[1,]), xlim=c(1, 3))


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
replicate(10, discrepancy(vazquenc[, sample(1:24)])) # varies with sequence!
replicate(10, vegan::nesteddisc(vazquenc[, sample(1:24)])$statistic) # pre-sorts the matrix


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

# extinction: see also second.extinct


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

# nestedness # DEPRECATED!!
#nestedness(Safariland, n.nulls=20)[c(4, 9:20)]
#nestedness(Safariland, null.models=FALSE)$temperature

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
networklevel(Safariland, index="modularity")
replicate(10, networklevel(Safariland[sample(1:9), sample(1:27)], index="discrepancy")) # despite reordering ties, this is not always achieved within the 200 iterations of nesteddisc
networklevel(Safariland, effective=T) # check changes relative to default for "interaction evenness" and "H2"

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
null.t.test(Safariland, index=c("generality", "vulnerability", "connectance","links per species","cluster coefficient"), nrep=4, N=20)
null.t.test(Safariland, index=c("cluster coefficient"), nrep=4, N=20)

# nullmodel
lapply(1:5, function(x) nullmodel(vazquenc, N=2, method=x))
lapply(c("r2d", "swap.web", "vaznull", "shuffle.web", "mgen"), function(x) nullmodel(vazquenc, N=2, method=x))

# PAC
PAC(vazquenc)
PAC(vazquenc>0)

# PDI
PDI(vazquenc, normalise=FALSE, log=TRUE)
PDI(vazquenc>0)

# plotmatrix
S <- sortmatrix(Safariland, topology = "nested", sort_by = "weights")
S <- sortmatrix(kato1990, topology = "nested", sort_by = "weights")
plotmatrix(S$matrix, binary=FALSE)
plotmatrix(S, binary=TRUE)

# plotModuleWeb
if (!exists("comp1")) comp1 <- computeModules(vazquenc)
plotModuleWeb(comp1)
plotModuleWeb(comp1, plotModules = TRUE, rank = TRUE, weighted = TRUE, displayAlabels = TRUE, displayBlabels = TRUE, labsize = 0.6, xlabel = "", ylabel = "", square.border = "lightgreen", fromDepth = 0, upToDepth = -1)

# plotPAC
plotPAC(PAC(vazquenc))
plotPAC(PAC(vazquenc), scaling = 2, plot.scale = 1.5, fill.col = rgb(0.2, 0.3, 0.4, 0.5), arrow.col = rgb(0.4, 0.3, 0.2, 0.5), outby = 0.5, label = F, text=TRUE, circles = T, radius = 1.5)
plotPAC(kevan1970, arrow.col=rainbow(30), text=F) # test multiple colours
plotPAC(kevan1970, arrow.col=rainbow(30), text=F, outby=.9) # test multiple colours

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


# restrictednull
Mod <- computeModules(Safariland)
Part <- module2constraints(Mod)
row.Part <- Part[1:nrow(Safariland)]
col.Part <- Part[(nrow(Safariland)+1):(nrow(Safariland) + ncol(Safariland))]
nulls <- restrictednull(web = Safariland, R.partitions = row.Part, C.partitions = col.Part)
nulls <- restrictednull(web = Safariland, Prior.Pij = "equiprobable", R.partitions = row.Part, C.partitions = col.Part)
nulls <- restrictednull(web = Safariland, conditional.level="matrix") # this should essential be vaznull, I think
nulls <- restrictednull(web = Safariland, Prior.Pij="degreeprob.byarea", conditional.level="areas", R.partitions = row.Part, C.partitions = col.Part)

# robustness

# second.extinct
bs <- second.extinct(Safariland, method="random", participant="both", details=T) 
slope.bipartite(bs) # should return an error with an explanation
bs <- second.extinct(Safariland, method="random", participant="both", details=F) 
slope.bipartite(bs) # should work
second.extinct(Safariland, participant="both", method="external", ext.row=1:9, ext.col=27:1) # should break!

web <- matrix(c(3, 2, 3, 0, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
second.extinct(web, participant = "lower", method = "abundance")
second.extinct(web, participant = "lower", method = "random")
second.extinct(web, participant = "higher", method = "abundance")


# slope.bipartite

# sortmatrix
sortmatrix(Safariland, topology="nested")
sortmatrix (Safariland, topology = "nested", sort_by = "weights")

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
rowSums(vaznull(1, Safariland)[[1]])
unname(rowSums(Safariland))
# often rather large discrepancies in the number of interactions!


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
web2edges(Safariland)
web2edges(Safariland, both.directions=TRUE)
web2edges(as.one.mode(Safariland, project="lower"), is.one.mode=T)

# webs2array
data(Safariland, vazquenc, vazquec)
allin1 <- webs2array()  # returns uninformative error
allin1 <- webs2array(Safariland) # returns informative error
allin1 <- webs2array(Safariland, vazquenc, vazquec) # just works
str(allin1)

otto <- list(Safariland, vazarr)
str(webs2array(otto, vazquec))
testfun <- function(x) webs2array(x)
str(testfun(otto)) # caused an error in the pre-2.18 version, which used only the ellipsis (...), not x.



## now we can compute distance between two webs:
vegdist(t(cbind(as.vector(allin1[,,2]), as.vector(allin1[,,3]))), method="jacc")
webinput <- substitute(list(Safariland, vazquenc, vazquec))
as.character(webinput)


# wine



#### OLD STUFF ####
## Below the testfile is the workflow on the old computer, where the compression did not work and thus several detours had to be taken.

# 0. sync vignette-folder for save-keeping with same directory one level higher (the one-level higher is the one I use for everything and only copy-paste the .Rnw into vignette before building the package).
#
# 1. build package normally (trying to compact the vignette)
R CMD build bipartite --compact-vignettes=gs+qpdf
# 2. go to some webpage and compact there the vignette from inst/doc in the built package, e.g. https://www.ilovepdf.com !!! Only if larger than 1MB !!!
...
# 3. put the thus compacted PDF into inst/doc of the development folder (not into the built!), put in there also the .R and .Rnw files from inst/doc of the just-built .tar.gz!
...
# 4. copy-paste the "build" folder of the .tar.gz of step 1 into the development folder; it contains the data and vignette reference .rdb and .rds needed for a working package!
...
# 5. Delete, in vignettes, the cache-folder and auxiliary LaTeX-files and Sweave.sty; keep only the .Rnw and the .bib, and the figures-folder (and the styles: mee.bst and Sweavel.sty) 
...
# 6. build package anew without rebuilding the vignettes (and without cleaning the docs, in case you are using devtools::build!)
R CMD build bipartite --no-build-vignettes --resave-data
# 7. check all is fine: first locally, then on win-builder (https://win-builder.r-project.org/upload.aspx)
R CMD check bipartite_2.17.tar.gz --as-cran
R CMD install bipartite_2.17.tar.gz # optional; check html of help and link to vignette in RStudio 
# if you get this error: Error in fetch(key) : lazy-load database '/Users/Carsten/Library/R/4.0/library/bipartite/help/bipartite.rdb' is corrupt
# re-start R (RStudio); this is just a point of the install not updating the central help pages (https://stackoverflow.com/questions/30424608/error-in-fetchkey-lazy-load-database).


## Comments on the workflow above:
## ad 0.: The vignette is a pain in the neck! The vignettes folder in the level of the top bipartite folder (in PDFetc) is the one to use for writing and processing the vignette! There are some problems, for example that in the betweenness comparison I call packages not listed in "Depends" of bipartite. Since I don't want to make bipartite dependent on packages that do not work well, I now compile the vignette in this folder with "eval=T", then copy the output for this since R-chunk from the .tex-file into the .Rnw file and set "eval=F" (there is a note to that effect in the .Rnw). Then I put the .Rnw into the vignettes-folder and the inst/doc of the package. What a mess!
## ad 1.: Somehow --compact-vignettes... does not compact at all. I tried all options (both, gs+qpdf, qpdf, gs; always without quotes!), nothing happened. I ran qpdf::pdf_compact and that did work, so qpdf is on my system(s); I have no idea what else to do. 
## ad 2.: You can try 
## tools::compactPDF("/Users/Carsten/Data/aktuell/Networks/bipartite/bipartite/inst/doc", qpdf=Sys.which(Sys.getenv("R_QPDF", "qpdf")), gs_quality = "ebook") ## but for me this did not yield any compression;
## ad 5.: There are other check options, e.g. rhub::check("bipartite_2.16.tar.gz", platform = "fedora-clang-devel") # requires validate_email() before first run; rhub misses some packages or package options (e.g. titlesec and nottoc in tocbibind and hidelinks in hyperref)
##
## Misc:
## * find non-UTF8 characters: find . | egrep [^a-zA-Z0-9_\.\/\-\s]
## * check link to external functions: https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Cross_002dreferences
## * if something goes wrong with the vignette build, fix and start over with the package building (it is difficult to move the right files into the right folders; R even checks the date of .Rnw files!)

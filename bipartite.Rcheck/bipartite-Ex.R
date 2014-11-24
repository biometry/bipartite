pkgname <- "bipartite"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bipartite')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("C.score")
### * C.score

flush(stderr()); flush(stdout())

### Name: C.score
### Title: Calculates the (normalised) mean number of checkerboard
###   combinations (C-score) in a matrix
### Aliases: C.score
### Keywords: package

### ** Examples

m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,TRUE)
C.score(m)
C.score(m, normalise=FALSE)
C.score(m, normalise=FALSE, FUN=print)




cleanEx()
nameEx("H2fun")
### * H2fun

flush(stderr()); flush(stdout())

### Name: H2fun
### Title: Specialisation of a bipartite web.
### Aliases: H2fun
### Keywords: package

### ** Examples

data(Safariland)
H2fun(Safariland)



cleanEx()
nameEx("ND")
### * ND

flush(stderr()); flush(stdout())

### Name: ND
### Title: Normalised degree, betweenness and closeness centrality
### Aliases: ND BC CC
### Keywords: package

### ** Examples

## example:
data(olesen2002flores)
(ndi <- ND(olesen2002flores))
(cci <- CC(olesen2002flores))
(bci <- BC(olesen2002flores))

cor.test(bci[[1]], ndi[[1]], method="spear") # 0.532
cor.test(cci[[1]], ndi[[1]], method="spear") # 0.403

cor.test(bci[[2]], ndi[[2]], method="spear") # 0.738
cor.test(cci[[2]], ndi[[2]], method="spear") # 0.827
## Not run: 
##D ## PLANTS:
##D bc <- bci[[1]]
##D cc <- cci[[1]]
##D nd <- ndi[[1]]
##D # CC:
##D summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) # lower RSE
##D summary(nls(cc ~ c*nd^d, start=list(c=0.072,d=0.2))) 
##D # BC:
##D summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
##D summary(nls(bc ~ c*nd^d, start=list(c=2,d=2))) # lower RSE
##D 
##D ## ANIMALS:
##D bc <- bci[[2]]
##D cc <- cci[[2]]
##D nd <- ndi[[2]]
##D # CC:
##D summary(nls(cc ~ a*nd+b, start=list(a=1,b=1)))  
##D summary(nls(cc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE 
##D # BC:
##D summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
##D summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
## End(Not run)



cleanEx()
nameEx("PAC")
### * PAC

flush(stderr()); flush(stdout())

### Name: PAC
### Title: Potential for Apparent Competition
### Aliases: PAC
### Keywords: package

### ** Examples

data(Safariland)
PAC(Safariland)



cleanEx()
nameEx("PDI")
### * PDI

flush(stderr()); flush(stdout())

### Name: PDI
### Title: Paired Differences Index
### Aliases: PDI
### Keywords: package

### ** Examples

data(Safariland)
PDI(Safariland) # for pollinators
PDI(t(Safariland), log=TRUE) # for plants



cleanEx()
nameEx("Safariland")
### * Safariland

flush(stderr()); flush(stdout())

### Name: Safariland
### Title: A pollination web from Argentina
### Aliases: Safariland
### Keywords: datasets

### ** Examples

data(Safariland)
plotweb(Safariland)



cleanEx()
nameEx("V.ratio")
### * V.ratio

flush(stderr()); flush(stdout())

### Name: V.ratio
### Title: Calculates the variance-ratio as suggested by Schluter (1984)
### Aliases: V.ratio
### Keywords: package

### ** Examples

data(Safariland)
V.ratio(Safariland)



cleanEx()
nameEx("as.one.mode")
### * as.one.mode

flush(stderr()); flush(stdout())

### Name: as.one.mode
### Title: Conversion of a network matrix
### Aliases: as.one.mode
### Keywords: package

### ** Examples

data(Safariland)
image(Safariland)
image(as.one.mode(Safariland))



cleanEx()
nameEx("as.tnet")
### * as.tnet

flush(stderr()); flush(stdout())

### Name: as.tnet
### Title: Ensures that networks conform to the tnet stardards
### Aliases: as.tnet

### ** Examples

## Load sample data
sample <- rbind(
c(1,2,4),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(4,2,1),
c(5,2,2),
c(5,6,1),
c(6,5,1))

## Run the programme
as.tnet(sample)




cleanEx()
nameEx("barrett1987")
### * barrett1987

flush(stderr()); flush(stdout())

### Name: barrett1987
### Title: Individuals caught in a pollination web in boreal Canada.
### Aliases: barrett1987
### Keywords: datasets

### ** Examples

data(barrett1987)



cleanEx()
nameEx("betweenness_w")
### * betweenness_w

flush(stderr()); flush(stdout())

### Name: betweenness_w
### Title: Betweenness centrality in a weighted network
### Aliases: betweenness_w

### ** Examples

## Load sample data
sampledata <- rbind(
c(1,2,1),
c(1,3,5),
c(2,1,1),
c(2,4,6),
c(3,1,5),
c(3,4,10),
c(4,2,6),
c(4,3,10))

## Run the programme
betweenness_w(sampledata)




cleanEx()
nameEx("bezerra2009")
### * bezerra2009

flush(stderr()); flush(stdout())

### Name: bezerra2009
### Title: Individuals observed in a flower-visitation network of
###   oil-collecting bees in a Brazilian steppe.
### Aliases: bezerra2009
### Keywords: datasets

### ** Examples

data(bezerra2009)



cleanEx()
nameEx("bipartite-package")
### * bipartite-package

flush(stderr()); flush(stdout())

### Name: bipartite-package
### Title: Analysis of bipartite ecological webs
### Aliases: bipartite-package bipartite
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D plotweb(Safariland)
##D visweb(Safariland)
##D networklevel(Safariland)
##D specieslevel(Safariland)
## End(Not run)



cleanEx()
nameEx("closeness_w")
### * closeness_w

flush(stderr()); flush(stdout())

### Name: closeness_w
### Title: Closeness centrality in a weighted network
### Aliases: closeness_w

### ** Examples

## Load sample data
sampledata <- rbind(
c(1,2,4),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(4,2,1),
c(5,2,2),
c(5,6,1),
c(6,5,1))

## Run the programme
closeness_w(sampledata)




cleanEx()
nameEx("compart")
### * compart

flush(stderr()); flush(stdout())

### Name: compart
### Title: Detects compartments
### Aliases: compart
### Keywords: package

### ** Examples

# make a nicely comparted web:
web <- matrix(0, 10,10)
web[1,1:3] <- 1 
web[2,4:5] <- 1 
web[3:7, 6:8] <- 1
web[8:10, 9:10] <- 1
web <- web[-c(4:5),] #oh, and make it asymmetric!
web <- web[,c(1:5, 9,10, 6:8)] #oh, and make it non-diagonal
compart(web)

# or, standard, use Safariland as example:
data(Safariland)
compart(Safariland)



cleanEx()
nameEx("computeModules")
### * computeModules

flush(stderr()); flush(stdout())

### Name: computeModules
### Title: Function "computeModules"
### Aliases: computeModules cM readModuleData deleteModuleData
### Keywords: Methods and Generic Functions Clustering

### ** Examples

	## Not run: 
##D 		data(small1976)
##D 		(res <- computeModules(small1976)) # takes several minutes!
##D 		plotModuleWeb(res)
##D 	
## End(Not run)



cleanEx()
nameEx("czvalues")
### * czvalues

flush(stderr()); flush(stdout())

### Name: czvalues
### Title: Computes c and z for network modules
### Aliases: czvalues
### Keywords: package

### ** Examples

## Not run: 
##D data(memmott1999)
##D set.seed(2)
##D mod <- computeModules(memmott1999, steps=1E4)
##D cz <- czvalues(mod)
##D plot(cz[[1]], cz[[2]], pch=16, xlab="c", ylab="z", cex=0.8, xlim=c(0,1), las=1)
##D abline(v=0.62) # threshold of Olesen et al. 2007
##D abline(h=2.5)   # dito
##D text(cz[[1]], cz[[2]], names(cz[[1]]), pos=4, cex=0.7)
## End(Not run)



cleanEx()
nameEx("degreedistr")
### * degreedistr

flush(stderr()); flush(stdout())

### Name: degreedistr
### Title: Fits functions to cumulative degree distributions of both
###   trophic levels of a network.
### Aliases: degreedistr
### Keywords: htest

### ** Examples

data(Safariland)
degreedistr(Safariland)



cleanEx()
nameEx("dfun")
### * dfun

flush(stderr()); flush(stdout())

### Name: dfun
### Title: Calculates standardised specialisation index d' (d prime) for
###   each species in the lower trophic level of a bipartite network.
### Aliases: dfun
### Keywords: htest

### ** Examples

data(Safariland)
dfun(Safariland) # gives d-values for the lower trophic level
# now using independent abundance estimates for higher trophic level:
dfun(Safariland, abuns=runif(ncol(Safariland)))

dfun(t(Safariland)) #gives d-values for the higher trophic level



cleanEx()
nameEx("discrepancy")
### * discrepancy

flush(stderr()); flush(stdout())

### Name: discrepancy
### Title: Calculates discrepancy of a matrix
### Aliases: discrepancy
### Keywords: package

### ** Examples

## Not run: 
##D #nulls <- replicate(100, discrepancy(commsimulator(Safariland, 
##D 		method="quasiswap")))
##D nulls <- simulate(vegan::nullmodel(Safariland, method="quasiswap"), nsim = 100)
##D null.res <- apply(nulls, 3, discrepancy)
##D hist(null.res)
##D obs <- discrepancy(Safariland)
##D abline(v=obs, lwd=3, col="grey")
##D c("p value"=min(sum(null.res>obs), sum(null.res<obs))/length(null.res))
##D # calculate Brualdi & Sanderson's Na-value (i.e. the z-score):
##D c("N_a"=(unname(obs)-mean(null.res))/sd(null.res))
## End(Not run)



cleanEx()
nameEx("elberling1999")
### * elberling1999

flush(stderr()); flush(stdout())

### Name: elberling1999
### Title: No. of visits in a pollination web of arctic-alpine Sweden
### Aliases: elberling1999
### Keywords: datasets

### ** Examples

data(barrett1987)
## maybe str(barrett1987) ; plot(barrett1987) ...



cleanEx()
nameEx("empty")
### * empty

flush(stderr()); flush(stdout())

### Name: empty
### Title: Deletes empty rows and columns from a matrix.
### Aliases: empty
### Keywords: package

### ** Examples

data(Safariland)
web <- Safariland
web[,5] <- 0
empty(web, count=TRUE)
attr(empty(web), "empty")




cleanEx()
nameEx("endpoint")
### * endpoint

flush(stderr()); flush(stdout())

### Name: endpoint
### Title: Computes end-point degrees for a bipartite network
### Aliases: endpoint
### Keywords: package

### ** Examples

# reproduces the example of Gitarranz et al. (2011):
data(memmott1999)
ends <- endpoint(memmott1999)
weights.mean <- tapply(memmott1999, ends, mean)
ends.weights <- tapply(ends, ends, mean)
plot(weights.mean, ends.weights, log="xy", pch=16)



cleanEx()
nameEx("extinction")
### * extinction

flush(stderr()); flush(stdout())

### Name: extinction
### Title: Simulates extinction of a species from a bipartite network
### Aliases: extinction
### Keywords: package

### ** Examples

## Not run: 
##D 	data(Safariland)
##D 	(w <- extinction(Safariland, participant="lower", method="abun"))
##D 	empty(w, count=TRUE)
## End(Not run)



cleanEx()
nameEx("fc")
### * fc

flush(stderr()); flush(stdout())

### Name: fc
### Title: Calculates the functional complementarity for the rows of a web
### Aliases: fd fc
### Keywords: package

### ** Examples

data(Safariland)
fc(Safariland)    
fc(t(Safariland), dist="canberra", method="complete")



cleanEx()
nameEx("frame2webs")
### * frame2webs

flush(stderr()); flush(stdout())

### Name: frame2webs
### Title: Converts a table of observations into a network matrix
### Aliases: frame2webs
### Keywords: package

### ** Examples

testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
frame2webs(testdata,type.out="array")
sapply(frame2webs(testdata,type.out="list"), networklevel, index=c("connectance", "C score"))



cleanEx()
nameEx("genweb")
### * genweb

flush(stderr()); flush(stdout())

### Name: genweb
### Title: Generate a random bipartite web
### Aliases: genweb
### Keywords: package

### ** Examples

genweb()



cleanEx()
nameEx("grouplevel")
### * grouplevel

flush(stderr()); flush(stdout())

### Name: grouplevel
### Title: Analysis of bipartite webs at the level of each of the two
###   levels (groups) of the network
### Aliases: grouplevel one.grouplevel
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D grouplevel(Safariland)
##D grouplevel(Safariland, level="lower", weighted=FALSE) #excludes degree distribution fits
## End(Not run)



cleanEx()
nameEx("inouye1988")
### * inouye1988

flush(stderr()); flush(stdout())

### Name: inouye1988
### Title: A pollination network from the Snowy Mountains of New South
###   Wales, Australia
### Aliases: inouye1988
### Keywords: datasets

### ** Examples

data(inouye1988)
plotweb(inouye1988)



cleanEx()
nameEx("junker2013")
### * junker2013

flush(stderr()); flush(stdout())

### Name: junker2013
### Title: Flower visitation network
### Aliases: junker2013
### Keywords: datasets

### ** Examples

data(junker2013)
## Not run: plotweb(junker2013)



cleanEx()
nameEx("kato1990")
### * kato1990

flush(stderr()); flush(stdout())

### Name: kato1990
### Title: No. of individuals caught in a pollination web of a Japanese
###   beech forest
### Aliases: kato1990
### Keywords: datasets

### ** Examples

data(kato1990)
## maybe str(kato1990) ; plot(kato1990) ...



cleanEx()
nameEx("kevan1970")
### * kevan1970

flush(stderr()); flush(stdout())

### Name: kevan1970
### Title: A pollination network from Northern Ellesmere Island, Canada
### Aliases: kevan1970
### Keywords: datasets

### ** Examples

data(kevan1970)



cleanEx()
nameEx("linklevel")
### * linklevel

flush(stderr()); flush(stdout())

### Name: linklevel
### Title: Indices of a bipartite network at the link level
### Aliases: linklevel
### Keywords: package

### ** Examples

data(Safariland)
linklevel(Safariland)



cleanEx()
nameEx("listModuleInformation")
### * listModuleInformation

flush(stderr()); flush(stdout())

### Name: listModuleInformation
### Title: Function "listModuleInformation"
### Aliases: listModuleInformation
### Keywords: Methods and Generic Functions Clustering

### ** Examples

## Not run: 
##D data(small1976)
##D 
##D moduleWebObject = computeModules(small1976);
##D moduleList = listModuleInformation(moduleWebObject);
## End(Not run)



cleanEx()
nameEx("memmott1999")
### * memmott1999

flush(stderr()); flush(stdout())

### Name: memmott1999
### Title: Flower visitation network from a meadow near Bristol, UK
### Aliases: memmott1999
### Keywords: datasets

### ** Examples

data(memmott1999)
## maybe str(memmott1999) ; plot(memmott1999) ...



cleanEx()
nameEx("mgen")
### * mgen

flush(stderr()); flush(stdout())

### Name: mgen
### Title: Generate simulated network according to a given probability
###   matrix
### Aliases: mgen

### ** Examples

## Not run: 
##D ## Generate simulated matrix from homogeneous probability matrix
##D probmat <- matrix(1/100, 10, 10)
##D mgen(web=probmat, n=100)
##D 
##D ## Generate binary matrix with probability proportional to degree
##D ## of an observed binary matrix m
##D obs.mat <- matrix(c(1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0), 5, 5)
##D rs <- rowSums(obs.mat)
##D cs <- colSums(obs.mat)
##D web <- rs ##D 
##D web <- web/sum(web)
##D n = sum(obs.mat)
##D mgen(web, n, keep.species=FALSE, rep.cell=FALSE) # Allowing zero marginal sums
##D mgen(web, n, keep.species=TRUE, rep.cell=FALSE) # Not allowing zero marginal sums
##D 
##D ## Generate quantitative matrix with probability proportional
##D ## to interaction frequency in an observed matrix m (which is the default of 'autotransform'):
##D mgen(mosquin1967, keep.species=FALSE, rep.cell=TRUE) # Allowing zero marginal sums
##D mgen(mosquin1967, keep.species=TRUE, rep.cell=TRUE) # Not allowing zero marginal sums
## End(Not run)



cleanEx()
nameEx("moduleWeb-class")
### * moduleWeb-class

flush(stderr()); flush(stdout())

### Name: moduleWeb-class
### Title: Class "moduleWeb"
### Aliases: moduleWeb-class
### Keywords: classes modules moduleWeb modularity

### ** Examples

showClass("moduleWeb")



cleanEx()
nameEx("mosquin1967")
### * mosquin1967

flush(stderr()); flush(stdout())

### Name: mosquin1967
### Title: Flower visitation network from Melville Island, Northwest
###   Territories, Canada
### Aliases: mosquin1967
### Keywords: datasets

### ** Examples

data(mosquin1967)
## maybe str(mosquin1967) ; plot(mosquin1967) ...



cleanEx()
nameEx("motten1982")
### * motten1982

flush(stderr()); flush(stdout())

### Name: motten1982
### Title: A spring flower visitation network from North Carolina, USA
### Aliases: motten1982
### Keywords: datasets

### ** Examples

data(motten1982)
## maybe str(motten1982) ; plot(motten1982) ...



cleanEx()
nameEx("nested")
### * nested

flush(stderr()); flush(stdout())

### Name: nested
### Title: Calculates any of several measures of nestedness
### Aliases: nested
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D nested(Safariland, "ALL")
##D nested(Safariland, "ALL", rescale=TRUE)
##D # illustration that non-normalised C.score and checker are the same:
##D nested(Safariland, c("C.score", "checker"), normalise=FALSE)
## End(Not run)



cleanEx()
nameEx("nestedcontribution")
### * nestedcontribution

flush(stderr()); flush(stdout())

### Name: nestedcontribution
### Title: Calculates the per-species contribution to nestedness (z-score
###   relative to null model)
### Aliases: nestedcontribution
### Keywords: package

### ** Examples

  data(Safariland)
  ## Not run: 
##D     nestedcontribution(Safariland)
##D   
## End(Not run)



cleanEx()
nameEx("nestedness")
### * nestedness

flush(stderr()); flush(stdout())

### Name: nestedness
### Title: Calculates nestedness temperature of presence/absence matrices
### Aliases: nestedness
### Keywords: package

### ** Examples

	## Not run: 
##D 	data(vazarr)
##D 	nestedness(vazarr) # null models are calculated
##D 	# no null models, much faster for bigger matrices:
##D 	nestedness(vazarr, null.models=FALSE) 
##D 	nestedness(vazarr, n.nulls=300, n.gen=300)
##D 	
## End(Not run)



cleanEx()
nameEx("nestedrank")
### * nestedrank

flush(stderr()); flush(stdout())

### Name: nestedrank
### Title: Calculates the rank of a species in a matrix sorted for maximum
###   nestedness
### Aliases: nestedrank
### Keywords: package

### ** Examples

## Not run: 
##D ranks <- sapply(c("nodf", "binmatnest", "wine", "sort"), function(x) 
##D   nestedrank(Safariland, method=x)[[2]])
##D cor(ranks) # high correlation between sort and other indicate that only abundance matters
## End(Not run)



cleanEx()
nameEx("networklevel")
### * networklevel

flush(stderr()); flush(stdout())

### Name: networklevel
### Title: Analysis of bipartite webs at the level of the entire network
### Aliases: networklevel .networklevel
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D networklevel(Safariland)
##D networklevel(Safariland, index="ALLBUTDD") #excludes degree distribution fits
## End(Not run)



cleanEx()
nameEx("nodespec")
### * nodespec

flush(stderr()); flush(stdout())

### Name: nodespec
### Title: Calculates the node-based specialisation index
### Aliases: nodespec
### Keywords: package

### ** Examples

data(Safariland)
nodespec(Safariland, inf.replace=Inf)



cleanEx()
nameEx("npartite")
### * npartite

flush(stderr()); flush(stdout())

### Name: npartite
### Title: Computes indices for a masked-one-mode network
### Aliases: npartite
### Keywords: package

### ** Examples

image(aomw <- as.one.mode(Safariland, fill=NA))
npartite(aomw)
networklevel(Safariland, index=c("connectance", "links per species"))



cleanEx()
nameEx("null.distr")
### * null.distr

flush(stderr()); flush(stdout())

### Name: null.distr
### Title: Null model based on fitted marginal distribution
### Aliases: null.distr
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D null.distr(N=2, Safariland)
##D null.distr(N=2, Safariland, distr="negbin")
##D 
##D round(networklevel(Safariland, "info"), 3)
##D sapply(null.distr(N=5, Safariland), function(x) networklevel(x, index="info"))
##D # highly connected
##D sapply(null.distr(N=5, Safariland, distr="negbin"), function(x) networklevel(x, 
##D 	index="info")[3])
##D # similarly highly connected
## End(Not run)



cleanEx()
nameEx("null.t.test")
### * null.t.test

flush(stderr()); flush(stdout())

### Name: null.t.test
### Title: Compares observed pattern to random webs.
### Aliases: null.t.test
### Keywords: package

### ** Examples

data(mosquin1967)
null.t.test(mosquin1967, index=c("generality", "vulnerability",
    "cluster coefficient", "H2", "ISA", "SA"), nrep=2, N=10)



cleanEx()
nameEx("nullmodel")
### * nullmodel

flush(stderr()); flush(stdout())

### Name: nullmodel
### Title: Generates null models for network analysis
### Aliases: nullmodel
### Keywords: package

### ** Examples

## Not run: 
##D 	data(Safariland)
##D 	nullmodel(Safariland, N=2, method=1)
##D 	nullmodel(Safariland>0, N=2, method=4)
##D 	# analysis example:
##D 	obs <- unlist(networklevel(Safariland, index="weighted nestedness"))
##D 	nulls <- nullmodel(Safariland, N=100, method=1)
##D 	null <- unlist(sapply(nulls, networklevel, index="weighted nestedness")) #takes a while ...
##D 	
##D 	plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))), 
##D 		main="comparison of observed with null model Patefield")
##D 	abline(v=obs, col="red", lwd=2)    
##D 	
##D 	praw <- sum(null>obs) / length(null)
##D 	ifelse(praw > 0.5, 1-praw, praw)    # P-value
##D 	
##D 	# comparison of null model 4 and 5 for binary:
##D 	nulls4 <- nullmodel(Safariland>0, N=100, method=4)
##D 	nulls5 <- nullmodel(Safariland>0, N=100, method=5)
##D 	null4 <- unlist(sapply(nulls4, networklevel, index="weighted nestedness"))
##D 	null5 <- unlist(sapply(nulls5, networklevel, index="weighted nestedness"))
##D 	
##D 	
##D 	plot(density(null4), xlim=range(c(null4, null5)), lwd=2, 
##D 		main="comparison of null models")
##D 	lines(density(null5), col="red", lwd=2)
##D 	legend("topright", c("shuffle", "mgen"), col=c("black", "red"), lwd=c(2,2), 
##D 		bty="n", cex=1.5)
##D 	abline(v=networklevel(Safariland>0, index="weighted nestedness"))
##D 	
## End(Not run)



cleanEx()
nameEx("olesen2002aigrettes")
### * olesen2002aigrettes

flush(stderr()); flush(stdout())

### Name: olesen2002aigrettes
### Title: A flower visitation network from the Azores
### Aliases: olesen2002aigrettes
### Keywords: datasets

### ** Examples

data(olesen2002aigrettes)
## maybe str(olesen2002aigrettes) ; plot(olesen2002aigrettes) ...



cleanEx()
nameEx("olesen2002flores")
### * olesen2002flores

flush(stderr()); flush(stdout())

### Name: olesen2002flores
### Title: Another flower visitation network from the Azores
### Aliases: olesen2002flores
### Keywords: datasets

### ** Examples

data(olesen2002flores)
## maybe str(olesen2002flores) ; plot(olesen2002flores) ...



cleanEx()
nameEx("ollerton2003")
### * ollerton2003

flush(stderr()); flush(stdout())

### Name: ollerton2003
### Title: ollerton2003
### Aliases: ollerton2003
### Keywords: datasets

### ** Examples

data(ollerton2003)
plotweb(ollerton2003)
## maybe str(ollerton2003) ; plot(ollerton2003) ...



cleanEx()
nameEx("plotModuleWeb")
### * plotModuleWeb

flush(stderr()); flush(stdout())

### Name: plotModuleWeb
### Title: Function "plotModuleWeb"
### Aliases: plotModuleWeb prepareWebForPlottingModules addEmptyRowToMatrix
###   addEmptyColToMatrix getModuleCoordinates isCorrectModuleWebObject
###   drawModules
### Keywords: Methods and Generic Functions Clustering

### ** Examples

## Not run: 
##D data(small1976)
##D 
##D moduleWebObject = computeModules(small1976);
##D plotModuleWeb(moduleWebObject);
## End(Not run)



cleanEx()
nameEx("plotPAC")
### * plotPAC

flush(stderr()); flush(stdout())

### Name: plotPAC
### Title: Function to draw a circular plot to visualise potential apparent
###   competition (PAC)
### Aliases: plotPAC
### Keywords: package

### ** Examples

## Not run: 
##D data(kevan1970)
##D plotPAC(kevan1970)
##D 
##D data(Safariland)
##D plotPAC(Safariland, plot.scale=1, fill.col="red", arrow.col="orange", 
##D 	circles=TRUE, radius=1)
## End(Not run)



cleanEx()
nameEx("plotweb")
### * plotweb

flush(stderr()); flush(stdout())

### Name: plotweb
### Title: Visualize a bipartite interaction matrix (e.g. a foodweb)
### Aliases: plotweb
### Keywords: package

### ** Examples

data(Safariland)
plotweb(Safariland)

# shorter labels
plotweb(Safariland, high.lablength=3, low.lablength=0, arrow="down")

# centered triangles for displaying interacions
plotweb(Safariland, text.rot=90, arrow="down.center", col.interaction="wheat2",
	y.lim=c(-1,2.5))

#orginal sequence, up arrows and different box width
plotweb(Safariland, method="normal", arrow="up", y.width.low=0.3, low.lablength=4)
# interactions as lines
plotweb(Safariland, arrow="both", y.width.low=0.05, text.rot=90, col.high="blue", 
	col.low="green")

# add an abundance vector for lower trophic species 
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction="purple", 
	y.width.low=0.05, y.width.high=0.05)

plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction ="red", 
	bor.col.interaction="red", arrow="down")

# now vectors for all colours can be given, to mark certain species or 
# interactions. Colour vectors are recycled if not of appropriate length
plotweb(Safariland,col.high=c("orange","green"))
plotweb(Safariland,col.low=c("orange","green"),col.high=c("white","grey","purple"),
	text.high.col=c("blue","red"), col.interaction=c("red",rep("green",26),rep("brown",242)),
	bor.col.interaction=c(rep("green",26),rep("brown",242)),method="normal", 
	text.rot=90, low.lablength=10, high.lablength=5)


#example one (tritrophic)
plotweb(Safariland,y.width.low=0.1, y.width.high=0.05,method="normal", 
	y.lim=c(0,3), arrow="up", adj.high=c(0.5,1.5), col.high="orange",
	high.lablength=3,high.lab.dis=0)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal",
	add=TRUE,low.y=1.5,high.y=2.5, col.low="green", text.low.col="red", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1),low.plot=FALSE)

#example two (4 trophic with abundance)
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, high.abun=low.abun, col.interaction="purple", 
	y.lim=c(0,4.5), high.lablength=0, arrow="up", method="normal", 
	y.width.high=0.05)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=1.7,high.y=2.7, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4, 
	high.lablength=0)

plotweb(Safariland,y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=2.95, high.y=3.95, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4)

# now some examples with the abuns.type-option:
plotweb(Safariland, abuns.type='independent',arrow="down.center")
plotweb(Safariland, abuns.type='additional',arrow="down.center")




cleanEx()
nameEx("printoutModuleInformation")
### * printoutModuleInformation

flush(stderr()); flush(stdout())

### Name: printoutModuleInformation
### Title: Function printoutModuleInformation
### Aliases: printoutModuleInformation
### Keywords: Methods and Generic Functions Clustering

### ** Examples

## Not run: 
##D data(small1976)
##D moduleWebObject = computeModules(small1976)
##D printoutModuleInformation(moduleWebObject)
## End(Not run)



cleanEx()
nameEx("projecting_tm")
### * projecting_tm

flush(stderr()); flush(stdout())

### Name: projecting_tm
### Title: Projecting binary and weighted two-mode networks onto weighted
###   one-mode networks.
### Aliases: projecting_tm

### ** Examples

#please download and investigate tnet for examples!



cleanEx()
nameEx("r2dexternal")
### * r2dexternal

flush(stderr()); flush(stdout())

### Name: r2dexternal
### Title: Generates null models for network analysis by considering
###   external abundances
### Aliases: r2dexternal
### Keywords: package

### ** Examples

## Not run: 
##D 	abun.lower <- c(15,5,2,7,4,8,6,0.01,6)
##D 	 set.seed(2)
##D 	(abun.higher <- rpois(27, lambda=4))
##D 	abun.higher[1] <- 0.001
##D 	sum(ext.polls)
##D 	## note: external abundances do not sum up; this is intentional!!
##D 	r2dexternal(2, Safariland, abun.higher=abun.higher, abun.lower=abun.lower)
##D 	r2dexternal(2, Safariland, abun.higher=abun.higher)
##D 
##D 	
## End(Not run)



cleanEx()
nameEx("robustness")
### * robustness

flush(stderr()); flush(stdout())

### Name: robustness
### Title: Robustness to species extinctions
### Aliases: robustness
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D ex <- second.extinct(Safariland, participant="lower", method="random", nrep=100, 
##D 	details=FALSE)
##D robustness(ex)
## End(Not run)



cleanEx()
nameEx("schemske1978")
### * schemske1978

flush(stderr()); flush(stdout())

### Name: schemske1978
### Title: A flower visitation network from Urbana, IL, USA
### Aliases: schemske1978
### Keywords: datasets

### ** Examples

data(schemske1978)
## maybe str(schemske1978) ; plot(schemske1978) ...



cleanEx()
nameEx("second.extinct")
### * second.extinct

flush(stderr()); flush(stdout())

### Name: second.extinct
### Title: Secondary extinctions in bipartite networks
### Aliases: second.extinct
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D (ex <- second.extinct(Safariland, participant="lower", method="random", nrep=50, 
##D 	details=TRUE))
##D (ex <- second.extinct(Safariland, participant="lower", method="random", nrep=50, 
##D 	details=FALSE))
## End(Not run)



cleanEx()
nameEx("shuffle.web")
### * shuffle.web

flush(stderr()); flush(stdout())

### Name: shuffle.web
### Title: Shuffle web entries
### Aliases: shuffle.web
### Keywords: package

### ** Examples


data(Safariland)

shuffle.web(Safariland, N=2)




cleanEx()
nameEx("slope.bipartite")
### * slope.bipartite

flush(stderr()); flush(stdout())

### Name: slope.bipartite
### Title: Slope of extinction simulation
### Aliases: slope.bipartite
### Keywords: package

### ** Examples

## Not run: 
##D data(Safariland)
##D ex <- second.extinct(Safariland, participant="lower", method="random", nrep=100,
##D details=FALSE)
##D slope.bipartite(ex)
## End(Not run)



cleanEx()
nameEx("small1976")
### * small1976

flush(stderr()); flush(stdout())

### Name: small1976
### Title: A flower visitation network from a peat bog in Ottawa, Canada
### Aliases: small1976
### Keywords: datasets

### ** Examples

data(small1976)
## maybe str(small1976) ; plot(small1976) ...



cleanEx()
nameEx("sortweb")
### * sortweb

flush(stderr()); flush(stdout())

### Name: sortweb
### Title: Function to sort bipartite webs
### Aliases: sortweb
### Keywords: package

### ** Examples

data(Safariland)
web <- Safariland

sortweb(Safariland, sort.order="dec") 
#rarest species first:
plotweb(sortweb(Safariland, sort.order="inc"), method="normal")   
visweb(sortweb(Safariland,sort.order="inc"), type="diagonal", 
	square="compartment", text="none", frame=TRUE)

# sorted by a given (here random) sequence
sequence <- list(seq.higher=sample(colnames(Safariland)), 
	seq.lower=sample(rownames(Safariland)))
web.ordered <- sortweb(web, sort.order="seq", sequence=sequence)




cleanEx()
nameEx("specieslevel")
### * specieslevel

flush(stderr()); flush(stdout())

### Name: specieslevel
### Title: Calculate various indices for network properties at the species
###   level
### Aliases: specieslevel
### Keywords: package

### ** Examples

data(Safariland)
## Not run: 
##D specieslevel(Safariland)
## End(Not run)
specieslevel(Safariland, index="ALLBUTD")[[2]]



cleanEx()
nameEx("strength")
### * strength

flush(stderr()); flush(stdout())

### Name: strength
### Title: Computes species strength according to either of two definitions
### Aliases: strength
### Keywords: package

### ** Examples

data(Safariland)
s1 <- strength(Safariland, type="Barrat")
s2 <- strength(Safariland, type="Bascompte")
plot(s1, s2, log="x")
cor.test(s1, s2, type="ken")
# for lower level:
strength(t(Safariland))




cleanEx()
nameEx("swap.web")
### * swap.web

flush(stderr()); flush(stdout())

### Name: swap.web
### Title: Creates null model for bipartite networks
### Aliases: swap.web
### Keywords: package

### ** Examples


swap.web(Safariland, N=2)




cleanEx()
nameEx("symmetrise_w")
### * symmetrise_w

flush(stderr()); flush(stdout())

### Name: symmetrise_w
### Title: Symmetrise_w
### Aliases: symmetrise_w

### ** Examples

## Load sample data
sample <- rbind(
c(1,2,2),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(5,2,2),
c(5,6,1))

## Run the programme
symmetrise_w(sample, method="MAX")



cleanEx()
nameEx("tnet_igraph")
### * tnet_igraph

flush(stderr()); flush(stdout())

### Name: tnet_igraph
### Title: Exports a tnet network to an igraph object
### Aliases: tnet_igraph

### ** Examples

## Load sample data
sample <- rbind(
c(1,2,4),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(4,2,1),
c(5,2,2),
c(5,6,1),
c(6,5,1))

## Run the programme
tnet_igraph(sample, type="weighted one-mode tnet")




cleanEx()
nameEx("togetherness")
### * togetherness

flush(stderr()); flush(stdout())

### Name: togetherness
### Title: Calculates the number of identical co-presences and co-absences
###   for species-on-islands
### Aliases: togetherness
### Keywords: package

### ** Examples

(m <- matrix(c(0,1,0,0,1,1,0,1,1,0), ncol=2, byrow=TRUE))
togetherness(m)
# or, with two togethernesses:
(n <- matrix(c(0,1,1,0,1,1,0,0,1,1,0,1,0,1), ncol=2, byrow=TRUE))
togetherness(n, normalise=FALSE)

data(Safariland)
togetherness(m)



cleanEx()
nameEx("vazarr")
### * vazarr

flush(stderr()); flush(stdout())

### Name: vazarr
### Title: A pollination network.
### Aliases: vazarr
### Keywords: datasets

### ** Examples

data(vazarr)
## maybe str(vazarr) ; plot(vazarr) ...



cleanEx()
nameEx("vazcer")
### * vazcer

flush(stderr()); flush(stdout())

### Name: vazcer
### Title: A pollination network.
### Aliases: vazcer
### Keywords: datasets

### ** Examples

data(vazcer)
## maybe str(vazcer) ; plot(vazcer) ...



cleanEx()
nameEx("vazllao")
### * vazllao

flush(stderr()); flush(stdout())

### Name: vazllao
### Title: A pollination network.
### Aliases: vazllao
### Keywords: datasets

### ** Examples

data(vazllao)
## maybe str(vazllao) ; plot(vazllao) ...



cleanEx()
nameEx("vazmasc")
### * vazmasc

flush(stderr()); flush(stdout())

### Name: vazmasc
### Title: A pollination network.
### Aliases: vazmasc
### Keywords: datasets

### ** Examples

data(vazmasc)
## maybe str(vazmasc) ; plot(vazmasc) ...



cleanEx()
nameEx("vazmasnc")
### * vazmasnc

flush(stderr()); flush(stdout())

### Name: vazmasnc
### Title: A pollination network.
### Aliases: vazmasnc
### Keywords: datasets

### ** Examples

data(vazmasnc)
## maybe str(vazmasnc) ; plot(vazmasnc) ...



cleanEx()
nameEx("vaznull")
### * vaznull

flush(stderr()); flush(stdout())

### Name: vaznull
### Title: Null model with constrained connectance and moderately
###   constrained marginal totals
### Aliases: vaznull
### Keywords: package

### ** Examples

	## Not run: 
##D 		data(Safariland)
##D 		networklevel(Safariland, index="info")
##D 		networklevel(vaznull(1, Safariland)[[1]], index="info")
##D 		system.time(vaznull(10, Safariland))
##D 		system.time(swap.web(10, Safariland))
##D 	
## End(Not run)



cleanEx()
nameEx("vazquec")
### * vazquec

flush(stderr()); flush(stdout())

### Name: vazquec
### Title: A pollination network.
### Aliases: vazquec
### Keywords: datasets

### ** Examples

data(vazquec)
## maybe str(vazquec) ; plot(vazquec) ...



cleanEx()
nameEx("vazquenc")
### * vazquenc

flush(stderr()); flush(stdout())

### Name: vazquenc
### Title: A pollination network.
### Aliases: vazquenc
### Keywords: datasets

### ** Examples

data(vazquenc)
## maybe str(vazquenc) ; plot(vazquenc) ...



cleanEx()
nameEx("vazquez.example")
### * vazquez.example

flush(stderr()); flush(stdout())

### Name: vazquez.example
### Title: Examples for some analyses
### Aliases: vazquez.example confint intasymm intereven mlik netstats
###   plotmat quant2bin sortmatr sortmatrext
### Keywords: package

### ** Examples

## Not run: 
##D 	data(Safariland)
##D 	
##D 	# confint:
##D 	N100 <- sapply(swap.web(100, Safariland), networklevel, index="nestedness")
##D 	quantile(unlist(N100), c(0.025, 0.975))
##D 	# intasymm: extract values for the asymmetry of interactions and the 
##D 	# dependency matrix for pollinators:
##D 	specieslevel(Safariland)$"higher trophic level"$"interaction push/pull"
##D 	specieslevel(Safariland)$"higher trophic level"$"dependence"
##D 	# for plants:
##D 	specieslevel(Safariland)$"lower trophic level"$"interaction push/pull"
##D 	specieslevel(Safariland)$"lower trophic level"$"dependence"
##D 	
##D 	#intereven
##D 	networklevel(Safariland, index="interaction evenness", intereven="sum")[2]
##D 	# or, as we recommend (see help on networklevel):
##D 	networklevel(Safariland, index="interaction evenness", intereven="prod")[2]
##D 		
##D 	# mlik:
##D 	# calculates the log-likelihood of observing a network, given a probability  
##D 	# matrix of the same size (pweb):
##D 	dmultinom(Safariland>0, prob=pweb, log=TRUE)
##D 	# AIC (the number of parameters is given by how many constraints are put onto the 
##D 	# null model; here, we constrain 9 rows and 27 columns, i.e. sum(dim(binweb))):
##D 	-2*dmultinom(Safariland>0, prob=pweb, log=TRUE) + 2*(sum(dim(binweb)))
##D 	
##D 	# netstats:
##D 	networklevel(Safariland, 
##D 	  index=c("connectance", "interaction evenness", "nestedness", "ISA"))
##D 	mean(specieslevel(Safariland)$"higher trophic level"$"interaction push/pull")
##D 	mean(specieslevel(Safariland)$"lower trophic level"$"interaction push/pull")
##D 	
##D 	#plotmat:
##D 	visweb(t(unname(Safariland)), circles=TRUE, boxes=FALSE)
##D 	
##D 	#sortmatr/sortmatrext:
##D 	sortweb(Safariland, sort.order="inc") #rares species first
##D 	plotweb(sortweb(Safariland, sort.order="dec"), method="normal")
##D 	plotweb(sortweb(web=Safariland, sort.order="seq", 
##D 	  sequence=list(seq.higher=sample(colnames(Safariland)), 
##D 	  seq.lower=sample(rownames(Safariland)))), 
##D 	  method="normal")
##D 	
## End(Not run)



cleanEx()
nameEx("visweb")
### * visweb

flush(stderr()); flush(stdout())

### Name: visweb
### Title: Plotting function to visualize a bipartite food web
### Aliases: visweb
### Keywords: package

### ** Examples


data(Safariland)
visweb(Safariland)
visweb(Safariland, type="diagonal", square="compartment", text="none", 
	frame=TRUE)
visweb(Safariland, type="nested", text="compartment")
visweb(Safariland, circles=TRUE,  boxes=FALSE,  labsize=1, circle.max=3, 
	text="no")
visweb(Safariland,square="b",box.col="green",box.border="red")

#define your colors here,length has to be the numbers of different entries
cols <-0:(length(table(Safariland))-1) 
visweb(Safariland, square="defined", def.col=cols) 



cleanEx()
nameEx("web2edges")
### * web2edges

flush(stderr()); flush(stdout())

### Name: web2edges
### Title: Conversion of a network matrix into a (weighted) edge list
### Aliases: web2edges
### Keywords: package

### ** Examples

data(Safariland)
web2edges(Safariland, return=TRUE)



cleanEx()
nameEx("webs2array")
### * webs2array

flush(stderr()); flush(stdout())

### Name: webs2array
### Title: Puts two or more webs into one array of webs
### Aliases: webs2array
### Keywords: package

### ** Examples

data(Safariland, vazquenc, vazquec)
allin1 <- webs2array(Safariland, vazquenc, vazquec)

# now we can compute distance between two webs:
vegdist(t(cbind(as.vector(allin1[,,1]), as.vector(allin1[,,2]), as.vector(allin1[,,3]))), 
  method="jacc")



cleanEx()
nameEx("wine")
### * wine

flush(stderr()); flush(stdout())

### Name: wine
### Title: Weighted-Interaction Nestedness Estimator
### Aliases: wine plot.wine
### Keywords: package

### ** Examples

data(Safariland, package="bipartite")
safariland.w <- wine(Safariland, 10)
plot.wine(safariland.w)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

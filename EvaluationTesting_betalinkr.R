# Test code for betalinkr ------------------


# load helper function (as long as not in bipartite)  -------------------------
library(bipartite)  # includes vegan
source("webs2array.R")  # existing bipartite function, now also accepting weblist input
source("array2linkmx.R")  # a true helper function; might also be a candidate for cleaning up / speeding up
source("betalinkr.R")  # the main function


# comparisons to original package (evaluation) ------
library(betalink)
betalink(prepare_networks(list(Safariland=Safariland, vazarr=vazarr))[[1]],prepare_networks(list(Safariland=Safariland, vazarr=vazarr))[[2]])
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), binary=T)
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), function.dist="betadiver",index=1)
#--> yep, all 3 the same

# how much do my modifications change the values?
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), partitioning = "adjusted")
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), partitioning = "adjusted", binary=F, proportions=T)


#---  following code also into bipartite testfile ---------------------------------------

# testdata of frame2webs (a case with low overlap and no shared links)
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
  lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
  webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
testarray <- frame2webs(testdata, type.out="array")

betalinkr(testarray) # returns NA for OS and ST
betalink(prepare_networks(list(testarray[,,1]))[[1]],  prepare_networks(list(testarray[,,2]))[[1]]) # betalink also returns NaN for this case
betalinkr(testarray, distofempty="zero") # fixed

# Vazquez-data  
testarray <- webs2array(Safariland, vazarr)  

# comparing various options of betalinkr (evaluation 2.0) ------
# some binary examples
betalinkr(testarray)
betalinkr(testarray, index="jaccard")
betalinkr(testarray, partitioning="adjusted", distofempty="zero")
betalinkr(testarray, partitioning="commondenom", binary=TRUE, distofempty="zero", index="sorensen")
betalinkr(testarray, partitioning="commondenom", binary=TRUE, distofempty="zero", index="jaccard")
# some quantitative examples
betalinkr(testarray, binary=F)
betalinkr(testarray, partitioning="adjusted", binary=FALSE, distofempty="zero")
betalinkr(testarray, partitioning="commondenom", binary=FALSE, distofempty="zero", index="sorensen")
betalinkr(testarray, partitioning="adjusted", binary=FALSE, proportions=TRUE, distofempty="zero")
betalinkr(testarray, partitioning="commondenom", binary=FALSE, proportions=TRUE, distofempty="zero", index="sorensen")
betalinkr(testarray, index="horn", binary=F)


# two fully connected  and completely shared webs
testarray <- array(1:24, dim=c(2, 3, 2))
betalinkr(testarray)

# mostly shared webs
testarray <- array(1:24, dim=c(2, 3, 4))
  testarray[sample(1:24, 10)] <- 0 # setting some entries of above matrix to zero
  testarray <- testarray[, , sample(1:4, 2)] # selecting two sites at random (for now, just developing the function for a 2row-matrix)
betalinkr(testarray)

# example from paper (Figure 1 in Poisot et al. 2012)
metaweb <- matrix(rep(0,25),nrow=5)
metaweb[as.matrix(data.frame(c(2,4,5,5),c(1,2,2,3)))] <- 1  # the adj. matrix for the metaweb
dimnames(metaweb) <- list(letters[1:5], letters[1:5])  # creating species names
web1 <- metaweb
web1[5,2] <- 0 # one link removed from metaweb
web2 <- web1[-1,-1] # top predator removed from web1
web3 <- metaweb[-1,-1] # top predator removed from metaweb

betalinkr_multi(webs2array(web1,web2,web3)) # OS>ST
betalinkr_multi(webs2array(web1,web2,web3), partitioning="adjusted") # OS=ST
betalinkr_multi(webs2array(web1,web2,web3), partitioning="commondenom") # also OS=ST



#----- not for bipartite after here ---------------------------------------------------------------

# output was wrong: bugsearch in webs2array (now fixed) -----------
  # ! problem: output was not same to betalink any more (found reason: webs2array vs. weblist2array; still works in file of 8may)
webs1 <- weblist2array(list(Safariland=Safariland, vazarr=vazarr))
webs2 <- webs2array(list(Safariland=Safariland, vazarr=vazarr))
webs3 <- webs2array(Safariland, vazarr)
sum(webs1!=webs2) # 151
sum(webs1>0) # only 82
sum(webs2>0) # also 82
sum(Safariland); sum(vazarr); sum(webs1); sum(webs2)
which(webs1!=webs2, arr.ind=T) # some differ!
webs1[1,2,1]
webs2[1,2,1]
all(webs1[rownames(Safariland), colnames(Safariland),"Safariland"] == Safariland)
all(webs2[rownames(Safariland), colnames(Safariland),"Safariland"] == Safariland)  # ok, webs2 is wrong!
all(webs3[rownames(Safariland), colnames(Safariland),"Safariland"] == Safariland)  # ok, webs2 
#--> webs2array is wrong!!!
# but why?
webs2array(list(Safariland[1:2,1:2], vazarr[1:2,1:2])) # shows the error: some values recycled, some names duplicated
# came from sapply instead of lapply in all.names assignment, which did wrong unique and thus later wrong indexing
# but the real bug was in indexing, which caused the values to be entered in wrong order!!
webs <- list(Safariland[1:2,1:2], vazarr[1:2,1:2])
webs <- list(Safariland[1:2,1:5], vazarr[1:2,1:4])
webs[[1]]



# calculating dissimilarity (some explorations) ---------------
  # most convenient would be to allow using vegdist, and alternatively betadiver for compatibility with Koleff-framework (24 numbers, used by Poisot)
  # Notes for vegdist (and betadiver) settings
  # NOT RUN
    # equivalents:
    vegdist(linkmx, method = "bray", binary=T) # aka Sorensen
      betadiver(linkmx, method=1)  # Whittaker
    vegdist(linkmx, method = "jaccard", binary=T)
      betadiver(linkmx, method=15)  # "ColwellCoddington"
    # Ruzicka  (preferred by Nacho; same as proportional overlap)
    vegdist(linkmx, method = "jaccard")
    # Bray-Curtis of proportions
    vegdist(decostand(linkmx, method="total") , method="bray")
  # END NOT RUN
  # Note: use index and binary.dist as arguments to pass to either function vegdist/betadiver
      # give warning (or break) if binary=F and betadiver is used
      # give "proportions" as argument (defaulting to TRUE) for binary=F (then using decostand with method="total" first on inputmx)
    

    
    
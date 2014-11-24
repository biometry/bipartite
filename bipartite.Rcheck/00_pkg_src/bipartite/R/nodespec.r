
nodespec <- function(web, inf.replace=NA){
  # an index to describe the functional specialisation of pollinators
  # proposed by Dalgaard et al. 2008 in Oikos
  # this is a purely qualitative measure at present, but could be extended to a
  # qualitative version in principle; in that case, distances should depend on
  # the number of interactions
  #
  # Note: unconnected parts of the network will yield infinite geodesic distances;
  # here we ignore these values (default), thus heavily underestimating distances!
  # The common alternative is to give unconnected paths the value of max(pathlength)+1, 
  # for no other reason than that it has to have some non-infinity value.

  # actually, nodespec is the inverse of closeness. Just apparently nobody has noticed ...


  rr <- nrow(web)
  cc <- ncol(web)

  g <- as.one.mode(web>0) # function of bipartite
  #  if (inf.replace==max1) inf.replace <- Inf
  d <- geodist(g, inf.replace=inf.replace)$gdist
  d[is.infinite(d)] <- max(d[!is.infinite(d)])+1

  meanwithoutself <- function(di, w=rep(1, length(di))){
    # helper function to calculate the mean path length without the path from a vertex to itself
    # di is a vector of path lengths take from the geodist-matrix
    # the 0 in this vector receives a weight of 0
    
    # divide by two! (Fig. 1c shows that they want the two steps (pollinator to plant than on to next pollinator) as one step.    
    zero.position <- which(di==0)                             
    w[zero.position] <- 0
    weighted.mean(di, w=w, na.rm=TRUE)/2
  }
  
  NSlower <- apply(d[1:rr, 1:rr], 2, meanwithoutself) 
  NShigher <- apply(d[(rr+1):(rr+cc), (rr+1):(rr+cc)], 2, meanwithoutself) 

  # minimum specialisation: all polls linked to all plants, hence mean(geodist(.))=1, except with more than 1 compartment!
  
  # maximum specialisation: each pollinator has only one plant to pollinate. For k plants and n pollinators, usually k < n.
  # Hence, we have to allocate the remaining n-k pollinator species to the k plants.
  # To achieve maximum specialisation, we dump all remaining interactions on one pollinator,
  # whose value hence is: (n-k)+1

  # In a maximally specialised (but not compartmented) network, the bipartite plot looks like a zigzag-line. Hence, each pollinator is connected to two others (to the right and left of it), except the two final ones, which have only one. Thus, there are (n-2)*2+2 links in total. The mean distance for any species is thus dependent on its position in the zigzag. Terminal species have mean distance (1+2+3+...+(n-1))/(n-1). The central species will have (1+2+3+..+(n-1)/2) * 2. 


#  if (quantitative) {  # the quantitative version according to CFD
#    # here the "solution" for compartments is to replace infinite path lengths by max.length+1
#    dquant <- geodist(g, inf.replace=inf.replace)$gdist
#    dquant[is.infinite(dquant)] <- max(dquant[!is.infinite(dquant)])+1
#    NSlower <- sapply(1:rr, function(x) meanwithoutself(d[x, 1:rr], w=dquant[x, 1:rr]))
#    NShigher <- sapply((rr+1):(rr+cc), function(x) meanwithoutself(d[x, (rr+1):(rr+cc)], w=1/dquant[x, (rr+1):(rr+cc)]))
#  }

  names(NSlower) <- rownames(web)
  names(NShigher) <- colnames(web)
  list(lower=NSlower, higher=NShigher)

}


#trying to find the maximum value by brute force didn't work. functspec-values are approx. normally distributed (mean=3.5) plus a lot of zeros.
# apply(replicate(100, nodespec(matrix(sample((web>0)*1), nrow=nrow(web), ncol=ncol(web)))$higher), 2, max)
#max(replicate(100, nodespec(matrix(sample((web>0)*1), nrow=nrow(web), ncol=ncol(web)))$higher))
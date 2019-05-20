# A new implementation of network dissimilarity (betalink and related) --------------------
# see helpfile for further info and word doc for Jochen's notes, ToDo and further justification

# main FUNCTION betalinkr --------------------
betalinkr <- function(webarray, index = "bray", binary=TRUE, partitioning="poisot", proportions=FALSE, function.dist="vegdist", distofempty="na"){
  # first, setting all the defaults to Poisot (for comparisons), later change to improved values

  if (class(webarray)=="list") {webarray <- webs2array(webarray)}
  if (dim(webarray)[[3]]!=2) warning("function is designed for a single pair of two webs; unclear output")

  webarray <- webarray[apply(webarray,1,sum)>0, apply(webarray,2,sum)>0, , drop=FALSE] # removing species not observed in either web: improves speed :-)
  
  # for "shared species subweb", set non-shared species links to zero
  array.sharedsp <- webarray
  array.sharedsp[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, , ] <- 0
  array.sharedsp[, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0

  # all links
  linkmx <- array2linkmx(webarray)
  # only links of shared species
  linkmx.sharedsp <- array2linkmx(array.sharedsp)  
  
  # removing links never observed in complete webarray; same subset for both linkmx!
    #links.observed <- colSums(linkmx)>0
    #linkmx <- linkmx[, links.observed]
    #linkmx.sharedsp <- linkmx.sharedsp[, links.observed]
  #-> not using this as it actually made the function slower! (in one test case)
  
  # species community matrix (combining upper and lower of bipartite web)
  specmx.lower <- apply(webarray, c(3,1), sum)
  specmx.higher <- apply(webarray, c(3,2), sum)
  specmx.all <- cbind(specmx.lower, specmx.higher)  # e.g. sites X (plants, pollinators)

  # standardizing to proportions if wanted
  if (proportions){
    if (binary){warning("standardizing to proportions for binary index; do you really want to do this?!?")}
    specmx.all <- decostand(specmx.all, method="total")
    linkmx <- decostand(linkmx, method="total")
    linkmx.sharedsp <- decostand(linkmx.sharedsp, method="total")
  }

  if (partitioning!="commondenom"){
    # alternative subsets of linkmx (partitioning ST and OS) --
      # my approach here is to set the species/links to zero instead of excluding them (which will be done in the index calculation anyways)
      # this makes it easier to match species / links (even without names)
    # shared links of shared species (only LINKS occurring in both sites)
    linkmx.sharedli <- linkmx   
    linkmx.sharedli[, colSums(linkmx.sharedli>0)==1] <- 0
    # varying links of shared species
    linkmx.rewiring <- linkmx.sharedsp - linkmx.sharedli
    linkmx.RewSha <- linkmx.rewiring + linkmx.sharedli  # all links excluding those from unique species
    # links of non-shared species
    linkmx.uniquesp <- linkmx - linkmx.sharedsp
    linkmx.UniSha <- linkmx.uniquesp + linkmx.sharedli # all links excluding rewiring links
  
    # standardizing to proportions if wanted
    if (proportions){
      linkmx.RewSha <- decostand(linkmx.RewSha, method="total")
      linkmx.UniSha <- decostand(linkmx.UniSha, method="total")
    }
  
    # calculating dissimilarity / the betalink components --
    if (function.dist=="vegdist"){
      b_s <- vegdist(specmx.all, method=index, binary=binary) # "S"
      b_wn <- vegdist(linkmx, method=index, binary=binary) # "WN"
      if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # set to conceptually correct value; avoids warning
        b_os.raw <- b_wn
        b_os.raw[] <- 0
      } else {
        b_os.raw <- vegdist(linkmx.RewSha, method=index, binary=binary) # "OS"
      }
      if (distofempty=="zero" & any(rowSums(linkmx.UniSha)==0)){ # set to conceptually correct value; avoids warning
        b_st.raw <- b_wn
        b_st.raw[] <- 0
      } else {
        b_st.raw <- vegdist(linkmx.UniSha, method=index, binary=binary) # "ST"
      }
    } 
    if (function.dist=="betadiver"){
      if (binary==FALSE) {
        warning("betadiver only uses binary data; for quantitative indices use vegdist")
      } else {
        b_s <- betadiver(specmx.all, method=index) # "S"
        b_wn <- betadiver(linkmx, method=index) # "WN"
        if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # adjustment (set to the conceptually correct value); avoids warning
          b_os.raw <- b_wn
          b_os.raw[] <- 0
        } else {
          b_os.raw <- betadiver(linkmx.RewSha, method=index) # "OS"
        }
        if (distofempty=="zero" & any(rowSums(linkmx.UniSha)==0)){ # adjustment (set to the conceptually correct value); avoids warning
          b_st.raw <- b_wn
          b_st.raw[] <- 0
        } else {
          b_st.raw <- betadiver(linkmx.UniSha, method=index) # "ST"
        }
      }
    }
  }
    
  if (partitioning=="commondenom"){
    # here, method must be explicitly specified as one of Sorensen or Jaccard (lower case)
      # quantitative equivalents available with binary=F; I actually use the quantitative formulae, which simplify to binary index if binary data given
    if (binary==TRUE){
      linkmx <- (linkmx > 0)
      linkmx.sharedsp <- (linkmx.sharedsp >0)
    }
    # A, B, C follows Legendre 2014
    # I can calculate all components with linkmx and linkmx.sharedsp, and their difference
    pmins <- pmin(linkmx[1,], linkmx[2,])
    A <- sum(pmins)
    B.tot <- sum(linkmx[1,] - pmins)
    C.tot <- sum(linkmx[2,] - pmins)
    B.rew <- sum(linkmx.sharedsp[1,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,])) # note: frequency-changes of shared interactions also included in rewiring here!
    B.uni <- B.tot - B.rew  # here it works with subtraction
    C.rew <- sum(linkmx.sharedsp[2,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,]))
    C.uni <- C.tot - C.rew
    if (index == "bray") {index <- "sorensen"}  # for convenience, this is also allowed (but avoiding this in first place to be consistent with Legendre-terminology)
    if (index == "sorensen"){denominator <- 2*A + B.tot + C.tot}
    if (index == "jaccard"){denominator <- A + B.tot + C.tot}
    b_wn <- (B.tot + C.tot) / denominator
    b_os <- (B.rew + C.rew) / denominator
    b_st <- (B.uni + C.uni) / denominator 
    b_s <- vegdist(specmx.all, method=switch(index, "jaccard"="jaccard", "sorensen"="bray"), binary=binary) # "S"
  }

  # output (and final steps of calculation)
  if (partitioning=="poisot") {
    b_st.minus <- b_wn - b_os.raw
    return(c(S=b_s, OS=b_os.raw, WN=b_wn, ST=b_st.minus))
  }
  if (partitioning=="adjusted") {
    b_os <- b_os.raw * b_wn / (b_os.raw + b_st.raw)  # the correction I propose to apply
    b_st <- b_st.raw * b_wn / (b_os.raw + b_st.raw)  # equivivalently for st
    return(c(S=b_s, OS=b_os, WN=b_wn, ST=b_st))
  }
  if (partitioning=="commondenom") {
    return(c(S=b_s, OS=b_os, WN=b_wn, ST=b_st))
  }
}

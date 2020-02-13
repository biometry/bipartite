# A new implementation of network dissimilarity (betalink and related) --------------------
# see helpfile for further info and word doc for Jochen's notes, ToDo and further justification

# main FUNCTION betalinkr --------------------
betalinkr <- function(webarray, index = "bray", binary=TRUE, partitioning="commondenom", proportions=!binary, function.dist="vegdist", distofempty="zero", partition.st=FALSE, partition.rr=FALSE){
  
  if (class(webarray)=="list") {webarray <- webs2array(webarray)}
  if (dim(webarray)[[3]]!=2) warning("function is designed for a single pair of two webs; unclear output")

  webarray <- webarray[apply(webarray,1,sum)>0, apply(webarray,2,sum)>0, , drop=FALSE] # removing species not observed in either web: improves speed :-)

  # standardizing to proportions if wanted (now on webarray with marginal totals; decostand only used for method!=commondenom)
  if (proportions){
    if (binary){warning("standardizing to proportions for binary index; do you really want to do this?!?")}
    array.of.sums <- webarray
    array.of.sums[,,1] <- sum(webarray[,,1])
    array.of.sums[,,2] <- sum(webarray[,,2])
    webarray <- webarray / array.of.sums
  }

  if (binary==TRUE){
    webarray[webarray>0] <- 1
  }  

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
    # now ensuring a single entry per species also in named foodwebs (adding ingoing and outgoing links, instead of counting them separately)
  specmx.lower <- apply(webarray, c(3,1), sum)
  specmx.higher <- apply(webarray, c(3,2), sum)
  specmx.higher.unique <- specmx.higher[, !(colnames(specmx.higher) %in% colnames(specmx.lower))]
  if (is.null(colnames(specmx.higher))) {specmx.higher.unique <- specmx.higher}  # if no species names are given, assuming bipartite webs
  specmx.all <- cbind(specmx.lower, specmx.higher.unique)  # e.g. sites X (plants, pollinators); 
  duplicolnames <- setdiff(colnames(specmx.higher), colnames(specmx.higher.unique))
  specmx.all[, duplicolnames] <- specmx.all[, duplicolnames] + specmx.higher[, duplicolnames]

  if (partitioning=="poisot"){
    if (partition.st | partition.rr){warning("further partitioning only available with method partitioning='commondenom'")}
    
    # alternative subsets of linkmx (partitioning ST and OS) --
      # my approach is to set the species/links to zero instead of excluding them (which will be done anyways when dissimilarity is calculated)
      # this makes it easier to match species / links (even without names)
    # shared links of shared species (only LINKS occurring in both sites)
    linkmx.sharedli <- linkmx   
    linkmx.sharedli[, colSums(linkmx.sharedli>0)==1] <- 0
    # varying links of shared species
    linkmx.rewiring <- linkmx.sharedsp - linkmx.sharedli
    linkmx.RewSha <- linkmx.sharedsp  # all links excluding those from unique species
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
      b_zero <- b_wn  # preparation for "distofempty": first get structure of distance matrix
      b_zero[] <- 0   # preparation for "distofempty": define a zero distance matrix
      if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # set to conceptually correct value; avoids warning
        b_os.raw <- b_zero # no shared species means zero contribution of OS
      } else {
        b_os.raw <- vegdist(linkmx.RewSha, method=index, binary=binary) # "OS"
      }
    } 
    if (function.dist=="betadiver"){
      if (binary==FALSE) {
        warning("betadiver only uses binary data; for quantitative indices use vegdist")
      } else {
        b_s <- betadiver(specmx.all, method=index) # "S"
        b_wn <- betadiver(linkmx, method=index) # "WN"
        if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # set to the conceptually correct value; avoids warning
          b_os.raw <- b_zero # no shared species means zero contribution of OS
        } else {
          b_os.raw <- betadiver(linkmx.RewSha, method=index) # "OS"
        }
      }
    }
    # output (and final steps of calculation)
    b_os <- b_os.raw
    b_st <- b_wn - b_os.raw
    return(c(S=b_s, OS=b_os, WN=b_wn, ST=b_st))
  }
    
  if (partitioning=="commondenom"){
    # here, index must be explicitly specified as one of Sorensen or Jaccard (lower case)
      # quantitative equivalents available with binary=F; I actually use the quantitative formulas, which simplify to binary index if given binary data 

    # A, B, C follows Legendre 2014; index "tot"=total differences, "rew"=rewiring differences, "uni"=unique species differences
    # I can calculate all components with linkmx and linkmx.sharedsp, and their difference
    pmins <- pmin(linkmx[1,], linkmx[2,])
    A <- sum(pmins)
    B.tot <- sum(linkmx[1,] - pmins)
    B.rew <- sum(linkmx.sharedsp[1,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,])) # note: frequency-changes of shared interactions also included in rewiring here!
    B.uni <- B.tot - B.rew  # here it works with subtraction
    C.tot <- sum(linkmx[2,] - pmins)
    C.rew <- sum(linkmx.sharedsp[2,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,]))
    C.uni <- C.tot - C.rew
    if (index == "bray") {index <- "sorensen"}  # for convenience, bray is also allowed (but avoiding this in first place to be consistent with Legendre-terminology)
    if (index == "sorensen"){denominator <- 2*A + B.tot + C.tot}
    if (index == "jaccard"){denominator <- A + B.tot + C.tot}
    b_wn <- (B.tot + C.tot) / denominator
    b_os <- (B.rew + C.rew) / denominator
    b_st <- (B.uni + C.uni) / denominator 
    b_s <- vegdist(specmx.all, method=switch(index, "jaccard"="jaccard", "sorensen"="bray"), binary=binary) # "S"
    if (partition.st){
      # partition ST further into ST.l, ST.h and ST.lh
      # this seems easier to do when calculating with the arrays first, converting to linkmx only later
      # lh = all - sharedhighORlow
      # l = all - sharedlow - sharedboth - lh
      # h = all - sharedhigh - sharedboth - lh
      array.sharedlow <- webarray
      array.sharedlow[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, , ] <- 0
      array.sharedhigh <- webarray
      array.sharedhigh[, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0
      array.sharedhighORlow <- webarray
      array.sharedhighORlow[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0
      linkmx.lh <- array2linkmx(webarray - array.sharedhighORlow)  # interactions with both partners unique to one of the webs
      linkmx.l <- array2linkmx(webarray - array.sharedlow) - linkmx.lh # interactions with lower sp unique to one of the webs (but higher species shared)
      linkmx.h <- array2linkmx(webarray - array.sharedhigh) - linkmx.lh # interactions with higher sp unique to one of the webs (but lower species shared)
      # Note: the pmin-part can be omitted here, as it will always be zero!
      B.l <- sum(linkmx.l[1,])
      B.h <- sum(linkmx.h[1,])
      B.lh <- sum(linkmx.lh[1,])
      # B.l + B.h + B.lh == B.uni  # check, must be true
      C.l <- sum(linkmx.l[2,])
      C.h <- sum(linkmx.h[2,])
      C.lh <- sum(linkmx.lh[2,])
      # C.l + C.h + C.lh == C.uni  # check, must be true
      b_st.l <- (B.l + C.l) / denominator 
      b_st.h <- (B.h + C.h) / denominator 
      b_st.lh <- (B.lh + C.lh) / denominator 
    }
    if (partition.rr){
      # partition WN and OS further into link Replacement component and link Richness difference component
        # for binary=FALSE, link richness difference component is the dissimilarity component due to differences in network totals
      if (proportions==TRUE){warning("partitionining into replacement and richness (abundance) difference components may be meaningless with proportions")}
      b_wn.repl <- 2*min(B.tot, C.tot) / denominator
      b_os.repl <- 2*min(B.rew, C.rew) / denominator
      b_wn.rich <- abs(B.tot - C.tot) / denominator
      b_os.rich <- abs(B.rew - C.rew) / denominator
    }
    # output (concetenated for secondary partitionings)
    output <- c(S=b_s, OS=b_os, WN=b_wn, ST=b_st)
    if (partition.st==TRUE){
      output <- c(output, ST.l=b_st.l, ST.h=b_st.h, ST.lh=b_st.lh)
    }
    if (partition.rr==TRUE){
      output <- c(output, WN.repl=b_wn.repl, OS.repl=b_os.repl, WN.rich=b_wn.rich, OS.rich=b_os.rich)
    }
    return(output)
  }
}

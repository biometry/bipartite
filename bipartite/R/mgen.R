mgen <- function(web, n=sum(web), keep.species=TRUE, rep.cell=TRUE, autotransform="sum_warn"){
	# function to generate a quantitative network based on a probability matrix
	# originally by Diego Vazquez, revised by CFD and greatly revised by JF (Jan 2021)
	# web 	a matrix with observation probabilities, emerging from some null model considerations external to this function; if an original network is used, this will be automatically converted to a probability matrix by dividing it by the number of interactions (CFD); ORIGINAL: a probability matrix
	# n     number of interactions to allocate into the new matrix
	# autotransform:    determines how a non-probability web is converted into probabilities; 
	#			option 1: "sum" or "sum_warn": simply divide each entry by the sum of interactions in the web
	#			option 2: "equiprobable": product of marginal probabilities (col/rowSums divided by sum(web) (in a sense this is the basis of the r2dtable null model, just without the 'turn into integers' bit)
	# keep.species: Random assignment of interactions may lead to empty columns or rows and hence reduce the dimensions of the simulated web. By default, this is prevented from happening, i.e. each row/column will receive at least one interaction. Setting keep.species to FALSE may (but need not) cause a loss of species.

  if (!isTRUE(all.equal(sum(web), 1)) & autotransform %in% c("sum", "sum_warn")) { # added by CFD, changed by JF (now tests for approximate equality of sum(web) to 1, as this may not be exactly true even after dividing by sum(web))
  	if (autotransform == "sum_warn") {
  	  warning(paste("This is not a probability matrix! I will proceed after transforming the entries as in 'autotransform'='sum'!"))
  	}
		m <- web/sum(web) # convert to probabilities
  } else m <- web
  if (autotransform == "equiprobable"){ # an odd name for expectation from marginal totals (nullmodelish)
    m <- (rowSums(web)/sum(web)) %*% t(colSums(web)/sum(web))
  }
  
  if (rep.cell == FALSE & n > sum(m>0)){ # JF: changed this to checking for n not exceeding no. of nonzero cells!
    stop("Argument n should not exceed the number of cells in matrix (with non-zero probability)!") # changed from message to error
  }
  mint <- matrix(0, nrow(m), ncol(m)) # Interaction matrix
  if (keep.species){
    transp.mx <- (nrow(mint) < ncol(mint))
    if (transp.mx) {
      # turn the problem into one where nrow is never smaller than ncol
      # (makes the first step more effective, and the removal of excess links easier to program)
      mint <- t(mint)
      m <- t(m)
    }
    for (i in 1:nrow(m)){
      mint[i, sample(ncol(m), size = 1, prob=m[i,])] <- 1 # changed JF
    }
    for (j in 1:ncol(m)){
      if(sum(mint[,j]) == 0){
        mint[sample(nrow(m), size=1, prob=m[,j]), j] <- 1 # changed JF
      }
    }
    if (any(rowSums(mint)>1) & any(colSums(mint)>1)){
      # new part added by JF: try to remove non-essential links; not very effective, but improves a bit
      # this criterion means too many links assigned; at least one of rs or cs should have only 1s
      # the awkward double transform (num->char->num avoids sampling from 1:n if which() returns a length 1 element)
      for (i.remove in 1:(sum(mint) - nrow(mint))){ # do exactly one attempt for each excess link
        myrow <- as.numeric(sample(as.character(which(rowSums(mint)>1)), size=1))
        validcols <- which(colSums(mint) > 1 & mint[myrow,] > 0)
        if (length(validcols) > 0){ 
          # Note: sometimes  problem cannot be solved by removing a link, swap would be needed: ignored
          mycol <- as.numeric(sample(as.character(validcols), size=1))
          mint[myrow, mycol] <- mint[myrow, mycol] - 1
        }
      }
    }
    if (transp.mx) { # back-transpose matrices, if transposed before in keep.species
      mint <- t(mint)
      m <- t(m)
    }
  }
  # new main sampling part: old while-loop replaced by using sample once
  if (n > sum(mint)){
    if (rep.cell == FALSE){m[mint == 1] <- 0}
    intvec <- sample(1:length(mint), size = n - sum(mint), replace = rep.cell, prob = m)
    for (i in intvec){mint[i] <- mint[i] + 1}
  }
  dimnames(mint) <- dimnames(web)  # added JF
  return(mint)
}
# mgen(web=Safariland)
# mgen(web=Safariland, autotransform="equiprobable")
# mgen(web=Safariland/sum(Safariland), n=sum(Safariland), keep.species=FALSE)
# mgen(web=Safariland/sum(Safariland), n=200, rep.cell=F)
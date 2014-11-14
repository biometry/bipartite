swap.web <- function(N, web, verbose=FALSE, c.crit=1e4){
  # function to generate null model webs under the following constraints:
  # 1. marginal totals are identical to those observed (as in r2dtable)
  # 2. connectance is as observed (as in shuffle.web)
  #
  # This is a rather constrained nullmodel!
  

  #--------- Helper functions -----------

  findmat.full <- function(web, c.crit=1e4){ 
    # helper function to find a 2 x 2 matrix with non-zero diagonal entries and at least one 0 on off-diagonal!
    #  1. while-loop: find ONE column with at least 2 entries
    #  2. while-loop: find a SECOND column with at least 2 entries
    #  3. if (sum(pmin(a>0, b>0)) != 2) go back to 1.
    #  4. return this matrix

    csums = 0
    counter = 0
    while (csums < 2 & counter < c.crit){
        c1sum = 0; c2sum = 0
        while (c1sum < 2 & counter < c.crit){  # find first column
          ca <- sample(1:ncol(web), size=1)
          c1 <- web[,ca]
          c1sum <- sum(c1>0)
          counter=counter+1
        }
        cb = ca
        while ( (c2sum < 2 | cb==ca) & counter < c.crit){ # find second column
          cb <- sample(1:ncol(web), size=1)
          c2 <- web[,cb]
          c2sum <- sum(c2>0)
          counter=counter+1
        }
        csums <- sum(pmin(web[,ca]>0, web[,cb]>0)) # check if the sum of parallel minima is at least 2
        counter=counter+1
    }
    if (counter >= c.crit){ # if no pair of columns can be found, set matrix to NULL
        mat <- NULL 
    } else {
      r <- as.numeric(which(pmin(web[,ca]>0, web[,cb]>0)>0))[sample(csums, 2, replace=FALSE)]
      # (last line:) from those columns with two entries select randomly 2 as rows
      mat <- web[r, c(ca, cb)]
      attr(mat, "rows") <- r
      attr(mat, "cols") <- c(ca, cb)
    }
#    previous version:
#    mat <- matrix(0, 2, 2)
#    while( sum(mat==0) > 0 ){
#      rselect <- sample(1:nrow(web), size=2, replace=FALSE)
#      cselect <- sample(1:ncol(web), size=2, replace=FALSE)
#      mat <- web[rselect, cselect]
#    }
#    attr(mat, "rows") <- rselect
#    attr(mat, "cols") <- cselect

    mat
  }

  findmat.empty <- function(web){
    # helper function to find a 2 x 2 matrix with non-zero off-diagonal entries
    mat <- matrix(0, 2, 2)
    while( any(diag(mat) == 0) | (all(c(mat[1,2], mat[2,1]) != 0)) ){
      rselect <- sample(1:nrow(web), size=2, replace=FALSE)
      cselect <- sample(1:ncol(web), size=2, replace=FALSE)
      mat <- web[rselect, cselect]
    }
    attr(mat, "rows") <- rselect
    attr(mat, "cols") <- cselect
    mat
  }


  downswap <- function(first, m, n, c.crit=1e4){
      second <- first  # to have starting conditions when the algorithm gets trapped
      while (m < n){
        mat <- findmat.full(first, c.crit=c.crit) # yields a 2 x 2 matrix with non-zero diagonal       
        while (is.null(mat)){
            first <- r2dtable(1, r=rowSums(second), c=colSums(second))[[1]]
            n <- sum(first>0)
            mat <- findmat.full(first)
            cat("New null matrix needed: old one sucked.\n")
        } 
        # cat("2x2 matrix found ")
        # swap:
        mat.new <- mat
        diag(mat.new) <- (diag(mat) - min(diag(mat)))
        diag(mat.new[2:1, ]) <- diag(mat[2:1,]) + min(diag(mat))
  		# check that we are not downswapping more than interactions still available:  
        trial <- first
        trial[attr(mat, "rows"), attr(mat, "cols")] <- mat.new
		if (sum(trial > 0) < m) next 
        if (sum(trial > 0) >= m) first[attr(mat, "rows"), attr(mat, "cols")] <- mat.new
	    n <- sum(first > 0)
        #cat(paste(n, "vs.", m, "\n", sep=""))
      }
      first
  }

  upswap <- function(first, m, n){
      if (m==n) {return(first); stop()}
      while (m > n){
        mat <- findmat.empty(first) # yields a 2 x 2 matrix with non-zero diagonal and at least one 0 on off-diagonal
        # swap:
        mat.new <- mat
        howmuchchange <- sample(min(diag(mat)), 1)
        diag(mat.new) <- diag(mat)-howmuchchange
        diag(mat.new[2:1, ]) <- diag(mat[2:1,]) + howmuchchange
        first[attr(mat, "rows"), attr(mat, "cols")] <- mat.new
        n <- sum(first > 0)
      }
      first
  }

  nullmaker <- function(web, verbose=TRUE, c.crit){ 
    m <- sum(web>0) # counts number of non-empty cells
    first <- r2dtable(1, r=rowSums(web), c=colSums(web))[[1]] # creates a null web with same marginal totals
    n <- sum(first>0)                         
    if (verbose) if (m > n) cat("Demands filling algorithm!\n") else cat("Demands emptying algorithm!\n")
    if (m < n)  null <- downswap(first, m=m, n=n, c.crit=c.crit)
    if (m >= n) null <- upswap(first, m=m, n=n)  
    # sum(null>0); sum(null) # for testing purposes
    null
  }

 #---------------------------s---------------------------------------------------
 # main part of the function:
#  nullmaker(barrett1987) #test
    
  nulls <- replicate(N, nullmaker(web, verbose=verbose, c.crit=c.crit), simplify=FALSE)
  
  nulls
}

#data(Safariland)
#swap.web(N=2, web=Safariland)
#system.time(swap.web(10, Safariland))
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


  downswap <- function(first, m, n, c.crit=1e4, swap.method="random", verbose=TRUE){
      second <- first  # to have starting conditions when the algorithm gets trapped
	  replace.counter <- 0
      while (m < n){
        mat <- findmat.full(first, c.crit=c.crit) # yields a 2 x 2 matrix with non-zero diagonal       
        while (is.null(mat)){
            first <- r2dtable(1, r=rowSums(second), c=colSums(second))[[1]]
            n <- sum(first>0)
            mat <- findmat.full(first)
            if (verbose) cat("New null matrix computed: old one had no manipulatable 2x2 entries.\n")
            replace.counter <- replace.counter + 1
            # if (replace.counter > 49) {
            	# swap.method <- "not.random"
            	# if (verbose) cat("Changed swapping algorithm due to cul-de-sac.\n")
            # }
            # print(replace.counter)
            if (replace.counter >= 10) return(NULL)
        } 
        if (replace.counter >= 10) return(NULL)
        # cat("2x2 matrix found ")
        # swap.mat <- function(mat, method="random"){
        	mat.new <- mat
        	 if (swap.method=="random"){
        		# new version, using always one entry to add to the other, irrespective of its size:
        		mat.new[1, 1] <- mat[1,1] + mat[2,2]
        		mat.new[2, 2] <- 0
        		# Potentially, this can cause problems, in small matrices, when a large cell entry is move to a small one, thus making further manipulations impossible. In this case, I guess it may be fair to move back to the old algorithm?
        	} else {
        		# old version, using always the smallest entry to add to the larger entry:
        		diag(mat.new) <- (diag(mat) - min(diag(mat)))
        		diag(mat.new[2:1, ]) <- diag(mat[2:1,]) + min(diag(mat))
        	}
        	#return(mat.new)
        #}
        #if (is.null(mat)){
        #	first <- r2dtable(1, r=rowSums(second), c=colSums(second))[[1]] 
        #	swap.method <- "not.random"
        #} else { mat.new <- swap.mat(mat, method=swap.method) }
        
  		# check that we are not downswapping more interactions than still available:  
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
    if (m < n)  {
    	null <- downswap(first, m=m, n=n, c.crit=c.crit, swap.method="random", verbose=verbose)
    	# for some small matrices, this will not work; then we go back to the non-random moving of cells in the downswap:
    	if (is.null(null))  null <- downswap(first, m=m, n=n, c.crit=c.crit, swap.method="not.random", verbose=verbose)
    }
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
#swap.web(N=2, web=Safariland, verbose=TRUE)
#system.time(swap.web(10, Safariland))
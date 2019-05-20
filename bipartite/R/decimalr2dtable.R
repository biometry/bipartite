decimalr2dtable <- function(N=10, web, steps=prod(dim(web))){
  # function to generate null models under the following ideas:
  # 1. ideally a model using the original data, rather than rates, should be used
  # 2. also works for non-integers
  # 3. no constraint on connectance
  # 4. test: should be very similar to r2dtable for integers

  # algorithm
  #1. randomly choose two rows and two columns
  #2. compute min on diagonal 
  #3. draw random value from uniform between 0 and min
  #4. subtract this value from diagonal elements and add to counter-diagonal
  #5. repeat prod(dim(web))*10 times
  
  if (any(as.vector(web) < 0)) stop("This function can only handle a non-negative matrix.")
  
  oneshifter <- function(web){
    # 1. 
    minval <- 0
    while(minval == 0){
      rr <- sample.int(nrow(web), size=2, replace=F)
      rc <- sample.int(ncol(web), size=2, replace=F)
      # 2.
      minval <- min(c(web[rr[1], rc[1]], web[rr[2], rc[2]]))
    }
    # 3. 
    shiftthis <- runif(1, 0, minval)
    # 4. 
    web[rr[1], rc[1]] <- web[rr[1], rc[1]] - shiftthis
    web[rr[2], rc[2]] <- web[rr[2], rc[2]] - shiftthis
    
    web[rr[1], rc[2]] <- web[rr[1], rc[2]] + shiftthis
    web[rr[2], rc[1]] <- web[rr[2], rc[1]] + shiftthis
    
    return(web)
  } 
  
  manyshifter <- function(web, steps=prod(dim(web))*10){
    for (i in 1:steps){
      web <- oneshifter(web)
    }
    return(web)
  }
  
  res <- replicate(N, manyshifter(web, steps=steps), simplify=FALSE)
  
  return(res)
}

#library(bipartite)
#library(profvis)
#profvis(decimalr2dtable(100, Safariland))

# system.time(nulls <- decimalr2dtable(100, Safariland)) 
# g.dec <- sapply(nulls, networklevel, index="generality")
# obs <- networklevel(Safariland, index="generality")
# system.time(nullsint <- nullmodel(Safariland, N=1000))
# g.int <- sapply(nullsint, networklevel, index="generality")
# plot(density(g.dec[1,]), xlim=c(1, 3))
# lines(density(g.int[1,]), col="red")
# abline(v=obs[1], col="green")


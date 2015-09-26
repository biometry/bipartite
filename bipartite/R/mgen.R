mgen <- function(web, n=sum(web), keep.species=TRUE, rep.cell=TRUE, autotransform="sum", trials=100){
	# function to generate a quantitative network based on a probability matrix
	# by Diego Vazquez (brushed up for a more bipartite-consistent use of names by CFD) 
	# web 	a matrix with observation probabilities, emerging from some null model considerations external to this function; if an original network is used, this will be automatically converted to a probability matrix by dividing it by the number of interactions (CFD); ORIGINAL: a probability matrix
	# n     number of interactions to allocate into the new matrix
	# autotransform:    determines how a non-probability web is converted into probabilities; 
	#			option 1: "sum": simply divide each entry by the sum of interactions in the web
	#			option 2: "equiprobable": product of marginal probabilities (col/rowSums divided by sum(web) (in a sense this is the basis of the r2dtable null model, just without the 'turn into integers' bit)
	# keep.species: Random assignment of interactions may lead to empty columns or rows and hence reduce the dimensions of the simulated web. By default, this is prevented from happening, i.e. each row/column will receive at least one interaction. Setting keep.species to FALSE may (but need not) cause a loss of species.
	# trials: allocating interactions (when rep.cell=TRUE) can be difficult or even impossible. When the number of attempts to allocate them exceeds trials * n it will be stopped and an error message will be returned. Default is 100. Setting 'trials' to a very large value may make this function run for hours and hours. Your choice!
  
  if (sum(web) != 1) { # added by CFD
  	message(paste("This is not a probability matrix! I will proceed after transforming the entries according to option 'autotransform':", autotransform, "!"))
  	if (autotransform == "sum") {
  		m <- web/sum(web)
  	} else {# equiprobable, or anything else
  		m <- (rowSums(web)/sum(web)) %*% t(colSums(web)/sum(web))
  	}
  } else m <- web
  
  if (rep.cell == FALSE & n > (nrow(m)*ncol(m))){
    message("Argument n should be smaller than the number of cells in matrix!")
  }
  else{
    mac <- matrix(cumsum(m),nrow(m),ncol(m)) #Cumulative probability matrix
    mint <- matrix(0,nrow(m),ncol(m)) #Interaction matrix
    if (keep.species){
      for (i in 1:nrow(m)){
        c1 <- sample(ncol(m), replace=TRUE, prob=colSums(m))
        c1 <- c1[1]
        mint[i, c1] <- 1
      }
      for (i in 1:ncol(m)){
        if(sum(mint[,i]) == 0){
          r1 <- sample(nrow(m), replace=TRUE, prob=rowSums(m))
          r1 <- r1[1] 
          mint[r1, i] <- 1
        }
      }
    }
    while.counter <- 0
    while (sum(mint) < n){
      rand <- runif(1, 0, 1)
      ri <- min(which(mac >= rand))
      if (rep.cell == TRUE) mint[ri] <- mint[ri] + 1
      if (rep.cell == FALSE) mint[ri] <- 1
      while.counter <- while.counter + 1
      if (while.counter >= trials*n) stop("Cannot allocate the requested interactions in a reasonable amount of time! \n Either increase 'trials' or decrease 'n'.")
    }
    mint
  }
}
# mgen(web=Safariland)
# mgen(web=Safariland, autotransform="equiprobab")
# mgen(web=Safariland/sum(Safariland), n=sum(Safariland), keep.species=FALSE)
# mgen(web=Safariland/sum(Safariland), n=200, rep.cell=F)
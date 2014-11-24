PDI <- function(web, normalise=TRUE, log=FALSE){
  # function implementing the Paired Difference Index proposed by Poisot et al. 2011:
  # log: this option was my own addition, since often the absolute values are extremely skewed
  # author: Carsten F. Dormann 10 Aug 2011
  web <- as.matrix(web)
  N <- NCOL(web)
  H <- NROW(web)
  PDIj <- web[1,] # to automatically get the names, too.
  
  if (log) web <- log(web+1)
  
  for (i in 1:N){
    Ps <- web[,i]
    if (normalise) Ps <- Ps/max(web[,i]) # compute relative abundances
    sort.Ps <- sort(Ps, decreasing=TRUE) # sort by value
    P1 <- sort.Ps[1] # P1 is maximum value
    PDIj[i] <- sum(P1 - sort.Ps[-1]) / (H - 1)
  }
  return(PDIj)
}

#data(Safariland)
#PDI(Safariland) # specialisation of pollinators
#PDI(t(Safariland)) # specialisation of plants
#PDI(t(Safariland), log=TRUE)
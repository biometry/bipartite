dfun <- function(web, abuns=NULL){     # abuns is external data on supply of "resources", e.g. abundance of higher trophic level
  if (is.null(abuns)) web <- empty(web)
  if (!is.null(abuns) & length(abuns)!= ncol(web)) stop("Length of abundance vector and number of higher trophic level species do not match!")
  if (is.null(abuns)) {q <- colSums(web) / sum(web)} else {q <- abuns / sum(abuns)}                 # q[j], the proportion of each "resource"
  cs <- colSums(web)

  # function for d
  d <- function(x, q=q){                      # q = normalized abuns
    p <- x/sum(x)                            # p'[ij], m cancels out of the formula
    sum(p[p!=0] * log(p[p!=0] / q[p!=0]))
  }

  # function for dmin
  dminfind <- function(x, q){               # x is the actual vector of resource utilization (one row of the web)
    expec <- floor(q * (sum(x)))                # fill in downrounded expected values
    restuse <- sum(x) - sum(expec)
    x.new <- expec                              # new vector/distribution of x
	  if (!is.null(abuns)) {i.vec <- 1:length(x)}   # i.vec will be used in the below (cells which are not "full")
    for (j in 1:restuse) {                      # add rest in single steps
      if (is.null(abuns)) {i.vec <- which(x.new < cs)}    # i.vec will be used in the next step (cells which are not "full" due to colSums reached)
      # looking for the cell where 1 can be added with the smallest increase in KLD = d
      d.check <- numeric(0)
      xsum <- sum(x.new)
      #p.0 <- x.new / xsum
      p.1 <- x.new / (xsum + 1)
      dstep.1 <- ifelse(x.new!=0, p.1 * log(p.1 / q), 0)         # all elements of new d except the one were 1 is added
      for (i in i.vec) {
        pi.1 <-  (x.new[i] + 1) / (xsum + 1)
        d.check[i] <- pi.1 * log(pi.1 / q[i]) + sum(dstep.1[-i])
      }
      i.best <- which.min(d.check)[1]
      x.new[i.best] <- x.new[i.best] + 1          # add 1 on the first possible cell with smallest KLD = d
    }
    return(d(x.new, q))
  }

  # function for dmax (two cases)
  dmaxfind <- function(x, q) {
  	  if (min(q) == 0) minq <- min(q[q > 0]) else minq <- min(q)
      dmax <- ifelse(is.null(abuns), log(sum(web)/sum(x)), log(1/minq))
      return(dmax)
  }
  # Here is the problem:
  # Consider a matrix n2 like this:
  # (n2 <- matrix(c(8,3,1,6,2,0,2,0,0,1,0,0,1,0,0,0,0,1), ncol=3, byrow=T)  )
  # with marginal totals
  # rowSums(n2); colSums(n2)
  # Now consider species [,1]: We need to allocate 18 interactions, but there is no combination of row totals to yield this number. That means, at least one interaction has to go into a row in which also another higher-level species has interactions. Thus, species [,1] cannot be perfectly specialised.
  # (Perfect specialisation would mean that it monopolises the species it interacts with.)
  # Another problem:
  # In the transposed matrix t(n2), species [,6] has only interaction, but the minimum marginal total is 2. Thus, again, [,6] doesn't have enough interactions to monopolise its partner.
  # As a consequence of both problems, the dmax calculate above does not hold! dmax is actually lower.
  # We could not come up with a heuristic way to create a maximally specialised web, so we leave dfun as it is.
  # (Nico found and Jochen explained this problem to me (CFD). Jochen also implemented the correct dmax when abundances are given.)
  

  d.raw <- apply(web, 1, d, q)
  d.min <- apply(web, 1, dminfind, q)
  d.max <- apply(web, 1, dmaxfind, q)

  # sometimes allocating integers will not find the minimal solution; in these cases, dmin will be replaced by the actual value of d.raw, should that be lower:
  d.min <- ifelse(d.raw < d.min, d.raw, d.min) # (this line added 31-Mar-2011)

  d.prime <- (d.raw - d.min)/(d.max - d.min)
  list("dprime"=d.prime, "d"=d.raw, "dmin"=d.min, "dmax"=d.max)
}

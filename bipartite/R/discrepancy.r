#P <- matrix(c(rep(1,10),0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,0,0,0,0), ncol=6, nrow=5, byrow=T)
#A <- matrix(c(rep(1,9),0,0,1, 1,1,0,1,0,0, 1,0,1,1,0,0, 1,1,0,0,0,0), ncol=6, nrow=5, byrow=T)

discrepancy <- function(mat){
    # discrepancy of a matrix sensu Brualdi & Sanderson 1999
    # there is already an implementation of this function in vegan (called: nesteddisc), but it does not deliver the "right" answer on the example above (taken from Brualdi & Sanderson 1999)
    mat <- unname(as.matrix(mat>0))*1
    # first, pack matrix by sorting, yielding matrix A: (code from nestedness.corso)
    A <- mat[order(rowSums(mat), decreasing=TRUE), order(colSums(mat), decreasing=TRUE), drop=FALSE]

    # second, pack matrix maximally, yielding P:
    P <- A
    P[,] <- 0
    rs <- rowSums(A)
    if (suppressWarnings(any(rs) != 0)){
      for (i in 1:nrow(A)){
          P[i,1:rs[i]] <- 1
      }
    }
    # this matrix P has the same row sums as A!
       
    c("discrepancy"=sum(P!=A)/2)
}
#discrepancy(A)
#commsimulator(A, method="quasiswap") # a null model for binary, maintaining row and column totals
# example:
#data(Safariland)
#nulls <- replicate(1000, discrepancy(commsimulator(Safariland, method="quasiswap")))
#hist(nulls)
#obs <- discrepancy(Safariland)
#abline(v=obs, lwd=3, col="grey")
#c("p value"=min(sum(nulls>obs), sum(nulls<obs))/length(nulls))
## calculate Brualdi & Sanderson's Na-value (i.e. the z-score):
#c("N_a"=(unname(obs)-mean(nulls))/sd(nulls))
## this value is t-distributed: pt(0.57, 1, lower.tail=F) # 0.335, i.e. not significant
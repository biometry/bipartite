vaznull <- function (N, web) 
{
  # This function produces a null model network with two constraints:
  # a) marginal totals are the same as in the original network (see also r2dtable)
  # b) connectance is the same as in the original network.
  # vaznull differs from swap.web both in the algorithm used as well as in the
  # null model it outputs.
  # vaznull is our implementation of the algorithm propose by Diego V·zquez, hence its name.
  # While vaznull is slower than swap.web, we regard it as better.
  #
  # implementor: Bernd Gruber <bernd.gruber@ufz.de> & Carsten Dormann

    web <- as.matrix(empty(web)) # otherwise we cannot compute the crossproduct later on. 
    
    vaznull.fast <- function(web) {
    	
        rs.p <- rowSums(web)/sum(web)
        cs.p <- colSums(web)/sum(web)
        P <- P1 <- tcrossprod(rs.p, cs.p)
        #create only one new mat
        finalmat <- matrix(0, nrow(web), ncol(web))
        n.int.finalmat <- 0
        while (n.int.finalmat < sum(dim(web))) {
            # randomly select a cell, according to probability matrix P:
            sel <- sample(1:length(web), 1, prob = P, replace = TRUE)
            #convert sel to coordinates:
            selc <- floor((sel - 1)/(dim(web)[1])) + 1
            selr <- ((sel - 1)%%dim(web)[1]) + 1
            #check if row or column is still free?
			if (sum(finalmat[, selc]) == 0 | sum(finalmat[selr, 
                ]) == 0) {
                finalmat[sel] <- 1
                P[sel] <- 0
            }
            n.int.finalmat <- sum(rowSums(finalmat) > 0) + sum(colSums(finalmat) > 
                0)
        }
        
      # Step 2. IV: fill up to full connectance:
        conn.remain <- sum(web > 0) - sum(finalmat > 0)
        if (conn.remain > 0) {
            add <- sample(which(finalmat == 0), conn.remain, 
                prob = P1[finalmat == 0])
            finalmat[add] <- 1
        }
        # FILL INTERACTIONS
        int.remain <- sum(web) - sum(finalmat)
        if (int.remain > 0) {
            add <- sample(which(finalmat > 0), int.remain, prob = P1[finalmat > 
                0], replace = TRUE)
            finalmat[as.numeric(names(table(add)))] <- finalmat[as.numeric(names(table(add)))] + 
                table(add)
        }
        finalmat
    }
    # vectorised call of vaznull.fast:
    replicate(N, vaznull.fast(web), simplify = FALSE)
}
#vaznull(1, Safariland)
#system.time(vaznull(10, Safariland))
#system.time(swap.web(10, Safariland))
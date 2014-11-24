C.score <- function(web, normalise=TRUE, FUN=mean, ...){
    # calculates the C-score for all pollinator species; the C-score represents
    # the average number of checkerboard units for each unique specis pair.
    # (Stone & Roberts 1990; here taken from Gotelli & Rohde 2002)
    # these data are extremely skewed! A mean is, hm, not exactly what I would
    # suggest, hence the FUN-option allows for other summaries (try, e.g., hist).
    # option normalise ranges the index between 0 (no complementarity) and 1 (perfect
    # distinctness).
    # ... to be passed on to FUN (If a species occurs on all sites, then its distance
    # to any other will be 0, and maxD will also be 0. Since the division by 0 will
    # lead to NAs, use the ellipses to pass on na.rm to mean or similar functions. )
    # Carsten F. Dormann, Dec. 2007
    
    web <- web>0 # this whole concept works only on binary data!
    D <- designdist(t(web), method="(A-J)*(B-J)", terms="minimum")
    
    # The minimum value for Ds is 0, for the special case were all species use the
    # hosts exactly co-occurringly.
    # The maximum value for Ds in each comparison is AB, when they are exactly
    # complementary and hence J=0. However, if (A+B)>length of vector(L), then there
    # will be some co-occurrences and hence J>0=(A+B-L). The general maximum
    # then becomes (A-A-B+L)(B-A-B+L)=(L-B)(L-A). For (A+B)<L, maximum is AB.

    if (normalise){
#      L <- ncol(t(web))
      maxD <- designdist(t(web), method="ifelse(P<(A+B),(P-A)*(P-B), (A*B))", terms="minimum")
      nonzeros <- which(maxD != 0)                                
      if (length(nonzeros) > 0){
          use <- seq_along(D)[nonzeros] # if maxD is 0, then D/maxD is non-sense!
          D[use] <- D[use]/maxD[use]
          # return to a triangular format:
		  if (length(D) == length(nonzeros)){
		  	Dnew <- maxD
		  	Dnew[] <- D # [] needed to keep the lower-triangular format
		  	D <- Dnew
		  } 	          
      } else { 
      	if (any(maxD == 0) & any(D == 0)) {# for the rare case of so small networks that no checkerboard will fit in
      	  	D <- 0
      	  } else { D <- D/maxD }
      }
    }
    FUN(D, ...)
}
# example:
#m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,TRUE)
#C.score(m)


# corrected 2 Aug 2009: if maxD contained 0s, then C.score failed!
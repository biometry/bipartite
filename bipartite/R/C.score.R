C.score <- function(web, normalise=TRUE, FUN=mean, ...){
    # calculates the C-score for all pollinator species; the C-score represents
    # the average number of checkerboard units for each unique specis pair.
    # (Stone & Roberts 1990; here taken from Gotelli & Rohde 2002)
    # these data are extremely skewed! A mean is, hm, not exactly what I would
    # suggest, hence the FUN-option allows for other summaries (try, e.g., hist).
    # option normalise ranges the index between 0 (no complementarity) and 1 (perfect
    # distinctness).
    # ... to be passed on to FUN
    # Carsten F. Dormann, Dec. 2007

    web <- (web>0)*1 # this whole concept works only on binary data! (the "*1" converts output to numeric)
    D <- designdist(t(web), method="(A-J)*(B-J)", terms="binary")
    out <- FUN(D, ...)

    if (normalise){
	    simvec <- function(a, b, n){
	      # This function computes, independent of the actual data, how many checkerboards
	      # one would find, given degrees a and b of the two species involved, on n sites.
	      #
	      # The minimum value for Ds is 0, for the special case were all species use the
	      # hosts exactly co-occurringly.
	      # The maximum value for Ds in each comparison is AB, when they are exactly
	      # complementary and hence J=0. However, if (A+B)>length of vector(L), then there
	      # will be some co-occurrences and hence J>0=(A+B-L). The general maximum
	      # then becomes (A-A-B+L)(B-A-B+L)=(L-B)(L-A). For (A+B)<L, maximum is AB.
	      
	      res <- matrix(0, ncol=2, nrow=n)
			  if (!(a == 0 | b == 0 | a == n | b ==n)){
				  if (a <= n/2) res[1:(2*a),1] <- rep(c(1,0), times=a)
				  if (a > n/2){
					  put1here <- rep(c(1,0), times=floor(n/2))
					  res[1:length(put1here),1] <- put1here
					  remains <- a - sum(res[,1])
					  res[rev(which(res[,1]==0))[1:remains],1] <- 1
				  }
				  if (b <= n/2) res[1:(2*b), 2] <- rep(c(0,1), times=b)
				  if (b > n/2){
					  put1here <- rep(c(0,1), times=floor(n/2))
					  res[1:length(put1here), 2] <- put1here
					  remains <- b - sum(res[, 2])
					  res[rev(which(res[, 2]==0))[1:remains], 2] <- 1
				  }
			  }
			  #print(res)
			  as.numeric(designdist(t(res), method="(A-J)*(B-J)", terms="binary"))
		  }
		#tests:
		#simvec(8, 9, 10)
		#simvec(sum(Safariland[,1]>0), sum(Safariland[,1]>0), 9)


	  ## CLAUDE-based revision (checked by CFD, line by line 21.09.2026) 
	  # simvec depends only on the two species' occupancies and the number of sites,
		# and it is symmetric in a and b. The old code called it once per ordered pair of
		# species (including the diagonal) and then used only the lower triangle, i.e. it
		# recomputed the same few distinct values thousands of times on a wide web.
		# Here each distinct (min(a,b), max(a,b)) is computed once and looked up.
		nsites <- NROW(web)
		nspec <- NCOL(web)
		occupancy <- colSums(web > 0) # how many non-zero entries ("degrees") per higher
		pairs <- which(lower.tri(matrix(NA, nspec, nspec)), arr.ind=TRUE) # all possible combinations
		# ...which is exactly the storage order of the "dist" object D
		a <- occupancy[pairs[, 1]] # lists all targets of highers of the lower.tri
		b <- occupancy[pairs[, 2]] # lists all origins of highers of the lower.tri
		key <- paste(pmin(a, b), pmax(a, b), sep="_") # could use without pmin/pmax and a unique after the lookup; this is more elegant, admittedly
		uniq <- !duplicated(key)
		lookup <- vapply(which(uniq), function(i) simvec(a[i], b[i], nsites), numeric(1))
		names(lookup) <- key[uniq]
		Cmax <- unname(lookup[key])

		# Cmax is 0 whenever one of the two species occupies all or none of the sites: no
		# checkerboard can be formed, so the normalised C-score is undefined for that pair
		# (the observed D is 0 there too, so the old code returned 0/0 = NaN and a single
		# such species turned the whole index into NaN).
		undefined <- Cmax == 0
		if (any(undefined)){
			warning(sum(undefined), " of ", length(Cmax), " species pairs have no possible checkerboard ",
					"(a species occupying all or no sites); these pairs are excluded from the normalised C-score.",
					call.=FALSE)
		}
		out <- FUN(as.vector(D)[!undefined] / Cmax[!undefined], ...)
	}

     return(out)
}
# example:
#m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,TRUE)
#C.score(m)
#C.score(t(Safariland))


# corrected 2 Aug 2009: if maxD contained 0s, then C.score failed!

# added new way to compute maxD, because the old was wrong: 25.12.2014

# 21.9.2026: memoised the simvec lookup and restricted it to the lower triangle (was O(nspec^2)
# calls to designdist); undefined pairs are now dropped instead of yielding NaN. Claude + CFD

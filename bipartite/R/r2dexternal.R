r2dexternal <- function(N, web, abun.higher=NULL, abun.lower=NULL){
	## null model with external abundance estimates
	##
	## basis is r2dtable
	## marginal totals will be based on externally estimated abundances
	##
	## author: Carsten F. Dormann, 9. Sept. 2013
	
	if (!is.null(abun.higher) & length(abun.higher) != ncol(web)) stop("Number of species in the higher level is not the same as in the vector of abundances!")
	if (!is.null(abun.lower) & length(abun.lower) != nrow(web)) stop("Number of species in the lower level is not the same as in the vector of abundances!")	
	
	if (any(abun.higher == 0) | any(abun.lower == 0)) warning("Your external abundances include a 0. This will effectively delete this species from the analyses. I hope that is what you want!")
	
	if (is.null(abun.higher)) abun.higher <- colSums(web)
	if (is.null(abun.lower)) abun.lower <- rowSums(web)
	
	total <- sum(web)
		
	rel.abun.higher <- abun.higher/sum(abun.higher)
	rel.abun.lower <- abun.lower/sum(abun.lower)

	## columns: 
	cc <- floor(rel.abun.higher * total)
	missing.cc <- sum(cc) - total
	if (missing.cc != 0 & length(which(cc == 0)) > 0){
		# allocate to missing species first
		cc[which(cc == 0)] <- cc[which(cc == 0)] + 1
		missing.cc <- sum(cc) - total
	}

	if (missing.cc != 0 & length(which(cc == 0)) == 0){
		here.cc <- sample(1:ncol(web), abs(missing.cc), prob=rel.abun.higher, replace=TRUE)
		move.cc <- table(here.cc)
		if (missing.cc < 0){# some left over:
			cc[as.numeric(names(move.cc))] <- cc[as.numeric(names(move.cc))] + move.cc
		} else {		# some missing :
			cc[here.cc] <- cc[here.cc] - 1
		}
	} 
	
	## rows: 
	rr <- floor(rel.abun.lower * total)
	missing.rr <- sum(rr) - total
	if (missing.rr != 0 & length(which(rr == 0)) > 0){
		# allocate to missing species first
		rr[which(rr == 0)] <- rr[which(rr == 0)] + 1
		missing.rr <- sum(rr) - total
	}

	if (missing.rr != 0 & length(which(rr == 0)) == 0){
		here.rr <- sample(1:nrow(web), abs(missing.rr), prob=rel.abun.lower, replace=TRUE)
		move.rr <- table(here.rr)
		if (missing.rr < 0){# some left over:
			rr[as.numeric(names(move.rr))] <- rr[as.numeric(names(move.rr))] + move.rr
		} else {		# some missing :
			rr[here.rr] <- rr[here.rr] - 1
		}
	} 
	## Problem remaining: what if 1 is subtracted from a cell with only 1 in it?
	## should be a rare event, but may happen ...
	## then put in a while-loop when drawing "here.rr/here.cc", so that the 
	#  index only points to cells that have >1 interaction.
	
	# Next: investigate, with real data, how conservative r2dtable is: can we really ignore effects of one species on the other? If so, r2dtable and r2dexternal should yield very similar results.
	
	r2dtable(n=N, r=rr, c=cc)	
	
}

## test: 
# abun.lower <- c(15,5,2,7,4,8,6,0.01,6)
# set.seed(2)
# (abun.higher <- rpois(27, lambda=4))
# abun.higher[1] <- 0.001
# sum(ext.polls)
# web <- Safariland
# ## note: external abundances do not sum up; this is intentional!!

# r2dexternal(2, Safariland, abun.higher=abun.higher, abun.lower=abun.lower)
# r2dexternal(2, Safariland, abun.higher=abun.higher)
# r2dexternal(2, Safariland, abun.lower=abun.lower)


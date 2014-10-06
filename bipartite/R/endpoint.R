endpoint <- function(web){
	# following Gilarranz et al. 2011, page 1156, first paragraph.
	degH <- colSums(web>0)
	degL <- rowSums(web>0)
	out <- tcrossprod(degL, degH)
	out * (web>0)*1
}
#data(memmott1999)
#ends <- endpoint(memmott1999)
#weights.mean <- tapply(memmott1999, ends, mean)
#ends.weights <- tapply(ends, ends, mean)
#plot(weights.mean, ends.weights, log="xy", pch=16)
# bears little resemblance to Fig. 5 of Gilarranz et al. (2011)
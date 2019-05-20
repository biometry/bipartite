# metafunction to call betalinkr for all pairwise comparisons of a webarray with more than 2 webs
betalinkr_multi <- function(webarray, ...){
  data.out <- data.frame(i=integer(0), j=integer(0), S=numeric(0), OS=numeric(0), WN=numeric(0), ST=numeric(0))
  webnames <- dimnames(webarray)[[3]]
  if (is.null(webnames)) webnames <- 1:dim(webarray)[3]
  for (i in 1:(dim(webarray)[3]-1)){
    for (j in (i+1):dim(webarray)[3]){
      data.out[nrow(data.out)+1, ]  <- c(webnames[i], webnames[j], betalinkr(webarray=webarray[, , c(i,j)], ...))
    }
  }
  return(data.out)
}

# example
# betalinkr_multi(testarray)
# betalinkr_multi(webs2array(Safariland, vazquenc, vazarr), index="jaccard")
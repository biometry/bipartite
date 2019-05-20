# a helper function for betalinkr

# FUNCTION to turn an array with sites as third dimension into a siteXlink matrix (for dissim. calculations) ----------
  # names are optional (already taken care of in frame2webs webs2array)
array2linkmx <- function(webarray){
  # assumes the third dimension is the webID in the array, which will be the first dimension in the linkmx (vegan community matrix)
  # to easily fill the array into the matrix, first creating the transpose of linkmx
  linkmx.transp <- matrix(nrow = dim(webarray)[1] *  dim(webarray)[2], ncol = dim(webarray)[3])
    colnames(linkmx.transp) <- dimnames(webarray)[[3]] # use webIDs as names
  linkmx.transp[] <- webarray
    rownames(linkmx.transp) <- as.vector(outer(dimnames(webarray)[[1]], dimnames(webarray)[[2]], paste, sep="__")) # give names to links
  linkmx <- t(linkmx.transp) 
  return(linkmx)
}
# linkmx <- array2linkmx(webarray)  # example how to use
# array2linkmx(testarray)

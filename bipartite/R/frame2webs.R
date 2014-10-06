frame2webs <- function(dframe, varnames = c("lower", "higher", "webID", "freq"), type.out = "list", emptylist = TRUE) {
  # author: Jochen Fruend
  if (length(varnames)==4) {
    if (any(is.na(dframe[,varnames[4]]))) warning(paste("NAs in", varnames[4], "converted to 0"))
    webarray <- tapply(dframe[,varnames[4]],dframe[,varnames[1:3]], sum)
  }
  if (length(varnames)==3) webarray <- tapply(rep(1,nrow(dframe)),dframe[,varnames[1:3]], sum)
  webarray[is.na(webarray)] <- 0   # needs to be done when using tapply: unobserved combinations always get a zero, even with na.rm=T
  if (type.out=="array") return(webarray)
  if (type.out=="list") {
    weblist <- list()
    for (i in dimnames(webarray)[[3]]) weblist[[i]] <- webarray[,,i]
    if (emptylist) weblist <- lapply(weblist,empty)
    return(weblist)
  }
}

## example:
#testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
#frame2webs(testdata,type.out="array")
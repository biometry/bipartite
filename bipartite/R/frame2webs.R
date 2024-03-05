frame2webs <- function(dframe,
                       varnames = c("lower", "higher", "webID", "freq"),
                       type.out = "list", emptylist = TRUE) {
  # author: Jochen Fruend
  if (length(varnames) == 4) {
    if (any(is.na(dframe[, varnames[4]])))
      warning(paste("NAs in", varnames[4], "converted to 0"))
    webarray <- tapply(dframe[, varnames[4]], dframe[, varnames[1:3]], sum)
  }
  if (length(varnames) == 3)
    webarray <- tapply(rep(1, nrow(dframe)), dframe[, varnames[1:3]], sum)
  # needs to be done when using tapply: unobserved combinations always get a zero, even with na.rm=T
  webarray[is.na(webarray)] <- 0
  if (type.out == "array")
    return(webarray)
  if (type.out == "list") {
    weblist <- list()
    for (i in dimnames(webarray)[[3]]) {
      weblist[[i]] <- webarray[, , i]
    }
    if (emptylist)
      weblist <- lapply(weblist, empty)
    return(weblist)
  }
}

## example:
#testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
#frame2webs(testdata,type.out="array")

frame2webs_fast <- function(dframe, 
                            varnames = c("lower", "higher", "webID", "freq"), 
                            type.out = "list", emptylist = TRUE) {
  if (length(varnames) == 4){
    INDEX <- lapply(dframe[, varnames[1:3]], as.factor)
    nI <- 3
    namelist <- lapply(INDEX, levels)
    extent <- lengths(namelist, use.names = FALSE)
    cumextent <- cumprod(extent)
    if (cumextent[nI] > .Machine$integer.max)
      stop("total number of levels >= 2^31")
    storage.mode(cumextent) <- "integer"
    ngroup <- cumextent[nI]
    group <- as.integer(INDEX[[1L]])
    if (nI > 1L){
      for (i in 2L:nI){
        group <- group + cumextent[i - 1L] * (as.integer(INDEX[[i]]) - 1L)
      }
    }
    # levels(group) <- as.character(seq_len(ngroup))
    # class(group) <- "factor"
    # print(group)
    # print(namelist)
    ans <- split(dframe[, varnames[4]], group, drop = TRUE)
    ans <- lapply(ans, sum)
    ansmat <- array(0, dim=extent, dimnames = namelist)
    # print(ans)
    # print(as.integer(names(ans)))
    # print(unlist(ans, recursive = FALSE, use.names = FALSE))
    ansmat[as.integer(names(ans))] <- unlist(ans, recursive = FALSE, use.names = FALSE)
    # print(ansmat)
    return(ansmat)
  }
}
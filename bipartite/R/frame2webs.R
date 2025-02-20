###
# @author: Tobias Bauer, Jochen Fruend
# Convenience function to convert a table of observations
# (i.e. an "edge list", typically compiled in a spreadsheet programm)
# into a network matrix for further use in bipartite.
#
# Optimization with partial reimplementation of tapply function by Tobias Bauer
frame2webs <- function(dframe,
                       varnames = c("lower", "higher", "webID", "freq"),
                       type.out = "list", emptylist = TRUE) {
  # If only 3 varnames are defined create a vector filled with 1s
  if (length(varnames) == 4) {
    values <- dframe[, varnames[4]]
  } else if (length(varnames) == 3) {
    values <- rep(1, nrow(dframe))
  }
  # !!! THE FOLLOWING CODE IS ALMOST AN EXACT COPY OF
  # THE STANDARD R FUNCTION TAPPLY WITH MINOR
  # CHANGES FOR OPTIMIZATION !!!
  # Convert all three the indices to factors
  INDEX <- lapply(dframe[, varnames[1:3]], as.factor)
  # Number of indices
  nI <- 3
  # Get the levels of all three indices
  namelist <- lapply(INDEX, levels)
  # Get the number of levels for all three indices
  extent <- lengths(namelist, use.names = FALSE)
  # Calculate the cumulative product of the index lengths
  # i.e n(index 1) * n(index 2) * n(index 3)
  # This value is equal to the total size of the output array.
  cumextent <- cumprod(extent)
  # Throw an error if the total size is larger than the maximum integer
  # size of the machine.
  if (cumextent[nI] > .Machine$integer.max)
    stop("total number of levels >= 2^31")
  storage.mode(cumextent) <- "integer"
  # > ngroup <- cumextent[nI]
  # Convert the index combination to unique integers
  group <- as.integer(INDEX[[1L]])
  if (nI > 1L) {
    for (i in 2L:nI) {
      group <- group + cumextent[i - 1L] * (as.integer(INDEX[[i]]) - 1L)
    }
  }
  # In the standard tapply function at this point a factor of all integers
  # up to ngroup is generated. In my opinion this is unnecessary and consumes
  # way to much memory.
  # > levels(group) <- as.character(seq_len(ngroup))
  # > class(group) <- "factor"
  # > print(group)
  # > print(namelist)
  # Split the observation data only regarding index combination,
  # that are actually present in the data.frame (drop = TRUE).
  # That way the calculation time for sparse matrices is reduced immensely.
  ans <- split(values, group, drop = TRUE)
  # Now apply the sum function on the sparse split list.
  ans <- lapply(ans, sum)
  # Create the resulting array with the full extent of all possible
  # index combination. Non-existing combinations are filled
  # with the default value 0. The dimension names are taken from the
  # factor levels of the indices.
  ansmat <- array(0, dim = extent, dimnames = namelist)
  # Now fill the array at the location for which index combination
  # existed in the data.frame with the calculated sum-values.
  ansmat[as.integer(names(ans))] <- unlist(ans, recursive = FALSE,
                                           use.names = FALSE)
  if (type.out == "array")
    return(ansmat)
  # Convert array to list and throw away empty rows or columns if emptylist == T
  if (type.out == "list") {
    weblist <- list()
    for (i in dimnames(ansmat)[[3]]) {
      # Because of the nonsensical behavior in R, that either all dimensions
      # of length 1 or no dimensions at all are dropped we need this weird
      # workaround...
      weblist[[i]] <- array(ansmat[, , i, drop = FALSE],
                            dim = extent[-3],
                            dimnames = namelist[-3])
    }
    if (emptylist)
      weblist <- lapply(weblist, empty)
    return(weblist)
  }
}

## example:
#testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
#frame2webs(testdata,type.out="array")

metaComputeModules <- function(moduleObject, N=5, ...){
    attempts <- list()
    for (i in 1:N){
        attempts[[i]] <- computeModules(moduleObject, ...)
    }
    likelik <- sapply(attempts, function(x) x@likelihood)
    return(attempts[[which.max(likelik)]])
}
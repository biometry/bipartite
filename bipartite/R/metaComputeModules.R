metaComputeModules <- function(moduleObject, N=5, method="Beckett", ...){
    attempts <- list()
    for (i in 1:N){
        attempts[[i]] <- computeModules(moduleObject, method=method, ...)
    }
    likelik <- sapply(attempts, function(x) x@likelihood)
    return(attempts[[which.max(likelik)]])
}
`clustering_tm` <-
function(net,subsample=1,seed=NULL){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet)) {
    if(ncol(net)==3) {
      net <- as.tnet(net, type="weighted two-mode tnet")
    } else {
      net <- as.tnet(net, type="binary two-mode tnet")
    }
  }
  if(attributes(net)$tnet!="binary two-mode tnet" & attributes(net)$tnet!="weighted two-mode tnet")
    stop("Network not loaded properly")
  if(!is.null(seed))
    set.seed(as.integer(seed))

  # 1-paths i1+p1
  paths <- net
  if(subsample!=1) {
    if(subsample<1) {
      index <- sample.int(nrow(paths), round(nrow(paths)*subsample))
    } else {
      index <- sample.int(nrow(paths), as.integer(subsample))
    }
    index <- index[order(index)]
    paths <- paths[index,]
  }
  if(attributes(net)$tnet=="binary two-mode tnet")  {
    dimnames(paths)[[2]] <- c("i1","p1")
    # 2-paths i1+p1+i2
    dimnames(net)[[2]] <- c("i2","p1")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i2"],]
    # 3-paths i1+p1+i2+p2
    dimnames(net)[[2]] <- c("i2","p2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"p1"] != paths[,"p2"],]
    # 4-paths i1+p1+i2+p2+i3
    dimnames(net)[[2]] <- c("i3","p2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i3"],]
    paths <- paths[paths[,"i2"] != paths[,"i3"],]
    denominator <- nrow(paths)
    # Find which 4-paths are part of 6-cycles
    dimnames(net)[[2]] <- c("i","p")
    paths <- paths[order(paths[,"i1"], paths[,"p1"], paths[,"p2"], paths[,"i3"]),c("i1","p1","p2","i3")]
    net.list <- split(net[,"p"], net[,"i"])
    numerator <- sum(apply(paths, 1, function(a) {ct <- c(net.list[[as.character(a[1])]],net.list[[as.character(a[4])]]); return(sum(duplicated(ct[ct!=a[2] & ct!=a[3]]))>0)}))
  } else {
    dimnames(paths)[[2]] <- c("i1","p1","w1")
    # 2-paths i1+p1+i2
    dimnames(net)[[2]] <- c("i2","p1","w2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i2"],]
    # 3-paths i1+p1+i2+p2
    dimnames(net)[[2]] <- c("i2","p2","w3")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"p1"] != paths[,"p2"],]
    # 4-paths i1+p1+i2+p2+i3
    dimnames(net)[[2]] <- c("i3","p2","w4")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i3"],]
    paths <- paths[paths[,"i2"] != paths[,"i3"],]
    pw <- data.frame(
      bi=rep(1, nrow(paths)),
      am=rowMeans(paths[,c("w1","w2","w3","w4")]), 
      gm=apply(paths[,c("w1","w2","w3","w4")], 1, function(a) sqrt(sqrt(prod(a)))), 
      ma=apply(paths[,c("w1","w2","w3","w4")], 1, max), 
      mi=apply(paths[,c("w1","w2","w3","w4")], 1, min))
    denominator <- colSums(pw)
    # Find which 4-paths are part of 6-cycles
    dimnames(net)[[2]] <- c("i","p","w")
    net <- net[,c("i","p")]
    paths <- paths[,c("i1","p1","p2","i3")]
    net.list <- split(net[,"p"], net[,"i"])
    index <- apply(paths, 1, function(a) {ct <- c(net.list[[as.character(a[1])]],net.list[[as.character(a[4])]]);sum(duplicated(ct[ct!=a[2] & ct!=a[3]]))>0})
    numerator <- colSums(pw[index,])
  }
  rm(paths)
  # Fraction
  return(numerator/denominator)
}

`symmetrise_w` <-
function(net,method="MAX"){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")

  # Join the net with it's reversed version
  net <- rbind(net, cbind(i=net[,"j"], j=net[,"i"], w=0))
  # Remove exact duplicates (i=i & j=j)
  net <- net[!duplicated(net[,c("i","j")]),]
  # Change ties so that i<j
  net[net[,"i"]>net[,"j"],c("i","j")] <- net[net[,"i"]>net[,"j"],c("j","i")];
  # Order ties (the greatest weight first)
  net <- net[order(net[,"i"],net[,"j"], -net[,"w"]),]
  # Create an index
  dup <- cumsum(rep.int(c(1,0), nrow(net)/2))
  # Create a weight vector
  w <- switch(method,
    MAX   = net[rep(c(TRUE,FALSE), length=nrow(net)),"w"],
    MIN   = net[rep(c(FALSE,TRUE), length=nrow(net)),"w"],
    AMEAN = tapply(net[,"w"], dup, mean),
    SUM   = tapply(net[,"w"], dup, sum),
    GMEAN = tapply(net[,"w"], dup, function(a) sqrt(a[1]*a[2])),
    PROD  = tapply(net[,"w"], dup, function(a) a[1]*a[2]),
    DIFF  = tapply(net[,"w"], dup, function(a) abs(a[1]-a[2])))
  # Extract only one entry per undirected tie
  net <- net[rep(c(TRUE,FALSE), length=nrow(net)),]
  # Add the weight vector to this list
  net[,"w"] <- w;
  # Only keep ties with a positive weight
  net <- net[net[,3]>0,]
  # Join this net with its reversed version
  net <- rbind(cbind(net[,1],net[,2],net[,3]), cbind(net[,2],net[,1],net[,3]))
  # Assign names to columns
  dimnames(net)[[2]]<-c("i","j","w")
  # Order net
  net <- net[order(net[,"i"],net[,"j"]),]
  row.names(net)<-NULL
  return(net)
}

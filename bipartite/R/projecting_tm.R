`projecting_tm` <-
function(net,method="sum"){
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
  net2 <- net
  # Add a weight column of 1 if the network is binary
  if(attributes(net)$tnet=="binary two-mode tnet")
    net2 <- cbind(net2, w=1)
  # Define support objects
  net2 <- net2[order(net2[,"i"], net2[,"p"]),]
  np <- table(net2[,"p"])
  net2 <- merge(net2, cbind(p=as.numeric(row.names(np)), np=np))
  # Create the one-mode topology
  net1 <- merge(net2, cbind(j=net2[,"i"],p=net2[,"p"]))
  net1 <- net1[net1[,"i"]!=net1[,"j"],c("i","j","w","np")]
  net1 <- net1[order(net1[,"i"],net1[,"j"]),]
  # Calculate weights in the one-mode network
  index <- !duplicated(net1[,c("i","j")])  
  w <- switch(method,
      binary = rep(1, sum(index)),
      sum    = tapply(net1[,"w"], cumsum(index), sum),
      Newman = tapply(1:nrow(net1), cumsum(index), function(a) sum(net1[a,"w"]/(net1[a,"np"]-1))))
  net1 <- cbind(net1[index,c("i","j")],w=as.numeric(w))
  return(as.tnet(net1, type="weighted one-mode tnet"))
}

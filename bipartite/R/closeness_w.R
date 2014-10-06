`closeness_w` <-
function(net, directed = NULL, gconly = TRUE, precomp.dist = NULL, alpha=1){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")
  
  # Transform weights accordingly to alpha (default 1 equal no change)
  net[,"w"] <- net[,"w"]^alpha

  # Compute distance matrix                
  if(is.null(precomp.dist)) {
    # Check if network is directed
    if(is.null(directed)) {
      tmp <- symmetrise_w(net, method = "MAX")
      directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
    }
    precomp.dist <- distance_w(net = net, directed = directed, gconly = gconly)
  }
  # Change algorithm based on gconly parameter
  if(!gconly) {
    precomp.dist <- 1/precomp.dist
  }
  # Sum up distances to all other nodes to get farness
  out <- cbind(
    node = attributes(precomp.dist)$nodes, 
    closeness = rowSums(precomp.dist, na.rm = TRUE),
    n.closeness=NaN)
  # If only gconly
  if(gconly) {
    out[, "closeness"] <- 1/out[, "closeness"]
  }
  # Normalise scores by N-1 (not always appropriate for weighted networks!)
  out[,"n.closeness"] <- out[,"closeness"]/(nrow(out)-1)
  # Return object
  return(out)
}


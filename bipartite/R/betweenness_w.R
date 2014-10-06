`betweenness_w` <-
function(net, directed = NULL, alpha=1){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))
    net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")
    stop("Network not loaded properly")

  # Check if network is directed
  if(is.null(directed)) {
    tmp <- symmetrise_w(net, method = "MAX")
    directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
  }
  
  # Transform weights -- key when high weights are strong ties!!
  net[,"w"] <- (1/net[,"w"])^alpha

  # Load in igraph
  g <- tnet_igraph(net, type="weighted one-mode tnet", directed=directed)

  # Create output object
  N <- length(V(g))
  out <- cbind(node = 1:N, betweenness = 0)
  out[,"betweenness"] <- igraph::betweenness(graph=g, directed = directed)
  return(out)                       
}


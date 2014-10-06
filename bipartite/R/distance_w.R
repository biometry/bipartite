`distance_w` <-
function(net, directed = NULL, gconly = TRUE, subsample = 1, seed = NULL){
  # Check whether confirms to tnet standard
  if (is.null(attributes(net)$tnet)) 
    net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet") 
    stop("Network not loaded properly")
  
  # Set seed
  if (!is.null(seed)) 
    set.seed(as.integer(seed))
    
  # Normalise weightes
  net[, "w"] <- mean(net[, "w"])/net[, "w"]
  
  # Check if network is directed
  if(is.null(directed)) {
    tmp <- symmetrise_w(net, method = "MAX")
    directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
  }
  
  # Do computation in igraph
  g <- tnet_igraph(net, type="weighted one-mode tnet", directed=directed)
  if(gconly) {
    gc <- igraph::clusters(g, mode="strong")
    gc <- which(gc$membership==names(sort(-table(gc$membership)))[1])
    g <- igraph::induced.subgraph(graph=g, vids=gc, "create_from_scratch")
  } else {
    gc <- as.integer(V(g))
  }
  d <- igraph::shortest.paths(g, mode="out", weights=igraph::get.edge.attribute(g, "tnetw"))
  diag(d) <- NA
  attributes(d)$nodes <- gc
  return(d)
}

fc <- function (web, dist="euclidean", method="average", weighted=TRUE) {
 # function to calculate the functional diversity for the rows of a web
 # follwing Devoto M., Bailey S., Craze P. & Memmott J. (2012) 
 # Understanding and planning ecological restoration of plant-pollinator networks
 # Ecology Letters, http://dx.doi.org/10.1111/j.1461-0248.2012.01740.x.
  web <- empty(web)
  if (!weighted) web <- (web>0)*1
  if (nrow(web)<2) {
    #warning ("nrows<2; FD for this web is 0")
    functdiv <- 0
  } else {
  dist.matrix <- dist(web, method=dist) #calculate dissimilarity matrix
  ddgram <- hclust(dist.matrix, method=method) #calculate dendrogram
  functdiv <- treeheight(ddgram) #calculate functional diversity
  }
  return (functdiv)
}

#data(Safariland)
#fc(Safariland)    
#fc(t(Safariland), dist="canberra", method="complete", weighted=FALSE)
    
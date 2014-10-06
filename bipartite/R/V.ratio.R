V.ratio <- function(web){
  # invented by Schluter (1984)
  # be careful with the interpretation!!!
  web <- web>0

  N <- ncol(web) #number of species
  M <- nrow(web) #number of "samples"
  Tj <- colSums(web) #species "density" per sample j
  ni <- rowSums(web) #nr. of individuals per species i
  p <- ni/N

  sigma.i.2 <- p*(1-p)
  S.T.2 <- 1/N*sum((Tj-mean(Tj))^2)
  V <- S.T.2 / sum(sigma.i.2)
  V
}
#example:
#data(Safariland)
#V.ratio(Safariland)



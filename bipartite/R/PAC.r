PAC <- function(web){
  # calculates the potential for apparent competition, following Morris et al. (2005)
  # web     a host parasitoid network, where the entries represent the sum of parasitoids emerging from the interactions between parasitoid and host (i.e. number of interactions * number of parasitoid individuals emerging from each host) 
  #
  # d_ij = sum_k( alpha_ik/sum_l(alpha_il) * alpha_jk/sum_m(alpha_mk))
  #
  #d_ij = effect of species j on i (parasitoids developed on j will attack i)
  #k    = parasitoid
  #i, j = host
  #alpha = "link strength" (i.e. entry in the network matrix)
  #
  # result is a m x m matrix, where m is the number of lower trophic level species; the diagonal is the intraspecific PAC; the upper triangle is the effect of j on i (d_ij), the lower triangle that of i on j (d_ji).
  # Reference: 
  # Mueller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370.
  # Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2005. Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449--462.

    m <- nrow(web)
    res <- matrix(NA, ncol=m, nrow=m)
    colnames(res) <- rownames(res) <- rownames(web)
 
    for (i in 1:nrow(web)){
        for (j in 1:nrow(web)){
            dd <- 1 # initialising a container to be summed over per i-j-pair
            for (k in 1:ncol(web)){
                alpha_ik <- web[i,k]
                alpha_jk <- web[j,k]
                sum.alpha.il <- sum(web[i,])
                sum.alpha.mk <- sum(web[,k])
                dd[k] <- alpha_ik/sum.alpha.il * alpha_jk/sum.alpha.mk
            } 
            res[i, j] <- sum(dd, na.rm=TRUE)   
        }
    }
    res
  
}
#data(Safariland)
#PAC(Safariland)
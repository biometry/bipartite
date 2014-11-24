compart <- function(web) {
# Author: Juan M. Barreneche
# E-mail: jumanbar@gmail.com
# replacing the slower version below 22 May 2012
  cweb <- as.matrix(web)
  rc <- which(web > 0, arr.ind=TRUE)
  r <- rc[,1]
  c <- rc[,2]
  p1 <- numeric(length(r)) + 1
  p2 <- p1
  comp <- 1
  while (any(p2 > 0)) {
    p2[which(p2 == 1)[1]] <- -comp
    while (any(p2 != p1)) {
      p1 <- p2
      p2[r %in% r[p2 == -comp]] <- -comp
      p2[c %in% c[p2 == -comp]] <- -comp
    }
    locs <- matrix(rc[p2 == -comp,], ncol=2)
    cweb[locs] <- -comp
    comp <- comp + 1
  }
  out <- list(cweb=cweb, n.compart=max(abs(p2)))
  out
}


# This (faster) version is affected by the occurrence of ties!
#compart <- function(web){
#    # finds compartments and assigns numbers to them:
#    # 1. organise web along the diagonal:
#    if (any(dim(web))<2){
#      out <- list(cweb=web, n.compart=1)
#    } else {
#      ca <- cca(web)
#      rindex <- order(summary(ca)$sites[,1], decreasing=TRUE)
#      cindex <- order(summary(ca)$species[,1], decreasing=TRUE)
#      web2 <- web[rindex, cindex]
#      # 2. analyse for compartments, first rows, then cols:
#      cr <- distconnected(vegdist(web2, "jaccard"), trace=FALSE)
#      cc <- distconnected(vegdist(t(web2), "jaccard"), trace=FALSE)
#      # 3. numbered sequences are now matching each other to describe the clusters:
#      clu.cr <- matrix(rep(cr, ncol(web2)), nrow=nrow(web2))
#      clu.cc <- matrix(rep(cc, nrow(web2)), ncol=ncol(web2), byrow=TRUE)
#      clu <- web
#      clu[clu.cr == clu.cc] <- clu.cr[clu.cr == clu.cc]
#      # 4. resort the web:
#      outweb <- clu[order(rindex), order(cindex)]
#      out <- list(cweb=outweb, n.compart=max(cr))
#    }
#    return(out)
#}
#
## example:
## make a nicely compartet web:
#web <- matrix(0, 10,10)
#web[1,1:3] <- 1 
#web[2,4:5] <- 1 
#web[3:7, 6:8] <- 1
#web[8:10, 9:10] <- 1
#web <- web[-c(4:5),] #oh, and make it asymmetric!
#web <- web[,c(1:5, 9,10, 6:8)] #oh, and make it non-diagonal
#compart(web)
#
## or, standard, use Safariland as example:
#data(Safariland)
#compart(Safariland)
#

# ## Bernd's version (not affected by lots of similarly linked species):
# compart <- function(web){
  # #recursive crossing out of horizontal and vertical neighbours
  # # can be simplified using vegan's distconnected on cca output ...
  # cross = function(web, start, comp)
  # {
    # n.r=nrow(web)
    # n.c=ncol(web)
    # r=start[1]
    # c=start[2]
    # web[r,c]=-comp     #assign a negative compartment number to interaction

    # for (i in 1:n.r) #schaue senkrecht
    # {
    # if (web[i,c]>0)  web<-cross(web,c(i,c),comp)
    # }

    # for (i in 1:n.c)    #schaue waagrecht
    # {
     # if (web[r,i]>0) web<-cross(web,c(r,i),comp)
    # }
    # return (web)
  # }

  # comp=1     #start with the first compartment
  # while (max(web)>0)  #any interactions left?
  # {
    # start=which(web==max(web),arr.ind=TRUE)[1,]  #start at the highest number of interactions (arbitrary)
    # web <- cross(web,start,comp) #start recursion until no more neighbours in one compartment are found
    # comp <- comp+1   #go to the next compartment
  # }
  # return(list(cweb=web, n.compart=max(abs(web))))
# }

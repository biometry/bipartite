library(bipartite)
library(tnet)

# gplot(as.one.mode(Safariland, project="higher")) # not very helpful ...
plotweb(Safariland, method="normal")

one.mode.tnet <- as.tnet(t(Safariland)) # to look at higher level
one.mode.bip <- web2edges(Safariland, return=T) # different labelling system

p.one.mode <- projecting_tm(one.mode.tnet, method="binary")
projecting_tm(one.mode.bip) # identical

igraph::betweenness(tnet_igraph(p.one.mode))

betweenness_w(p.one.mode)
sna::betweenness(as.one.mode(Safariland, project="higher"), cmode="undirected")
sna::betweenness(Safariland, cmode="undirected")

!! sequence of species differs: tnet by rows, sna by columns!

web <- Safariland 

one.mode.tnet <- as.tnet(web)
p.one.mode <- projecting_tm(one.mode.tnet, method="sum")
betweenness_w(p.one.mode)
sna::betweenness(as.one.mode(web, project="lower"), cmode="undirected")
igraph::betweenness(tnet_igraph(p.one.mode))

# 3 ways, three values ...

gplot(as.one.mode(Safariland, project="lower"), label=1:9)
specieslevel(web, index="betweenness")

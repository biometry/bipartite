library(bipartite)
library(tnet)

gplot(as.one.mode(Safariland, project="lower"), label=rownames(Safariland))# not very helpful ...
plotweb(Safariland, method="normal")

one.mode.tnet <- as.tnet(Safariland) # to look at higher level
one.mode.bip <- web2edges(Safariland, return=T) # different labelling system; higher trophic level starts at nrow+1

p.one.mode <- projecting_tm(one.mode.tnet, method="sum")
projecting_tm(one.mode.bip) # identical, if method = "sum"

# all the same, i.e. no effect of "type":
igraph::betweenness(tnet_igraph(p.one.mode))
igraph::betweenness(tnet_igraph(p.one.mode, type="weighted two-mode tnet"))
igraph::betweenness(tnet_igraph(p.one.mode, type="weighted one-mode tnet"))
igraph::betweenness(tnet_igraph(p.one.mode, type="binary two-mode tnet"))

# betweenness:
bipartite::betweenness_w(p.one.mode) # length 9
bipartite::BC(Safariland, rescale=F)$lower
# RBGL::betweenness.centrality.clustering() # format of graph unclear
DiagrammeR::get_betweenness(DiagrammeR::from_igraph(tnet_igraph(p.one.mode)))
igraph::estimate_betweenness(tnet_igraph(p.one.mode), cutoff=0)
# igraph::cluster_edge_betweenness(tnet_igraph(p.one.mode)) # something else, really
influenceR::betweenness(tnet_igraph(p.one.mode))
sna::betweenness(Safariland) # length 36
sna::betweenness(as.matrix(p.one.mode)) # length 27!!
sna::betweenness(one.mode.tnet) # same as before, by a factor




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

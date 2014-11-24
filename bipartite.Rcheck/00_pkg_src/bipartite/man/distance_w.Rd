\name{distance_w}
\alias{distance_w}
\title{ Distance in a weighted network }
\description{ The shortest path length, or geodesic distance, between two nodes in a binary network is the minimum number of steps you need to make to go from one of them to the other. This distance is the quickest connection between nodes when all ties are the same. However, in a weighted network, all ties are not the same. See http://toreopsahl.com/2009/01/09/average-shortest-distance-in-weighted-networks/ for more deatails. }
\usage{ distance_w(net, directed=NULL, gconly=TRUE, subsample=1, seed=NULL) }
\arguments{
  \item{net}{A weighted edgelist}
  \item{directed}{ logical, whether the network is directed or undirected. Default is NULL, this means that the function checks whether the edgelist is directed or not.}
  \item{gconly}{ logical, whether the function should only be calculated for the giant component. Default is TRUE.}
  \item{subsample}{Whether a only a subset of starting nodes should we used when calculating the measure. This is particularly useful when running out of memory analysing large networks. If it is set to 1, all distances are analysed. If it set to a value below one, this is roughly the proportion of starting noes that will be analysed. If it is set to an interger greater than 1, this number of starting nodes that will be analysed.}
  \item{seed}{If a subset of starting nodes is analysed, by setting this parameter, the results are reproducable.}
}
\value{ Returns a distance matrix. }
\references{ http://toreopsahl.com/2009/01/09/average-shortest-distance-in-weighted-networks/ }
\author{ Tore Opsahl; http://toreopsahl.com }
\note{ version 1.0.0, taken, with permission, from package tnet }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
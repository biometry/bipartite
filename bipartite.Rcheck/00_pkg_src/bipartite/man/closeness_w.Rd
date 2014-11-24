\name{closeness_w}
\alias{closeness_w}
\title{ Closeness centrality in a weighted network }
\description{ This function calculates closeness scores for nodes in a weighted network based on the distance_w-function. }
\usage{ closeness_w(net, directed=NULL, gconly=TRUE, precomp.dist=NULL, alpha=1) }

\arguments{
  \item{net}{ A weighted edgelist }
  \item{directed}{ Logical: whether the edgelist is directed or undirected. Default is NULL, then the function detects this parameter. }
  \item{gconly}{ Logical: whether to calculate closeness only on the main component (traditional closeness). Default is TRUE. If this parameter is set to FALSE, a closeness measure for all nodes is computed. For details, see \url{http://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/}}
  \item{precomp.dist}{ If you have already computed the distance matrix using distance_w-function, you can enter the name of the matrix-object here. }
  \item{alpha}{ sets the alpha parameter in the generalised measures from Opsahl, T., Agneessens, F., Skvoretz, J., 2010. Node Centrality in Weighted Networks: Generalizing Degree and Shortest Paths. Social Networks. If this parameter is set to 1 (default), the Dijkstra shortest paths are used. The identification procedure of these paths rely simply on the tie weights and disregards the number of nodes on the paths. }
}
\value{ Returns a data.frame with three columns: the first column contains the nodes' ids, the second column contains the closeness scores, and the third column contains the normalised closeness scores (i.e., divided by N-1). }
\references{ http://toreopsahl.com/2009/01/09/average-shortest-distance-in-weighted-networks/ }
\author{ Tore Opsahl; http://toreopsahl.com }
\note{ version 1.0.0, taken, with permission, from package tnet  }
\examples{
## Load sample data
sampledata <- rbind(
c(1,2,4),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(4,2,1),
c(5,2,2),
c(5,6,1),
c(6,5,1))

## Run the programme
closeness_w(sampledata)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
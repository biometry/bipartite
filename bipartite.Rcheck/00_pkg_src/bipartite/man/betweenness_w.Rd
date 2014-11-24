\name{betweenness_w}
\alias{betweenness_w}
\title{ Betweenness centrality in a weighted network }
\description{ This function calculates betweenness scores for nodes in a weighted network based on the distance_w-function. \cr
Note: This algorithm relies on the igraphs package's implementation of Dijkstra's algorithm. Currently, it does not find multiple shortest paths if two exist.}
\usage{ betweenness_w(net, directed=NULL, alpha=1) }

\arguments{
  \item{net}{ A weighted edgelist }
  \item{directed}{ logical, whether the network is directed or undirected. Default is NULL, this means that the function checks whether the edgelist is directed or not.}
  \item{alpha}{ sets the alpha parameter in the generalised measures from Opsahl, T., Agneessens, F., Skvoretz, J., 2010. Node Centrality in Weighted Networks: Generalizing Degree and Shortest Paths. Social Networks. If this parameter is set to 1 (default), the Dijkstra shortest paths are used. The length of these paths rely simply on the tie weights and disregards the number of nodes on the paths. }
}
\value{ Returns a data.frame with two columns: the first column contains the nodes' ids, and the second column contains the nodes' betweenness scores. }
\references{ http://toreopsahl.com/2009/02/20/betweenness-in-weighted-networks/ }
\author{ Tore Opsahl; http://toreopsahl.com }
\note{ version 1.0.0, taken, with permission, from package tnet  }
\examples{
## Load sample data
sampledata <- rbind(
c(1,2,1),
c(1,3,5),
c(2,1,1),
c(2,4,6),
c(3,1,5),
c(3,4,10),
c(4,2,6),
c(4,3,10))

## Run the programme
betweenness_w(sampledata)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
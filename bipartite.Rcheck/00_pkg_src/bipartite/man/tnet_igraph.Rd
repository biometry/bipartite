\name{tnet_igraph}
\alias{tnet_igraph}
\title{ Exports a tnet network to an igraph object }
\description{
  The tnet_igraph function creates an igraph object from a tnet network.
}
\usage{tnet_igraph(net,type=NULL, directed=NULL)}
\arguments{
  \item{net}{A tnet network }
  \item{type}{ type of tnet network, see as.tnet. }
  \item{directed}{ if a one-mode networks, this can be set to avoid testing whether the network is directed. }
}
\value{Returns the igraph object.}
\references{ http://toreopsahl.com/ }
\author{ Tore Opsahl; http://toreopsahl.com }
\note{ version 1.0.0, taken, with permission, from package tnet  }
\examples{
## Load sample data
sample <- rbind(
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
tnet_igraph(sample, type="weighted one-mode tnet")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
\name{symmetrise_w}
\alias{symmetrise_w}
\title{ Symmetrise_w }
\description{ The symmetrise_w-function creates an undirected one-mode network from a directed one-mode network. }
\usage{symmetrise_w(net, method="MAX")}
\arguments{
  \item{net}{A one-mode network }
  \item{method}{ the method used to decide the weight of the undirected edge. It can be: 
  "MAX" sets the weight to the maximum of the weight(s) of the arc(s)
  "MIN" sets the weight to the minimumof the weight(s) of the arc(s)
  "AMEAN" sets the weight to the average (arithmetic mean) of the weight(s) of the arc(s)
  "SUM" sets the weight to the sum of the weight(s) of the arc(s)
  "PROD" sets the weight to the product of the weight(s) of the arc(s)
  "DIFF" sets the weight to the absolute difference between the weight(s) of the arc(s) }
}
\value{Returns the undirected network }
\references{ \url{https://toreopsahl.com/2008/11/28/network-weighted-network/} }
\author{ Tore Opsahl; \url{https://toreopsahl.com/} }
\note{ version 1.0.0, taken, with permission, from package tnet }
\examples{
## Load sample data
sample <- rbind(
c(1,2,2),
c(1,3,2),
c(2,1,4),
c(2,3,4),
c(2,4,1),
c(2,5,2),
c(3,1,2),
c(3,2,4),
c(5,2,2),
c(5,6,1))

## Run the programme
symmetrise_w(sample, method="MAX")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
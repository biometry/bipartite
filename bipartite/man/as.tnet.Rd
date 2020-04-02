\name{as.tnet}
\alias{as.tnet}
\title{ Ensures that networks conform to the tnet stardards }
\description{
  Checks that a network conforms to the tnet stardards, and attaches a label. If the type parameter is not set, the network is assumed to be a binary two-mode network, a weighted one-mode network, or a longitudinal network if there are 2, 3, or 4 columns respectively. Moreover, if a matrix is entered (more than 4 columns and rows), it is assumed to be a weighted one-mode network if square or a two-mode network if non-square.
}
\usage{as.tnet(net, type=NULL)}
\arguments{
  \item{net}{A network in an edgelist or matrix format. It can be a weighted one-mode network, a binary two-mode network, a weighted two-mode netork, or a longitudinal network. If the data-object has two-columns, it is assumed to be a binary two-mode network; three columns, weighted one-mode network; four columns, longitudinal; five or more and the same number of rows and columns, weighted one-mode network; five or more and --not-- the same number of rows and columns, it is assumed to be a two-mode network. }
  \item{type}{If you would like to specify the type of network. This could be "weighted one-mode tnet", "binary two-mode tnet", "weighted two-mode tnet", or "longitudinal tnet". }
}
\value{Returns the network with an attached label.}
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
as.tnet(sample)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
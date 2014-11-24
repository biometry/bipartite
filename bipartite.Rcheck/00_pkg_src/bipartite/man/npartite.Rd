\name{npartite}

\alias{npartite}

\title{Computes indices for a masked-one-mode network}

\description{
 Computes network indices for networks that are represented as stacks of bipartite networks (``npartite''), thus having many impossible links.
}

\usage{
npartite(as.one.mode.web, index=c("connectance", "links per species", "mean degree"))
}

\arguments{
  \item{as.one.mode.web}{a one-mode representation of a network; if this is based on several bipartite networks, i.e. species in levels NOT interacting within that level OR without interactions between separated levels (e.g. 1st and 3rd) then such forbidden links must be represented as NA! For bipartite networks this can easily be achieved using as.one.mode(., fill=NA), but obviously for bipartite networks this function is obsolete.}
  \item{index}{Vector of names of indices to be computed for the npartite network. Currently only \option{connectance} and \option{links per species} (=\option{mean degree}) are available.}
}

\note{
An attempt is made to ensure that these indices converge to the same value as returned for a bipartite network!
}

\value{
A named list of indices.
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} }

\seealso{\code{\link{as.one.mode}}.}

\examples{
image(aomw <- as.one.mode(Safariland, fill=NA))
npartite(aomw)
networklevel(Safariland, index=c("connectance", "links per species"))
}

\keyword{package}


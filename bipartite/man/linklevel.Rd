\encoding{UTF-8}
\name{linklevel}
\alias{linklevel}

\title{Indices of a bipartite network at the link level}

\description{
Computes various indices of a network at the link-level, i.e. for each cell of the network matrix}

\usage{
linklevel(web, index=c("dependence", "endpoint"))
}

\arguments{
  \item{web}{ A matrix with pollinators as columns and plants as rows. Alternatively, when used on e.g. species occurrences across islands, rows are islands. }
  \item{index}{Vector of indices to be computed at the link level; options are: 
  	\itemize{
		\item \option{dependence} to compute dependence-matrix for each group level;
		\item \option{endpoint} to compute end-point degrees following Barratt et al. (2004).
	}
   }
 }

\details{
	For summaries of such indices see \code{\link{networklevel}}. It's still early days for this function ...
}

\value{
Returns a list of indices, each entry being a matrix of the same dimensions as the input web.
}

\references{ 
Barrat, A., Barthélemy, M., Pastor-Satorras, R. & Vespignani, A. (2004) The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences of the USA} \bold{101}, 3747-–3752
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} }

\note{ This function aims to facilitate analyses at the link level. So far, most studies at this level did not correct for the fact that observations per cell are clearly non-independent, nor for the abundances of the species, which also greatly affect indices. Room for improvement, but little room for new findings ...
}

\seealso{
	\code{\link{endpoint}}, \code{\link{specieslevel}}
}

\examples{
data(Safariland)
linklevel(Safariland)
}

\keyword{package}

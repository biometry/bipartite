\encoding{UTF-8}
\name{endpoint}
\alias{endpoint}

\title{Computes end-point degrees for a bipartite network}

\description{
Computes end-point degrees for a bipartite network, following the suggestion of Barrat et al. (2004)
}

\usage{
endpoint(web)
}

\arguments{
  \item{web}{ A matrix with pollinators as columns and plants as rows. Alternatively, when used on e.g. species occurrences across islands, rows are islands. }
 }

\details{Computation follows the outline of Gitarranz et al. (2004): ``the product k_i k_j of the degree of the two nodes connected by that link''. We then set additionally endpoint degrees for all non-existing links to 0! Thus, only for existing links endpoint degrees are computed. This is (to me) not obvious from the description in Gitarranz et al. (2004).
}

\value{
A matrix of end-point degrees
}

\references{ 
Barrat, A., M. Barthélemy, R. Pastor-Satorras, and A. Vespignani. 2004. The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences of the USA} \bold{101}, 3747-–3752. doi: 10.1073/pnas.0400087101.

Gilarranz, L. J., J. M. Pastor, and J. Galeano. 2011. The architecture of weighted mutualistic networks. \emph{Oikos} \bold{121}, 1154–-1162. doi: 10.1111/j.1600-0706.2011.19592.x.
}

\author{ Carsten F. Dormann }

\note{ This approach is, AFAIK, not tested by simulation; whether it is useful has still to be shown.
}

\examples{
# reproduces the example of Gitarranz et al. (2011):
data(memmott1999)
ends <- endpoint(memmott1999)
weights.mean <- tapply(memmott1999, ends, mean)
ends.weights <- tapply(ends, ends, mean)
plot(weights.mean, ends.weights, log="xy", pch=16)
}
\keyword{package}

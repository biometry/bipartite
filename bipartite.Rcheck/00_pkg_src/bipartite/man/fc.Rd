\encoding{UTF-8}
\name{fc}
\alias{fd} %its old name
\alias{fc}

\title{Calculates the functional complementarity for the rows of a web}
\description{
A community-level measure of ecological niche complementarity measured as the total branch length of a functional dendrogram based on qualitative differences in visitor assemblages between plants.
}
\usage{
fc(web, dist="euclidean", method="average", weighted=TRUE)
}

\arguments{
  \item{web}{A bipartite interaction web, i.e.~a matrix with higher (cols) and lower (rows) trophic levels.}
  \item{dist}{A distance accepted by the function \code{\link{dist}}.}
  \item{method}{The linkage method to be used by \code{\link{hclust}}. Any option accepted by \code{\link{hclust}} is allowed; defaults to \option{"average"}.}
  \item{weighted}{Option to analyse the web as binary or as weighted network. Default is \option{weighted=TRUE}, but analysis presented in Devoto et al. (2012) uses \option{weighted=FALSE}.}
}

\details{
fc measures community-level ecological niche complementarity as the total branch length of a functional dendrogram based on qualitative differences in visitor assemblages between plants. For details see Devoto et al. (2012).
}

\value{
The value of fc, which is not standardised and lies anywhere between 0 and a large number.
}


\references{ 
Devoto M., Bailey S., Craze P., and Memmott J. (2012) Understanding and planning ecological restoration of plant-pollinator networks.
\emph{Ecology Letters} \bold{15}, 319--328
}

\author{Mariano Devoto \email{mdevoto@agro.uba.ar}}

\seealso{\code{\link{networklevel}}, which uses this function.}

\examples{
data(Safariland)
fc(Safariland)    
fc(t(Safariland), dist="canberra", method="complete")
}

\keyword{ package }


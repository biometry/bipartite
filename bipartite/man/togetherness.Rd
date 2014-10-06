\name{togetherness}

\alias{togetherness}

\title{ Calculates the number of identical co-presences and co-absences for species-on-islands }

\description{
  Togetherness, or T-score, describes the level of similarity in the distributional pattern of two species. Originally proposed by Stone & Roberts (1992) for biogeographical situations can it also be applied e.g. to pollinators on different host plants.
}

\usage{
togetherness(web, normalise=TRUE, FUN = mean, ...)
}

\arguments{
  \item{web}{ A matrix with binary or counted interactions/links, where the higher trophic level is represented by columns.}
  \item{normalise}{Logical; shall index be normalised to a range of 0-1?}
  \item{FUN}{ The function to summarise species-pair T-scores with; defaults to mean.}
  \item{\dots}{ Arguments passed on to \code{FUN}, especially \code{na.rm=T} or colours for \code{hist}. }
}

\value{
  Returns the average (default) togetherness of all species combinations.   
}

\references{ 
Stone, L. and Roberts, A. (1992) Competitive exclusion, or species aggregation? An aid in deciding. \emph{Oecologia} \bold{91}, 419--424
}
\author{ Carsten F. Dormann}

\seealso{ \code{\link{C.score}} for another of Stone & Roberts' indices. }

\examples{
(m <- matrix(c(0,1,0,0,1,1,0,1,1,0), ncol=2, byrow=TRUE))
togetherness(m)
# or, with two togethernesses:
(n <- matrix(c(0,1,1,0,1,1,0,0,1,1,0,1,0,1), ncol=2, byrow=TRUE))
togetherness(n, normalise=FALSE)

data(Safariland)
togetherness(m)
}
\keyword{ package}


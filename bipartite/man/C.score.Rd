\name{C.score}
\alias{C.score}

\title{ Calculates the (normalised) mean number of checkerboard combinations (C-score) in a matrix }

\description{
Calculates the C-score for all pollinator species; the C-score represents the average number of checkerboard units for each unique species pair.
}

\usage{
C.score(web, normalise = TRUE, FUN = mean, ...)
}

\arguments{
  \item{web}{ A matrix with pollinators as columns and plants as rows. Alternatively, when used on e.g. species occurrences across islands, rows are islands. }
  \item{normalise}{ Logical; if \code{TRUE} (default), the C-score is ranged between 0 (no checkerboards) and 1 (only checkerboards). For \code{FALSE} the standard value of mean number of checkerboard pairs is returned. This is somewhat awkward for comparing different data sets, that's what the normalisation is for.}
  \item{FUN}{ Function to use when summarising the C-scores for each pairwise comparison. Defaults to \code{mean}, but other useful functions could be \code{median} (because C-scores are rather skewed) or \code{hist} (for a nice graph). }
  \item{\dots}{ Options to be passed on to FUN, e.g. \option{na.rm=T} for matrices with many zeros and \option{normalise=TRUE}. }
}

\details{
 As a first step, any quantitative matrix is converted to a binary matrix of presences and absences.

 Then, the formula given in Stone and Roberts (1990) is calculated for all species combinations, by calling \code{designdist} from the package \pkg{vegan}. See code for details.
}

\value{
Returns whatever the \option{FUN} produces as output. Default would be a single value, i.e. the mean C-score of the web.
}

\references{ 
Gotelli, N.J. and Rohde, K. (2002) Co-occurrence of ectoparasites of marine fishes: a null model analysis. \emph{Ecology Letters} \bold{5}, 86--94

Stone, L. and Roberts, A. (1990) The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79
}

\author{ Carsten F. Dormann }

\note{ 
The normalisation uses the following logic: let A and B be the number of occurrences of each of the two species, respectively. Then J is the number of sites where both occur and P the total number of sites.
The maximum value for C in each comparison is A*B, when they are exactly complementary and hence J=0. 
However, if (A+B)>length of vector(P), then there will be some co-occurrences and 
hence J>0=(A+B-P). The general maximum then becomes (A-A-B+P)(B-A-B+P)=(P-B)(P-A). 
For (A+B)<P, maximum is AB.

The minimum is set to 0.

}

\examples{
m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,TRUE)
C.score(m)
C.score(m, normalise=FALSE)
C.score(m, normalise=FALSE, FUN=print)

}
\keyword{package}

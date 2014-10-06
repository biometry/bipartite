\name{V.ratio}
\alias{V.ratio}

\title{ Calculates the variance-ratio as suggested by Schluter (1984) }
\description{
  A of species association is provided by the ratio of the variance in total species number (or total density of individuals) in samples to the sum of the variances of the individual species.
}
\usage{
V.ratio(web)
}

\arguments{
  \item{web}{ A matrix with pollinators in columns and plants in rows. For biogeographical applications: rows are islands (or sites). }
}
\details{
  This is a rather straight-forward index, which is described and evaluated extensively in Schluter (1984). He also warns against overinterpretation of the value. In principle, V-ratios larger than 1 indicate positive, smaller than 1 negative associations. Ecologically, competition can lead to small or large values, depending on their exact effects (see discussion in the Schluter paper).
}
\value{
Returns the V-ratio, i.e. a single value representing the ratio of variance in species number and variance in individual numbers within species.
}
\references{ 
Gotelli, N.J. and Rohde, K. (2002) Co-occurrence of ectoparasites of marine fishes: a null model analysis. \emph{Ecology Letters} \bold{5}, 86--94

Schluter, D. (1984) A variance test for detecting species associations, with some example applications. \emph{Ecology} \bold{65}, 998--1005

 }
\author{ Carsten F. Dormann}
\note{ 
Any quantitative matrix is first transformed into a binary (presence-absence) matrix!

Do not interpret without first reading the paper! It's worth it! See also applications in other studies, such as Gotelli and Rohde (2002).
}
\seealso{ \code{\link{C.score}} for another measure of species associations. }

\examples{
data(Safariland)
V.ratio(Safariland)
}
\keyword{ package }

\encoding{UTF-8}

\name{shuffle.web}
\alias{shuffle.web}

\title{ Shuffle web entries  }
\description{
Shuffles (= relocates) entries in a web matrix whilst maintaining the dimensionality.
}
\usage{
shuffle.web(web, N, legacy=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{web}{ An interaction matrix.}
  \item{N}{ Number of desired shuffled matrices. }
  \item{legacy}{Logical; use the old or new algorithm? Defaults to TRUE, i.e. the old algorithm. The new algorithm was written by Paul Rabie and is about 3 times faster (due to avoiding a loop). For consistency reasons, the old, slow algorithm remains the default.}
}
\details{
 This function is designed to behave similar to \code{r2dtable}, i.e. it returns a list of randomised matrices. 
 In contrast to \code{r2dtable} is does not keep marginal sums constant! 
 
 This function is thought of as a null model for the analysis of bipartite webs. It keeps two web properties
 constant: The number of interactions and the number of links (and hence connectance). A comparison of
 \code{shuffle.web}- and \code{r2dtable}-based webs allows to elucidate the effect of marginal sums.
}
\value{
 A list of N randomised matrices with the same dimensions as the initial web.
}

\references{
 This null model can be thought of as a quantitative version of Fortuna & Bascompte (2006) \dQuote{null model 1}:

 Fortuna, M. A., and J. Bascompte. 2006. Habitat loss and the structure of plant-animal mutualistic networks. Ecology Letters 9: 281-286.

 For a very nice and thorough overview of null models in general see:

 Gotelli, N. J., and G. R. Graves. 1996. Null Models in Ecology. Smithsonian Institution Press, Washington D.C.
 
 For null models and their application to webs/networks see, e.g.:

 Vázquez, D. P., and M. A. Aizen. 2003. Null model analyses of specialization in plant-pollinator interactions. Ecology 84: 2493-2501.

 Vázquez, D. P., C. J. Melian, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction 
strength in ecological networks. Oikos 116: 1120-1127.

}

\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }

\note{ \code{shuffle.web} is not an extremely intelligent nullmodel. You may want to think of a better one for your specific application!
}

\seealso{ \code{\link{r2dtable}} }

\examples{

data(Safariland)

shuffle.web(Safariland, N=2)

}
\keyword{ package }


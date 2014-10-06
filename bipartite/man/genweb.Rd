\encoding{UTF-8}
\name{genweb}
\alias{genweb}

\title{ Generate a random bipartite web }
\description{
 Generates a random bipartite web, based on \code{r2dtable} and lognormal marginal distributions.
}
\usage{
genweb(N1 = 10, N2 = 30, dens = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{ Number of species in the lower trophic level; or a vector of length 2 giving number of lower and higher trophic level species. Defaults to 10. Large values (>70) will take much longer!}
  \item{N2}{ Number of species in the higher trophic level. Ignored if N1 a vector of length 2. Defaults to 30. Large values (>70) will take much longer!}
  \item{dens}{ Interaction density, i.e. how many interactions there shall be, on average, for each \emph{link}. Defaults to 2 (the median observed interaction density in the NCEAS pollination webs). Large values (> 10) will take much longer to find matrices for.}
}
\details{
  This function can be used to create simple, but not necessarily realistic, bipartite webs for given dimensionality and interaction density. Marginal distributions are  assumed to be lognormal, mean and standard deviation are calculated from \option{N1}, \option{N2} and 
  \option{dens} (see code for details).
}
\value{
  A matrix with N1 x N2 species.
}
\author{ Jochen Fründ and Carsten F. Dormann}

\note{ Can be a bit time-consuming for large webs, because the absolute values for both dimensions have to match perfectly. This involves a rather inelegant \code{while}-loop.}

\examples{
genweb()
}
\keyword{ package}


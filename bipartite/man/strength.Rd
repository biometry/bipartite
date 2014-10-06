\encoding{UTF-8}
\name{strength}

\alias{strength}

\title{Computes species strength according to either of two definitions}

\description{
Computes species strength of the higher level species as a measure of how important a species is in the network
}

\usage{ strength(web, type="Bascompte")
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species as columns and number of interactions as entries.}
  \item{type}{Which definition of species strength should be used, \option{Bascompte} (default) or \option{Barrat}? See Details for definitions.}
}

\details{
There are two definitions of species strength, that of Bascompte et al. (2006) as the sum of dependencies of a species, and that of Barrat et al. (2004) as the weighted sum of links. As a consequence, Bascompte et al.'s strength sums to the number of species in the other group, while Barrat et al.'s strength is simply the number of interactions, a trivial measure of a species importance. 

In contrast to the claim of Gilarranz et al. (2012, p. 1155), this definition of strength gives no information of the centrality of a species within a network structure (and neither does Bascompte et al.'s).}

\value{
A vector of species strengths for the higher level. Employ this function on the transpose of the web to compute the strengths of the lower level (see example).
}

\references{
Barrat, A., Barthélemy, M., Pastor-Satorras, R. & Vespignani, A. (2004) The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences of the USA} \bold{101}, 3747-–3752

Bascompte, J., Jordano, P. & Olesen, J.M. (2006) Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431-–433

Gilarranz, L.J., Pastor, J.M. & Galeano, J. (2012) The architecture of weighted mutualistic networks. \emph{Oikos} \bold{121}, 1154-–1162

}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} }

\seealso{
  \code{\link{specieslevel}} which could (but doesn't yet) call \code{strength} (instead it uses the default always)
}

\examples{
data(Safariland)
s1 <- strength(Safariland, type="Barrat")
s2 <- strength(Safariland, type="Bascompte")
plot(s1, s2, log="x")
cor.test(s1, s2, type="ken")
# for lower level:
strength(t(Safariland))

}

\keyword{package}


\encoding{UTF-8}
\name{czvalues}
\alias{czvalues}

\title{Computes c and z for network modules}

\description{
Function to compute c and z values of module members according to Guimerà & Amaral (2005), with formulae taken from Olesen et al. (2007)}
\usage{
czvalues(moduleWebObject, weighted=FALSE, level="higher")
}

\arguments{
  \item{moduleWebObject}{A \code{moduleWeb-class} object as created by \code{\link{computeModules}}.}
  \item{weighted}{logical; if TRUE computes c and z from quantitative (=weighted) data; in this case, it will compute strength, rather than degrees for each species.}
  \item{level}{"higher" or "lower" trophic level to compute c and z values for; defaults to "higher"}
}

\details{
c = 1 - sum( (k.it/k.i)^2)    # among-module connectivity = participation coefficient P in Guimerà & Amaral

z = (k.is - ks.bar) / SD.ks  # within-module degree

k.is = number of links of i to other species in its own module s; 
ks.bar = average k.is of all species in module s; 
SD.ks = standard deviation of k.is of all species in module s; 
k.it = number of links of species i to module t; 
k.i = degree of species i

Note that for any species alone (in its level) in a module the z-value will be NaN, since then SD.ks is 0. This is a limitation of the way the z-value is defined (in multiples of degree/strength standard deviations).

Olesen et al. (2006) give critical c and z values of 0.62 and 2.6, respectively. Species exceeding these values are deemed connectors or hubs of a network. The justification of these thresholds remains unclear to me. They may also not apply for the quantitative version.
}

\value{
A list with two vectors, c and z, for all species of the selected trophic level. 
}

\note{
These indices were developed for one-mode networks; we'll have to see whether they make sense for bipartite networks, too! In particular, note that this function is based on a higher trophic level perspective. While the modules are identified using both trophic levels, c and z are computed through the strengths/degrees of only one trophic level. It would be desirable to have a truly two-level version. Since there are usually very different numbers of species in the two trophic levels, simply averaging the values of each trophic level won't do. But maybe a weighted average?

Consider the following problem for computing c and z for the higher trophic level: For modules with only one species from the lower trophic level, the z-values will be NaN, since SD.ks is 0!  I decided to SET these values to 0, since they only occur when all species in that module will have the same number of links (which is obviously the case when there is only one lower-level species). Then the numerator is also 0. Thus, the value of 0 indicates that this species has no deviation from the rest of the module members (which is what I think z is supposed to represent).

Since \code{computeModules} is experimental, also \code{czvalues} may not always work (i.e. if object \code{mod} is corrupted or ill-formed).
}

\references{ 
Guimerà, R. and Amaral, L.A.N. (2005) Functional cartography of complex metabolic networks. \emph{Nature} \bold{433}, 895--900.

Olesen, J.M., Bascompte, J., Dupont, Y.L. and Jordano, P. (2007) The modularity of pollination networks. \emph{Proceedings of the National Academy of Sciences of the USA} \bold{104}, 19891-19896.
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}, 20 Mar 2012}

\examples{
\dontrun{
data(memmott1999)
set.seed(2)
mod <- computeModules(memmott1999, steps=1E4)
cz <- czvalues(mod)
plot(cz[[1]], cz[[2]], pch=16, xlab="c", ylab="z", cex=0.8, xlim=c(0,1), las=1)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5)   # dito
text(cz[[1]], cz[[2]], names(cz[[1]]), pos=4, cex=0.7)
}
}

\keyword{ package }


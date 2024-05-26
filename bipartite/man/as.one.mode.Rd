\name{as.one.mode}

\alias{as.one.mode}

\title{Conversion of a network matrix}

\description{
 This helper function converts a bipartite matrix into a one-mode matrix.
}

\usage{
as.one.mode(web, fill = 0, project="full", weighted=TRUE, legacy=FALSE)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species 
    as columns and number of interactions as entries.}
  \item{fill}{What shall unobserved combinations be represented as in the one-mode matrix (see below)? Defaults to 0. Set to NA if links not possible for bipartite networks should be masked (i.e. those within a level).}
  \item{project}{There are different ways to convert a two-mode (bipartite) network into one-mode networks. The most common is to focus on one set (e.g. the n pollinators) and compute a n x n matrix with entries between species that pollinate the same plant (\option{"higher""}). Similarly, one can compute a k x k matrix for the k plant species (\option{"lower"}). Or, finally and the default, one can compute an (n+k) x (n+k) matrix in which only the observed interactions are present (\option{"full"}). This is in fact a near-trivial, symmetric matrix with 0s between species of the same trophic level.}
  \item{weighted}{Logical; shall the strength of links be included in the one-mode output? Defaults to TRUE, but can be set to FALSE to turn a weighted two-mode into a binary one-mode network.}
  \item{legacy}{Logical, set to FALSE by default. Until version 2.21, projections to either higher or lower were based on reasoning what the correct weights should be from the perspective of an ecologist: If, say, species A and B interact with species a-z, then the link weight would be the sum of the minima of the two vectors A with a-z and B with a-z. In other words, the minimal overlap. In network science, the projection is defined differently, however, namely as the maximum! Thus, bringing \pkg{bipartite}'s \code{as.one.mode} closer to common practice, we now use the proper definition as default. The old definition remains available by setting \option{legacy=TRUE}. The new definition is: \eqn{AA^T-\textrm{diag}(A)}, follwing Thurner et al. (2018, p. 149).}
}

\details{
In bipartite (or: two-mode) networks, participants are of different types (e.g. pollinators and plants, actors and parties in social research). Hence, a party cannot connect to another party except through actors. A pollinator interacts with another pollinator only through the host plant.
  
Much network theory, however, is based on one-mode networks, where all participants are listed in one vector, i.e. plants and pollinators alike, actors together with events. This function here transforms the more condensed bipartite representation into a one-mode-representation, filling the unobserved type of interactions (i.e. plants with plants and pollinators with pollinators) with 0 (unless you specify it differently in \option{fill}).

The lower trophic level (e.g. plants or rows) is listed first, then the higher trophic level (e.g. pollinators or columns). Hence, pollinator 2 becomes species number r+2, where r is the number of rows of the network matrix.

In addition to the "full" projection, there are "inner" projections, yielding a network \bold{only} of the lower or higher level (hence the argument \option{project="lower"/"higher"/"full"}). Such an inner projection inevitably loses information: if two pollinators pollinate three plant species, then they are connected in such a projection through 3 links. The weight of each link will be different, but in the projection only one weight can be given. This is where the information is lost. Several indices (betweenness, centrality) depend on one-mode projections of this kind. Still, the user should always ask herself, whether the projection might not have unintended consequences!

If \option{weighted=TRUE}, then the returned one-mode network contains the parallel maximum (or minimum, if \option{legacy=TRUE}) of the observed interactions between two species. That means, if two species A and B interact with species 1 to 5 in the other group, then the two interaction vectors for A with 1 to 5 and B with 1 to 5 are placed next to each other, and for every species 1 to 5 the maximum for each of these 5 values for the two vectors is retained (the parallel maximum). The idea is that the similarity between A and B is driven by their largest communality in interactions. Next, the five parallel maximal values are added to yield the final weight for this link.
  
The benefit of this conversion is access to the wonderful R-package Social Network Analysis (\pkg{sna}), with its many one-mode indices (such as \code{betweenness}, \code{closeness}, \code{centralization}, \code{degree}, \code{kpath.census} and so forth). Furthermore, \code{gplot} in that package also provides cool network depictions well worth checking out.
  
With respect to \pkg{bipartite}, \code{as.one.mode} is employed in the function \code{\link{nodespec}}, which itself uses the \pkg{sna}-function \code{geodist}.
}

\value{
  A matrix of dimension (n+k) x (n+k), where n and k are the dimensions of the input web. Both dimensions are given the names of the original web (first the lower, then the higher trophic level).
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} }

\seealso{Function \code{projecting_tm} in package \pkg{tnet} provide an analogous ways of converting two-modes into one-modes. This function can be accessed after transforming the web-matrix into an edge list using \code{\link{web2edges}}.}

\references{
Thurner, S., Hanel, R. and Klimek, P. (2018) Introduction to the Theory of Complex Systems. Oxford University Press, Oxford
}

\examples{
data(Safariland)
image(Safariland)
image(as.one.mode(Safariland))
par(xpd=TRUE, mar=c(0, 6, 0, 6))
gplot(as.one.mode(Safariland, project="lower"), 
	label=rownames(Safariland), gmode="graph", 
	label.cex=0.6, vertex.cex=2, vertex.col="green")
}

\keyword{package}


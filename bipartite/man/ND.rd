\encoding{UTF-8}
\name{ND}

\alias{ND}
\alias{BC}
\alias{CC}

\title{Normalised degree, betweenness and closeness centrality}

\description{
 Calculates normalised degrees, and two measures of centrality, betweenness and closeness. These two are based on one-mode representations of the network and invoke functions from \pkg{sna}.
}

\usage{
ND(web, normalised=TRUE)
BC(web, rescale=TRUE, cmode="undirected", weighted=TRUE, ...)
CC(web, cmode="suminvundir", rescale=TRUE, ...)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species
    as columns and number of interactions as entries.}
  \item{normalised}{Shall the degrees be normalised? If so (default), the degree for a species is divided by the number of species in the other level (see, e.g., Martín González et al. 2010).}
  \item{rescale}{If TRUE (default), centrality scores are rescaled such that they sum to 1.}
  \item{cmode}{String indicating the type of betweenness/closeness centrality being computed (directed or undirected geodesics, or a variant form - see help for \code{closeness} and \code{betweenness} in \pkg{sna} for details). The default, \option{"suminvundir"} for \code{CC} and \option{"undirected"} for \code{BC}, uses a formula that can also be applied to disconnected (=compartmented) graphs. Other cmodes may not.}
  \item{weighted}{Logical; if TRUE, bipartite projection will include edge weights, i.e. number of interactions. Defaults to TRUE.}
  \item{...}{Options passed on to \code{betweenness} and \code{closeness}, respectively. Notice that in particular the option \option{ignore.eval=FALSE} will yield VERY different values than the default. BC and CC use defaults of sna::betweenness and sna::closeness, respectively, but that does not imply that these settings are per se the best! (Thanks to Michael Pocock for drawing my attention to this issue!)}
}

\details{
 These functions are convenience functions to enable easy reproduction of the type of analyses by Martín González et al. (2010). BC and CC are wrappers calling two functions from \pkg{sna}, which uses one-mode, rather than bipartite data. 
  
 One-mode projections of two-mode networks are carried out by assigning a link to two species that share a interaction with a member of the other set (plant in case of pollinators, or pollinators in case of plants). There are different ways to do this (see \code{\link{as.one.mode}}), and many authors do not communicate well, which approach they have taken.
 
 If the network is fully connected, all species of the same level will be linked to each other through only one step and hence have the same betweenness. This leads to values of 0.

BC reflects the number of shortest paths going through the focal node. CC is the inverse of the average distance from the focal node to all other nodes.
% BE AWARE that there are two definitions of closeness centrality, one being the inverse of the other! The networkX homepage defines CC as the inverse of the average distance from the focal node to all other nodes (\url{http://networkx.lanl.gov/reference/generated/networkx.closeness_centrality.html#networkx.closeness_centrality}), while Wikipedia defines CC simply as the average distance itself (\url{http://en.wikipedia.org/wiki/Centrality}). In \pkg{sna} the first definition is implemented, and this makes also more sense to me: closeness should be higher for central nodes.

Both BC and CC can be normalised so that they sum to 1 (using \option{rescale=TRUE}). This only affects the absolute values, but not the qualitative results.

The interested user may want to also have a look at the networkX homepage (\url{http://networkx.lanl.gov}) for an excellent, open, Python-based tool to analyse, depict and manipulate (one-mode) networks. It is not specifically meant for bipartite networks such as this package, though.
}

\value{
  A list with two entries, ``lower'' and ``higher'', which contain a named vector of normalised degrees, betweenness centrality and closeness centrality, respectively. The lower-entry contains the lower trophic level species, the higher analogously the higher trophic level species.
}

\references{
 Martín Gonzáles, A.M., Dalsgaard, B. and Olesen, J.M. 2010. Centrality measures and the importance of generalist species in pollination networks. \emph{Ecological Complexity} \bold{ 7}, 36--41
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} }

\note{ 
Experimental. Should work most of the time, but not necessarily always. Also, on trials with the same data as those of Martín González et al. (2010), numerical values differed. Whether this is due to rounding errors, different non-linear least square fits in JMP and R or whatever I cannot tell. See example for my attempt to reproduce their values for the network ``Azores'' (aka \code{\link{olesen2002flores}}).
}

\seealso{
  \code{centralization}, \code{betweenness} and \code{closeness} in \pkg{sna}; \code{\link{specieslevel}} which calls them
}

\examples{
## example:
data(olesen2002flores)
(ndi <- ND(olesen2002flores))
(cci <- CC(olesen2002flores))
(bci <- BC(olesen2002flores))

cor.test(bci[[1]], ndi[[1]], method="spear") # 0.532
cor.test(cci[[1]], ndi[[1]], method="spear") # 0.403

cor.test(bci[[2]], ndi[[2]], method="spear") # 0.738
cor.test(cci[[2]], ndi[[2]], method="spear") # 0.827
\dontrun{
## PLANTS:
bc <- bci[[1]]
cc <- cci[[1]]
nd <- ndi[[1]]
# CC:
summary(nls(cc ~ a*nd+b, start=list(a=1,b=1))) # lower RSE
summary(nls(cc ~ c*nd^d, start=list(c=0.072,d=0.2))) 
# BC:
summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
summary(nls(bc ~ c*nd^d, start=list(c=2,d=2))) # lower RSE

## ANIMALS:
bc <- bci[[2]]
cc <- cci[[2]]
nd <- ndi[[2]]
# CC:
summary(nls(cc ~ a*nd+b, start=list(a=1,b=1)))  
summary(nls(cc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE 
# BC:
summary(nls(bc ~ a*nd+b, start=list(a=1,b=1)))
summary(nls(bc ~ c*nd^d, start=list(c=0.2,d=2))) # lower RSE
}
}

\keyword{package}


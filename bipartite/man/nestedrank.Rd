\encoding{UTF-8}

\name{nestedrank}
\alias{nestedrank}

\title{
Calculates the rank of a species in a matrix sorted for maximum nestedness
}

\description{
Ranks species according to their generality, which is measured as the position in the nestedness matrix. A generalist will interact with more species and thus have a rank closer to 1, while specialists (and rare species) will have ranks with higher values.
}

\usage{
nestedrank(web, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=FALSE)
}

\arguments{
  \item{web}{A matrix with elements of a set (e.g., plants) as rows, elements of a second set (e.g., pollinators) as columns and number of interactions as entries.}
  \item{method}{One or more of the following: \option{NODF}, \option{nodf}, \option{binmatnest}, \option{wine}, \option{sort}. See details for details on each method.}
  \item{weighted}{For NODF and wine only: should the number of interactions per link be used as weights? See help of \code{nestednodf} in \pkg{vegan} for details.}
  \item{normalise}{Logical; defaulting to TRUE. Divides the rank-1 by the number of species -1, thereby ranging it between 0 (most generalist) and 1 (most specialised).}
  \item{return.matrix}{Logical, defaulting to FALSE. Should the matrix resulting from the nestedness-sorting be returned as well?}
}

\details{
The idea is to re-arrange the network matrix according to its nestedness, so that the most ``generalist'' species with most links will be in the first row/column and decreasing from there. The nestedness matrix can be computed in different ways. There are four different methods currently available:
  \describe{
  \item{NODF (or nodf)}{will use \pkg{vegan}'s  \code{nestednodf}-function to arrange the matrix. With \option{weighted=TRUE}, which is the default, it will use the actual number of interactions, rather than the number of links}
  \item{binmatnest}{will use the \pkg{vegan}'s \code{nestedtemp}-function to arrange the matrix. This is only using binary information, so weighting has no effect.}
  \item{wine}{will use the \code{\link{wine}}-function to arrange the matrix. When \option{weighted=FALSE}, wine will be applied to a binary matrix.}
  \item{sort}{will simply sort the matrix by marginal totals (i.e. by number of interactions per species when \option{weighted=TRUE} or by number links (=degree) when \option{weighted=FALSE}. In this case the rank simply represents the abundance of the species in this network.}
    }

}

\value{
A list of nestedness ranks vectors for the lower and higher trophic level (smallest value for the most generalist). If \option{return.matrix=TRUE}, a third list entry will contain the nested matrix.
}

\note{
Since nestedness is itself not a straight-forward measure of something ecologically meaningful, also these ranks may or may not be. At least there is a high chance that they represent merely abundance of each species. See example for an idea on how to check for the effect of abundance as such.
}

\references{
Alarcon, R., Waser, N.M. and Ollerton, J. 2008. Year-to-year variation in the topology of a plant-pollinator interaction network. \emph{Oikos} \bold{117}, 1796--1807
}

\author{
Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}
}


\seealso{
\code{\link{nested}}; \code{\link{nestedrank}} is called by \code{\link{specieslevel}}
}
\examples{
\dontrun{
ranks <- sapply(c("nodf", "binmatnest", "wine", "sort"), function(x) 
  nestedrank(Safariland, method=x)[[2]])
cor(ranks) # high correlation between sort and other indicate that only abundance matters
}
}

\keyword{ package }


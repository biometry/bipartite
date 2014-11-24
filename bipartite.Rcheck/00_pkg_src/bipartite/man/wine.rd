\encoding{latin1}

\name{wine}
\alias{wine}
\alias{plot.wine}

\title{
Weighted-Interaction Nestedness Estimator
}
\description{
Calculates the nestedness of a network taking into account the weight of the interactions, according to the method proposed by Galeano et al. (2008).
}
\usage{
wine(web, nreps = 1)
\method{plot}{wine}(x, ...)
}

\arguments{
  \item{web}{A matrix with elements of a set (e.g., plants) as rows, 
           elements of a second set (e.g., pollinators) as columns and number 
           of interactions as entries.}
  \item{nreps}{Number of replicates for constructing random networks.}
  \item{x}{ An object resulting of applying wine function}
  \item{\dots}{Additional graphical parameters to \code{\link[fields]{image.plot}}}
}
\details{
Nestedness estimators use presence-absence (binary) adjacency matrices as the basis for calculating nestedness,
as they provide a simple description and characterization of the topology of the network. However, networks are
specified not only by their topology but also by the heterogeneity in the weight (or the intensity) of the connections
(Barrat et al., 2004). Characterizing links just with presence-absence data does not take into account the possible
differences in intensity among links. WINE (Weighted-Interaction Nestedness Estimator) is a new nestedness estimator
that takes into account the weight or intensity of each interaction (e.g., in a plant-pollinator network, the number
of registered visits of a particular interaction). Thus, instead of using presence-absence matrices, 
WINE calculates nestedness from quantitative data matrices that include the number of events of each interaction.
This is the first estimator that allows for the characterization of weighted nestedness. 
WINE calculates a nestedness value that approaches zero when the nestedness pattern of the original data matrix
is close that of equivalent random matrices, and it approaches one as it gets closer to the nestedness 
of the maximal nestedness matrix. Thus, this estimator evaluates the relative position of the data matrix
between the corresponding random matrices and the maximal nestedness matrix. Negatives values for this estimator
can be found in some synthetic matrices that have been described as `anti-nestedness' matrices (Almeida-Neto et al. 2007).

The calculation of the weighted-interaction nestedness estimator starts with the matrix containing the number of events
of each interaction, Mij. The matrix is packed by arranging rows and columns from top to bottom and from left to right,
respectively, in ascending order according to their marginal totals. Nestedness is related to the proximity of existing 
links to one another in the packed matrix, so that the most nested matrix is the one that after packing shows a minimum 
mixing of filled cells (links) with empty cells (no links) (Corso et al., 2008, Ulrich et al., 2009).
WINE is based on the concept of estimating nestedness through the calculation of a Manhattan distance from each of the 
matrix cells containing a link to the cell corresponding to the intersection of the row and columns with the lowest 
marginal totals (number of links). This concept resembles in a way the one used by Corso et al. (2008), although the 
distances are measured to the opposite corner of the packed matrix. Additionally, in WINE, the Manhattan distance is 
replaced by a weighted Manhattan distance. 
The statistical significance of this nestedness index value is tested against a null model that constrains matrix 
fill to observed values, retains the distribution of number of events in the links but does not constrain marginal totals. 
Further details can be found in Galeano et al. (2008).

}
\value{
wine returns an object of class \code{wine}, basically a list with the following components:
\item{win}{Weighted-interaction nestedness of dataset (WIN)}
\item{wine}{Weighted-interaction nestedness estimator (WINE): The weighted-interaction nestedness estimator value. It will be 0 for random distribution and 1 for maximal nestedness}
\item{zscore}{z-score of the weighted-interaction nestedness}
\item{pvalue}{probability of having a z-value equal to or greater than Z (from the tabulated value of the cumulative function). Values of p<0.05 indicate that the dataset is significantly nested.}
\item{dmax}{Weighted-interaction nestedness of the maximal nestedness matrix.}
\item{drnd}{Average weighted-interaction nestedness of random replicates}
\item{dij.w}{Matrix of dijw values.  These values provide a measure of the contribution of each interaction (link) to total nestedness}
\item{dij.max}{Maximal nestedness matrix}

The S3 plot method for wine displays dij.w in a coloured image plot where red cells have high weights in the network and blue cells have minimum weights.

}
\references{

Barrat, A., Barthelemy, M., Pastor-Satorras, R., and Vespignani, A. (2004)  The architecture of complex weighted networks.  \emph{PNAS} \bold{101}, 3747--3752

Corso G, de Araujo AIL, de Almeida AM (2008) A new nestedness estimator in community networks.
     \emph{arXiv} 0803.0007v1 [physics.bio-ph]

Galeano J, Pastor JM, Iriondo JM (2008) Weighted-Interaction Nestedness Estimator (WINE): A new 
   estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v2 [physics.bio-ph] 


Ulrich, W., Almeida-Neto, M., and Gotelli, N.J. (2009) A consumer's guide to nestedness analysis. \emph{Oikos} 118, 3-17 


Ulrich, W. and Gotelli, N.J. (2007) Null model analysis of species nestedness patterns. \emph{Ecology} 88, 1824-1831 

}
\author{
Marcelino de la Cruz \email{marcelino.delacruz@upm.es}, Juan M. Pastor \email{juanmanuel.pastor@upm.es}, 
Javier Galeano \email{javier.galeano@upm.es} and Jose M. Iriondo \email{jose.iriondo@urjc.es}

}
\note{
This is the first approach to a weighted nestedness and a full ecological interpretation of  its meaning is still lacking. It is not possible to perform a systematic comparison between this and other nestedness indices because the latter rely just on presence absence data  whereas the former feeds on a quantitative data matrix. For a well-performed comparison of other nestedness indices see Ulrich & Gotelli (2007). 

\code{wine} may return NaN for different parameters essentially for two different
reasons: 
a) if \option{nreps} is not specified, \code{wine} adopts \code{nreps=1} by default and NaN is
returned for z-score and p value. This is due to the fact that with \code{nreps=1} the variance of drnd is zero and z-score becomes infinite. The same outcome may occur in some cases with very low values of nreps. To ensure proper values of z-score and p-values \code{nreps=100} or higher is suggested.
b) if dw = drnd = dmax  \code{wine} equals 0/0, and if drnd = dmax \code{wine} tends to infinity. In both cases, NaN is returned by \code{wine}. This is more likely to occur in cases where the dimensions of the matrix are very low (e.g,  (\code{dim < c(4,4)}) because in those cases the number of possible values of dw, drnd and dmax is also reduced.

This is WINE version 3.2, available also in Matlab and C++ at \url{http://hypatia.agricolas.upm.es/WINE/WINE.html}.
}


\seealso{
\code{\link{nestedness}} and \code{\link{discrepancy}}.
}
\examples{
data(Safariland, package="bipartite")
safariland.w <- wine(Safariland, 10)
plot.wine(safariland.w)
}

\keyword{package}

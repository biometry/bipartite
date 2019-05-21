\encoding{UTF-8}

\name{decimalr2dtable}

\alias{decimalr2dtable}

\title{ Generates matrix with same marginal totals for non-negative decimal values}

\description{
Function to generate null model matrices for the cases when the matrix entries are non-negative numbers, possibly non-integers. It maintains marginal totals (as does \code{r2dtable}), but "smears" them out over all cells.
}

\usage{decimalr2dtable(N=10, web, steps=prod(dim(web)))}

\arguments{
  \item{N}{ Number of desired null model matrices. }
  \item{web}{An interaction matrix.}
  \item{steps}{Number of successive alterations of matrix entries. Defaults to 10 times the number of cells. See details.}
}

\details{
 This function is a cross between \code{\link{r2dtable}} and \code{\link{swap.web}}. Its output are N matrices of the same dimension as the input, and with same marginal totals, but with different allocation of values to the cells. 
 Here is what the algorithm does:
 \enumerate{
 \item Index a 2 x 2 submatrix by select randomly two rows and two columns.
 \item Draw a random value between 0 and the minimum of the diagonal entries in the submatrix.
 \item Subtract this value from the diagonals and add it to the counter-diagonal entries.
 \item Repeat \option{steps} times
 }
 The result is a decimal-numbered matrix, typically with values > 0 in each cell. 
 
 Indication: This function may be useful in some rare constellations. Imagine you sampled a plant-pollinator network and instead of counting the number of flower visits you recorded the nectar extracted by each pollinator. Then the marginal totals would indicate nectar production (plus confounding nectar attractiveness) and consumption potential for plants and pollinators, respectively. So, given that species differ in nectar production and consumption, what would you expect the network to look like? Enter \code{decimalr2dtable}.
 
 If external abundances (even in funny units such as biovolume in ml) are available, this function can easily provide the respective null models. See examples.
}

\note{
  The output will typically be a fully filled matrix! Computing any index sensitive to matrix filling (such as connectance, degree, nestedness) for such a matrix is non-sensical! 
  
  Also, if used as a null model, an implicit assumption is that the values in the original matrix are meaningful as marginal totals. This may often not be the case, for example if entries are rates. Thus, probably this function is of very limited usefulness in the context of network analyses!
}

\value{
 A list of N randomised matrices with the same dimensions as the initial web, all probably filled completely.
}


\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }

\seealso{ \code{\link[stats:r2dtable]{r2dtable}}, \code{\link{vaznull}}, \code{\link{shuffle.web}} and \code{\link{swap.web}}}

\examples{
obs <- networklevel(Safariland, index="generality")

nulls <- decimalr2dtable(10, Safariland)
g.dec <- sapply(nulls, networklevel, index="generality")
nullsint <- nullmodel(Safariland, N=10)
g.int <- sapply(nullsint, networklevel, index="generality")
plot(density(g.dec[1,]), xlim=c(1, 3))
lines(density(g.int[1,]), col="red")
abline(v=obs[1], col="green")


## If you want to use external abundances to set up your null model:
set.seed(1)
ext.rows <- runif(9) # imagine these are your external abundances for Safariland
ext.cols <- runif(27)
# standardise to sum = 1:
ext.rows <- ext.rows/sum(ext.rows)
ext.cols <- ext.rows/sum(ext.cols)
web <- tcrossprod(ext.rows, ext.cols) * sum(Safariland) 
# (to get to the same interaction density as original web)
#
# this can now be used as input for decimalr2dtable:
image(decimalr2dtable(N=1, web)[[1]]) # remember: white are high values!
}
\keyword{ package }


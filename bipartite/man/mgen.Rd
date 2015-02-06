\encoding{UTF-8}

\name{mgen}
\alias{mgen}

\title{Generate simulated network according to a given probability matrix}

\description{
Generic network simulating algorithm based on a probability matrix and a desired number of interactions.
}
\usage{
mgen(web, n=sum(web), keep.species=TRUE, rep.cell=TRUE, autotransform="sum", trials=100)
}
\arguments{
	\item{web}{a matrix with observation probabilities, emerging from some null model considerations external to this function; if an original network is used, this will be automatically converted to a probability matrix by dividing it by the number of interactions (CFD); ORIGINAL: a probability matrix}
	\item{n}{number of interactions to allocate into the new matrix}
	\item{keep.species}{Random assignment of interactions may lead to empty columns or rows and hence reduce the dimensions of the simulated web. By default, this is prevented from happening, i.e. each row/column will receive at least one interaction. Setting keep.species to FALSE may (but need not) cause a loss of species.}
\item{rep.cell}{whether cells can be assign interactions repeatedly; defaults to TRUE}
\item{autotransform}{determines how a non-probability web is converted into probabilities; option 'sum' simply divides each entry by the sum of interactions in the web;  option 'equiprobable' is the product of marginal probabilities (col/rowSums divided by sum(web) (in a sense this is the basis of the r2dtable null model, just without the 'turn into integers' bit)}
\item{trials}{allocating interactions (when rep.cell=TRUE) can be difficult or even impossible. When the number of attempts to allocate them exceeds trials * n it will be stopped and an error message will be returned. Default is 100. Setting 'trials' to a very large value may make this function run for hours and hours. Your choice!}
}

\details{
This is a generic function to build null models for mutualistic networks, used by Vázquez et al. (2009). It is general in the sense that it allows any type of probability matrix to be used for constructing the simulated matrices. It does not, however, constrain rown and column totals, nor does it constrain connectance.

If \option{rep.cell=TRUE}, repeated interactions are added, thus generating a quantitative matrix with cell values as positive integers. If \option{rep.cell=FALSE}, no repeated assignment of interactions is allowed, thus generating a binary matrix of zeros and ones. Note that when \option{rep.cell=FALSE} the number of interactions to be assigned must be equal or lower than the number of cells in the matrix.
}

\references{
Vázquez, D. P.; Chacoff, N. P. & Cagnolo, L. 2009. Evaluating multiple determinants of the structure of mutualistic networks Ecology. 90: 2039-2046
}

\examples{
\dontrun{
## Generate simulated matrix from homogeneous probability matrix
probmat <- matrix(1/100, 10, 10)
mgen(web=probmat, n=100)

## Generate binary matrix with probability proportional to degree
## of an observed binary matrix m
obs.mat <- matrix(c(1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0), 5, 5)
rs <- rowSums(obs.mat)
cs <- colSums(obs.mat)
web <- rs %*% t(cs)
web <- web/sum(web)
n = sum(obs.mat)
mgen(web, n, keep.species=FALSE, rep.cell=FALSE) # Allowing zero marginal sums
mgen(web, n, keep.species=TRUE, rep.cell=FALSE) # Not allowing zero marginal sums

## Generate quantitative matrix with probability proportional
## to interaction frequency in an observed matrix m (which is the default of 'autotransform'):
mgen(mosquin1967, keep.species=FALSE, rep.cell=TRUE) # Allowing zero marginal sums
mgen(mosquin1967, keep.species=TRUE, rep.cell=TRUE) # Not allowing zero marginal sums
}
}
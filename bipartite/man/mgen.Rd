\encoding{UTF-8}

\name{mgen}
\alias{mgen}

\title{Generate simulated network according to a given probability matrix}

\description{
Generic network simulating algorithm based on a probability matrix and a desired number of interactions.
}
\usage{
mgen(web, n=sum(web), keep.species=TRUE, rep.cell=TRUE, autotransform="sum_warn")
}
\arguments{
	\item{web}{a matrix with observation (interaction) probabilities, emerging from some considerations (model) external to this function; if an original network is used, it will be converted to a probability matrix by dividing it by its sum.}
	\item{n}{number of interactions to allocate into the new matrix. Note that in rare cases, the number of interactions assigned may be larger than n due to \option{keep.species=TRUE}.}
	\item{keep.species}{Ensure each species has at least one interaction? Random assignment of interactions may lead to empty columns or rows and hence reduce the dimensions of the simulated web. By default, this is prevented from happening, i.e. each row/column will receive at least one interaction. Setting \option{keep.species} to FALSE may (but need not) cause a loss of species.}
\item{rep.cell}{Whether cells can be assigned interactions repeatedly; defaults to TRUE (quantitative webs); use FALSE for binary webs}
\item{autotransform}{determines how a non-probability web is converted into probabilities; option \option{"sum"} simply divides each entry by the sum of interactions in the web (i.e. probabilities will be proportional to web entries); option \option{"sum_warn"} (the default) does the same, but with a warning;  option \option{"equiprobable"} is an odd name for making probabilities proportional to the product of marginal totals}
}

\details{
This is a generic function to simulate mutualistic networks, based on the original function with the same name used by Vázquez et al. (2009). This function can be used for different purposes, as it allows any type of probability matrix to be used for constructing the simulated matrices. However, this probability matrix must be derived from some other model, either fully synthetic or based on empirical data (if an observed network is given, it will sample from observed interactions, but this is not the typical use case). This function can thus be used for implementing various types of null models, but also highly structured or specialised webs.

The use of an externally generated probability matrix is a key difference to \code{genweb} and functions (methods) in \code{nullmodel}. Nevertheless, \option{autotransfrom="equiprobable"} turns \code{mgen} into a null model function, with interaction probabilities according to marginal totals (row and column sums). This null model does not fix row and column totals (otherwise it is similar to the \code{r2dtable} null model) nor connectance (in difference to \code{vaznull}); also, if \option{keep.species=FALSE}, species may be lost (i.e. have no interactions) by chance.

If \option{rep.cell=TRUE}, repeated interactions are added, thus generating a quantitative matrix with cell values as positive integers. In this case, connectance is not fixed or constrained.
If \option{rep.cell=FALSE}, no repeated assignment of interactions is allowed, thus generating a binary matrix of zeros and ones. Note that when \option{rep.cell=FALSE} the number of interactions to be assigned must be equal or lower than the number of (nonzero) cells in the matrix. If also \option{keep.species=TRUE}, connectance is fixed (\code{n / prod(dim(web))}).
}

\references{
Vázquez, D. P.; Chacoff, N. P. & Cagnolo, L. 2009. Evaluating multiple determinants of the structure of mutualistic networks. \emph{Ecology} \bold{90} 2039--2046
}

\examples{
## Generate simulated matrix from homogeneous probability matrix
probmat <- matrix(1/100, 10, 10)
mgen(web=probmat, n=100)

## Generate binary matrix with probability proportional to degree
## of an observed binary matrix m
obs.mat <- matrix(c(1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0), 5, 5)
rs <- rowSums(obs.mat)
cs <- colSums(obs.mat)
web <- rs \%*\% t(cs)
web <- web/sum(web)
n = sum(obs.mat)
# Allowing zero marginal sums (but there will be none here):
mgen(web, n, keep.species=FALSE, rep.cell=FALSE) 
# Not allowing zero marginal sums:
mgen(web, n, keep.species=TRUE, rep.cell=FALSE) 

## Generate quantitative matrix with probability proportional
## to interaction frequency in an observed matrix m:
# Allowing zero marginal sums:
mgen(mosquin1967, keep.species=FALSE, rep.cell=TRUE) 
# Not allowing zero marginal sums: 
mgen(mosquin1967, keep.species=TRUE, rep.cell=TRUE) 
}
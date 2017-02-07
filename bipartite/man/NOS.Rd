\encoding{UTF-8}
\name{NOS}
\alias{NOS} 


\title{Calculates the node overlap and separation according to Strona \& Veech (2015)}
\description{
This index computes a variation of nestedness, called node overlap and segragation, as well as a modularity measure
}
\usage{
NOS(web, keep.Nij=FALSE, keep.diag=FALSE)
}

\arguments{
  \item{web}{A bipartite interaction web, i.e.~a matrix with higher (cols) and lower (rows) trophic levels.}
  \item{keep.Nij}{Shall node overlap matrix for each trophic level be returned, too? Logical, defaults to FALSE.}
  \item{keep.diag}{Shall the diagonal of Nij be kept at a value of 1, or rather be omitted from computations? If the diagnoal is kept, the index will not be centred on zero as described in the paper. However, the paper is unclear about this point, and hence the option to keep the diagonal (rather than setting omitting it from computation) is offered. Logical, defaults to FALSE.}
}

\details{
According to the authors, NOS is "a new statistical procedure to measure both [the tendency of network nodes to share interaction partners] and the opposite one (i.e. species’ tendency against sharing interacting partners) that we call ‘node segregation’. In addition, our procedure provides also a straightforwardmeasure of modularity, that is, the tendency of a network to be compartmented into separated clusters of interacting nodes."
}

\value{
  \item{Nbar}{the NOS-value (referred to by the authors as a funny N (which with the fc-package in LaTeX is coded as \verb{\m{J}}, apparently).} 
  \item{mod}{modularity, computed as mean standard deviation of the Nij-values for each trophic level.}
  \item{Nbar_higher}{the NOS-value for the higher trophic level.}
  \item{Nbar_lower}{the NOS-value for the lower trophic level.}
  \item{mod_higher}{the mod-value for the higher trophic level.}
  \item{mod_lower}{the mod-value for the lower trophic level.}
  \item{N_ij_higher}{Optional; the matrix of NOS-values for each species pair, for the higher trophic level.}
  \item{N_ij_lower}{Optional; the matrix of NOS-values for each species pair, for the lower trophic level.}
}


\references{ 
Strona, G., and Veech, J.A. (2015) A new measure of ecological network structure based on node overlap and segregation.
\emph{Methods in Ecology & Evolution} \bold{15}, 319--328
}

\author{Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\seealso{\code{\link{grouplevel}}, which eventually shall use this function.}

\examples{
data(Safariland)
# illustrate difference between keeping/removing the diagonal:
NOS(Safariland)
NOS(Safariland, keep.diag=TRUE)
}

\keyword{ package }


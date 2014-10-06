\encoding{UTF-8}

\name{vaznull}
\alias{vaznull}

\title{ Null model with constrained totals and connectance  }
\description{
Implements Diego Vazquez proposal of a null model for pollination networks
}
\usage{
vaznull(N, web)
}
\arguments{
  \item{N}{ Number of desired null model webs. }
  \item{web}{ An interaction matrix.}
}
\details{
 This function produces a null model network with two constraints: a) marginal totals are the same as in the original network (see also r2dtable); b) connectance is the same as in the original network. \code{vaznull} is our implementation of the algorithm propose by Diego Vazquez, hence its name.
 \code{vaznull} differs from \code{\link{swap.web}} both in the algorithm used as well as in the null model it outputs. While \code{vaznull} is slower, we regard it as the better algorithm.
 
 The algorithm was described as follows:
 "The algorithm randomized the total number of individual interactions observed in the original interaction matrix, F. To this end, the algorithm first created a binary matrix, assigning interspecific  interactions  according  to  species-specific probabilities, requiring that each species had at least one interaction. As in Vazquez et al. (2005b), the species-specific probabilities were proportional to species' relative  abundances  (probabilities  are  in  fact  approximately proportional and not equal to relative abundances because of the requirement that each species receives at  least one interaction; this requirement causes probabilities to deviate from relative abundances, especially for rare  species).  Once  the  number  of  filled  cells  in  the original matrix was reached, the remaining interactions were  distributed  among  the  filled  cells,  so  that  connectance in the original and randomized matrices was the same." (Vazquez et al. 2007, page 1122-1123).
 
 Since this leaves a little leeway to the EXACT implementation, check the code for details.
}
\value{
 A list of N randomised matrices with the same dimensions and connectivity as the initial web.
}

\references{
  
 Vázquez, D. P., C. J. Melian, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecological networks. Oikos 116: 1120-1127.
}

\author{ Bernd Gruber <bernd.gruber@canberra.edu.au> & Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }

\note{ 
  It is clearly difficult to decide when a null model is appropriate. It should, in any case, be correctly implemented. Thus, if a null model produces a systematic bias, it should not be used. This seems to be the case for \code{\link{swap.web}}, which yields more high values than necessary. \code{\link{vaznull}} is currently the best alternative.
}

\seealso{ \code{\link{r2dtable}}, \code{\link{swap.web}} }

\examples{

data(Safariland)
networklevel(Safariland, index="info")
networklevel(vaznull(1, Safariland)[[1]], index="info")
system.time(vaznull(10, Safariland))
system.time(swap.web(10, Safariland))
}
\keyword{ package }


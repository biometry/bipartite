\encoding{UTF-8}
\name{r2dexternal}
\alias{r2dexternal}
\alias{vaznullexternal}

\title{Generates null models for network analysis by considering external abundances}
\description{
An extension of r2dtable (and vaznull, respectively) which rescales marginal totals according to independent data
}
\usage{
r2dexternal(N, web, abun.higher=NULL, abun.lower=NULL)

vaznullexternal(N, web, abun.higher=NULL, abun.lower=NULL)
}

\arguments{
  \item{N}{number of null models to be generated.}
  \item{web}{Web is a matrix representing the interactions observed between higher level species (columns) and lower level species (rows). }
  \item{abun.higher}{Optional vector of externally measured abundances of the higher level. If missing (NULL) it will be replaced by column totals.}
  \item{abun.lower}{Optional vector of externally measured abundances of the lower level. If missing (NULL) it will be replaced by row totals.}
}

\details{
The underlying functions are \code{r2dtable} and \code{\link{vaznull}}, which require a vector of row and column totals or a web, respectively. In function \code{\link{nullmodel}}, these marginal totals are computed from the observed interaction matrix. Here, external abundances can be provided. These will be rescaled and combined with the observed marginal total to construct new row and column vectors for \code{r2dtable}.
 
If neither row nor column abundances are provided this function will be identical to \code{r2dtable} and \code{\link{vaznull}}, respectively. 
}

\value{
  Returns a list of \code{N} null model-generated networks. Species names are (obviously) dropped.
}

\note{
Since the function contains a rounding operation, it also has to include a re-distribution of 1s to keep all species in the system. That means, if one species has been observed in the external abundances extremely rarely (compared to the others), it will be overrepresented in this null model, because otherwise it would have to be dropped altogether! If you have a better solution, please let me know.

Also note that \code{\link{vaznull}} itself does not return a matrix with exactly the same marginal totals as the input, so don't expect any different from \code{vaznullexternal}!
}


\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\seealso{ \code{\link{nullmodel}}, \code{\link{vaznull}}, \code{\link{mgen}}, \code{vegan::simulate} and \code{r2dtable}
}

\examples{
\dontrun{
	abun.lower <- c(15,5,2,7,4,8,6,0.01,6)
	 set.seed(2)
	(abun.higher <- rpois(27, lambda=4))
	abun.higher[1] <- 0.001
	sum(ext.polls)
	## note: external abundances do not sum up; this is intentional!!
	r2dexternal(2, Safariland, abun.higher=abun.higher, abun.lower=abun.lower)
	r2dexternal(2, Safariland, abun.higher=abun.higher)

  vaznullexternal(2, Safariland, abun.higher=abun.higher, abun.lower=abun.lower)
	}
}

\keyword{ package }



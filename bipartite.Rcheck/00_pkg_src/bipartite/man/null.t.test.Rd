\name{null.t.test}

\alias{null.t.test}

\title{Compares observed pattern to random webs.}

\description{
A little null-model function to check, if the observed values actually are much different to what one would expect under random numbers given the observed row and column totals (i.e.~information in the structure of the web, not only in its species' abundances). Random matrices are based on the function \code{\link{r2dtable}}. The test itself is a t-test (with all its assumptions).
}

\usage{
null.t.test(web, N = 30, ...)
}

\arguments{
  \item{web}{A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows).}
  \item{N}{Number of null models to be produced; see \sQuote{Note} below!}
  \item{\dots}{Optional parameters to be passed on to the functions
      \code{\link{networklevel}} and \code{\link{t.test}}.}
}
\details{
  This is only a very rough null-model test. There are various reasons why one may consider \code{r2dtable} as an incorrect way to construct null models (e.g.~because it yields very different connectance values compared to the original). It is merely used here to indicate into which direction a proper development of null models may start off. Also, if the distribution of null models is very skewed, a t-test is obviously not the test of choice.
  
 Finally, not all indices will be reasonably testable (e.g.~number of species is fixed), or are returned by the function \code{\link{networklevel}} in a form that \code{null.t.test} can make use of (e.g.~degree distribution fits).
}
\value{
   Returns a table with one row per index, and columns giving
   \item{obs}{observed value}
   \item{null mean}{mean null model value}
   \item{lower CI}{lower 95\% confidence interval (or whatever level is specified in the function's call)}
   \item{upper CI}{upper 95\% confidence interval (or whatever level is specified in the function's call)}
   \item{t}{t-statistic}
   \item{P}{P-value of t statistic}
}

\author{Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\note{This function is rather slow. Using large replications in combination with iterative indices (degree distribution, compartment diversity, extinction slope, H2) may lead to rather long runtimes!
}

\examples{
data(mosquin1967)
null.t.test(mosquin1967, index=c("generality", "vulnerability",
    "cluster coefficient", "H2", "ISA", "SA"), nrep=2, N=10)
}
\keyword{package}


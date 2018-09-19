\encoding{UTF-8}
\name{H2fun}

\alias{H2fun}

\title{ Specialisation of a bipartite web. }

\description{
  Calculates the overall level of specialisation of all interacting species in a bipartite web.
}
\usage{
H2fun(web, H2_integer=TRUE)
}

\arguments{
  \item{web}{ A matrix representing the interactions observed between higher trophic  level species (columns) and lower trophic level species (rows). Usually this will  be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{H2_integer}{logical; indicates whether web entries are integer numbers. If set to false, H2fun can be used also on rates, percentages and other non-integer values indicating the intensity of an interaction.}
}

\details{
  H2' is an index describing the level of \dQuote{complementarity specialisation} (or should one say: selectiveness?) of an entire bipartite network (Blüthgen et al. 2006). It describes to which extent observed interactions deviate from those that would be expected given the species marginal totals. The more selective a species, the larger is H2' for the web. 
  
  H2' is an extension of d' (see \code{\link{dfun}}) for the entire network. %Check there for details on finding maximum and minimum H2 values for a given web.
  
For non-integer values, H2 max can be readily computed and is thus more reliable and much faster.
 
}

\value{
  \item{H2}{The H2'-value for the web matrix.}
  \item{H2min}{Heuristic minimum H2-value for the web matrix.}
  \item{H2max}{Heuristic maximum H2-value for the web matrix.}
  \item{H2uncorr}{Uncorrected H2-values (before ranging between min and max), rounded  to three digits.}
}

\references{ 
Blüthgen, N., Menzel, F. and Blüthgen, N. (2006) Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 9.
 }

\author{ Carsten F. Dormann and Jochen Fründ }

\seealso{ \code{\link{dfun}} following the same idea for each species in the web matrix. }

\examples{
data(Safariland)
H2fun(Safariland)
}

\keyword{ package }


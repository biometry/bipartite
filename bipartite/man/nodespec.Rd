\encoding{UTF-8}

\name{nodespec}

\alias{nodespec}

\title{Calculates the node-based specialisation index }

\description{
 Calculates a specialisation index based on the node positions for all species 
 in a bipartite network, separately for the higher and lower trophic level.
}

\usage{
nodespec(web, inf.replace = NA)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic level species as columns and number of interactions as entries.}
  \item{inf.replace}{What should infinite geodesic distances (e.g. between compartments)  be represented as? Defaults to \option{NA}; only currently implemented alternative is \option{inf.replace=Inf}, which replaces infinite distances by the maximum path length plus 1.}
}

\details{
 This index aims to describe the functional specialisation of pollinators and was proposed by Dalgaard et al. (2008). It is a purely qualitative measure. 
 
 After calculating the geodesic distances between species, i.e. the minimum number of steps from one species to another, these values are averaged for each species. 
This mean geodesic distance is interpreted as functional specialisation (Dalgaard et al. 2008).
 
 Notice that this ``new'' index is in fact little else than the inverse of (unscaled) closeness centrality in (poor) disguise.
}

\value{
  A list with two components, names \dQuote{higher} and \dQuote{lower}, both containing the node specialisation index for each species.
}

\references{ 

Dalsgaard, B., Martín González, A. M., Olesen, J. M., Timmermann, A., Andersen, L. H. and Ollerton, J. (2008) Pollination networks and functional specialization: a test using Lesser Antillean plant-hummingbird assemblages. \emph{Oikos} \bold{117}, 789--793

 }

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\note{
This index is as yet unevaluated. We don't know how it responds to true specialisation at all. In fact, it is a rather good example of how to get a new thing published 
without even having demonstrated in which way it differs from existing indices of specialisation (such as standardised d included in the function \code{\link{dfun}}), 
or how it performs on artificial data with known properties.

One major disadvantage of any index based on path lenghts is its difficulty with compartments, i.e. species not linked to the rest of the network. There are, generally 
speaking, three ways to handle this: Firstly, ignore it (that is, set infinite distances to \code{NA}; our default). Secondly, leave it as it is (that is, leave infinite distances 
as infinite). This is not really an option, since then ALL species would have infinite specialisation values. Thirdly, replace infinite by the largest distance plus one (see 
comments in \code{geodist} in \pkg{sna}). That would probably be a plausible thing to do, since we could argue that with a little bit extra observation we might have 
found a species linking a compartment to the rest of the network. However, this solution is \dQuote{not canonical}, as put in \code{geodist} and hence biased to an 
unknown extent. To use this option, specify \option{inf.replace=Inf}.

}

\seealso{ See also as \code{\link{specieslevel}}, which calls \code{nodespec}. }

\examples{
data(Safariland)
nodespec(Safariland, inf.replace=Inf)
}

\keyword{package}

\encoding{UTF-8}
\name{massaction}

\alias{massaction}

\title{Adjusts a network matrix with regard to mass action}

\description{
  Implementation of mass action calculations based on 
  section 6.2 "Effective abundance, pseudo-inverse and mass action" in
  the Supplementary Information of 
  "The ghost of nestedness in ecological networks" 
  by Staniczenko et al. (2013).
}

\usage{
  massaction(web)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic 
             level species as columns and number of interactions as entries.}
}

\details{
  
}

\value{
  The network matrix adjusted for mass action effects.
}

\references{
  Staniczenko, P., Kopp, J. & Allesina, S. (2013) 
  The ghost of nestedness in ecological networks. 
  \emph{Nat. Commun.}  \bold{4}, 1391. \doi{10.1038/ncomms2422}
}

\author{ Tobias Bauer \email{tobias.bauer-2@studium.uni-hamburg.de}}

\examples{
  data(bezerra2009, package="bipartite")
  mass.action.norm.bezerra2009 <- massaction(bezerra2009)
}

\keyword{ package }


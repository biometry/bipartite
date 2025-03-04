\encoding{UTF-8}
\name{spectral.radius}

\alias{spectral.radius}
\alias{massaction}

\title{Calculates the spectral radius as a nestedness measurement of a bipartite matrix, possibly with regard to mass action}

\description{
  This function calculates the spectral radius (i.e. the maximum absolute value of the Eigenvalues) 
  of the square adjacency matrix corresponding tp a given bipartite matrix. This may be used as a measure of 
  nestedness as described by Staniczenko et al. (2013).
}

\usage{
  spectral.radius(web, mass.action.norm = FALSE)
  massaction(web)
}

\arguments{
  \item{web}{A matrix with lower trophic level species as rows, higher trophic 
    level species as columns and number of interactions as entries.}
  \item{mass.action.norm}{If TRUE the matrix is normalized for mass action 
    using the massaction() function before calculating 
    the spectral radius.}
}

\section{Massaction}{
  Implementation of mass action calculations based on 
  section 6.2 "Effective abundance, pseudo-inverse and mass action" in
  the Supplementary Information of Staniczenko et al. (2013).
}

\value{
  The spectral radius (i.e. the largest Eigenvalue) of the square adjacency matrix
  corresponding to the given bipartite matrix.
}

\references{
  Staniczenko, P., Kopp, J. & Allesina, S. (2013) The ghost of nestedness in ecological networks. \emph{Nat. Commun.} \bold{4}, 1391. \doi{10.1038/ncomms2422}
}

\author{ Tobias Bauer \email{tobias.bauer-2@studium.uni-hamburg.de}}

\examples{
  data(bezerra2009, package="bipartite")
  spectral.radius(bezerra2009, mass.action.norm = TRUE)
}

\keyword{ package }
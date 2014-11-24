\name{extinction}

\alias{extinction}

\title{ Simulates extinction of a species from a bipartite network }

\description{
  Following (how I remember) the paper of Memmott et al. (2004), this function deletes a column (e.g. pollinator) or row (e.g. plant). Only a helper function for \code{\link{second.extinct}}, really.
}

\usage{
extinction(web, participant = "both", method = "random", ext.row=NULL, 
	ext.col=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{web}{ A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{participant}{Which level of participant to remove: \option{lower} removes a row, \option{higher} removes a row, \option{both} randomly picks either row or column.}
  \item{method}{Determines sequence of extinctions: \option{random} removes a random participant, while \option{abundance} removes the least abundant species first. \option{external} will use the externally provided vector to determine extinction sequence.}
   \item{ext.row}{Optional vector giving the sequence in which lower-level species are to be deleted.}
   \item{ext.col}{Optional vector giving the sequence in which higher-level species are to be deleted.}
}

\details{
  In itself rather useless. Called repeatedly by \code{\link{second.extinct}} to build an extinction sequence and accordingly a sequence of secondary extinctions.
}

\value{
Returns the same matrix that was given as object, just with one row or column being turned into zeros.
}
\references{ 
Memmott, J., Waser, N. M. and Price, M. V. 2004 Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

 }

\author{ Carsten F. Dormann }

\seealso{ \code{\link{second.extinct}} }

\examples{
\dontrun{
	data(Safariland)
	(w <- extinction(Safariland, participant="lower", method="abun"))
	empty(w, count=TRUE)
}
}

\keyword{ package }


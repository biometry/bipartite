\encoding{UTF-8}
\name{as.bipartite_web}
\alias{as.bipartite_web}
\alias{as.matrix.bipartite_web}

\title{Construct bipartite web and bipartite web array objects}
\description{
\code{bipartite_web()} creates an S3 object of class \code{"bipartite_web"} from a bipartite adjacency matrix.
\code{bipartite_webarray()} creates an S3 object of class \code{"bipartite_webarray"} from one or more bipartite webs,
combining them into a three-dimensional array.

The \code{"bipartite_web"} class represents a single bipartite network with node-level attributes,
while \code{"bipartite_webarray"} represents a collection of such webs.
}

\usage{
as.bipartite_web(x, ...)

\method{as.matrix}{bipartite_web}(x, ...)

}

\arguments{
    \item{x}{A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
}

\value{
\code{bipartite_web()} returns an object of class \code{"bipartite_web"}.

\code{bipartite_webarray()} returns an object of class \code{"bipartite_webarray"}.
}

\examples{
data(Safariland, vazquenc, vazquec)
Safariland.web <- bipartite_web(Safariland)
print(Safariland.web)
}

\keyword{ package }
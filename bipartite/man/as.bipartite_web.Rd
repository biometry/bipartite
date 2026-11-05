\encoding{UTF-8}
\name{as.bipartite_web}
\alias{as.bipartite_web}

\title{Coerce to a bipartite_web}
\description{
    Function to coerce matrices, data.frames or similar objects to \code{bipartite_web} objects. 
}

\usage{
as.bipartite_web(x)
}

\arguments{
    \item{x}{
        A matrix, data.frame or other object, that can be coerced into a matrix.
    }
}

\value{
    Returns a \code{bipartite_web} object.
}

\examples{
    m <- matrix(1:12, ncol = 3)
    as.bipartite_web(m)
    df <- data.frame(insect1 = 1:4, insect2 = 4:1)
    as.bipartite_web(df)
}

\keyword{ package }
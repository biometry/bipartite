\name{empty}

\alias{empty}

\title{Deletes empty rows and columns from a matrix.}

\description{
Gets rid of empty columns and rows in a matrix. Optionally counts removed rows and columns, and returns these values as attribute.
}

\usage{
empty(web, count=FALSE)
}

\arguments{
  \item{web}{A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{count}{Logical. Shall be counted how many columns and rows were removed? Numbers are returned in attribute. Defaults to FALSE.}
}

\details{
  Helper function to remove empty (i.e. all-zero or all-NA) rows and columns, thereby concentrating the matrix. This function is also invoked for its side effect by \code{\link{extinction}} to investigate the effect of removing a species from the network.
}

\value{
  Returns matrix without empty rows or columns. Its attribute `out' (if \code{count=TRUE}) contains a named vector with the number of rows removed and the number of columns removed.
}

\author{ Carsten F. Dormann }

\seealso{ \code{\link{extinction}} and \code{\link{second.extinct}}, which repeatedly call \code{empty}. }

\examples{
data(Safariland)
web <- Safariland
web[,5] <- 0
empty(web, count=TRUE)
attr(empty(web), "empty")

}

\keyword{ package }


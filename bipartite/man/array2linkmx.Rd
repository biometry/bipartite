\encoding{UTF-8}
\name{array2linkmx}
\alias{array2linkmx}

\title{Reshape a webarray to a web X link matrix}

\description{
Function to turn an array with sites as third dimension into a web by link matrix (e.g. sites X links). This is the "link community matrix" to use for dissimilarity calculations. Mostly just a helper function for \code{betalinkr}. 
}

\usage{
array2linkmx(webarray)
}


\arguments{
  \item{webarray}{An array of two or more networks. Assumes the third dimension is the webID in the array, which will become the first dimension in the link-matrix (output).}
}

\details{
This function converts the two-dimensional adjacency matrices (i.e., bipartite webs) to one-dimensional vectors (similar to what network people call ``edgelists''). These vectors become rows of a matrix (which has one row per web). This makes the data available to community ecology methods, e.g. those offered by the \code{vegan} package. Links are treated equivalently to species in a ``normal'' community analysis (note that it makes no difference whether one or both partner of an interaction differ, in both cases the link has a different identity). The function is used here mostly as a helper for interaction dissimilarity calculations with \code{betalinkr}.

\code{dimnames} are optional but recommended for the webarray. Names of the third dimension will become rownames in the output. Colnames (i.e. names of links) will be created from dimnames[[1]] (``lower species'') and dimnames[[2]] (``higher species''), separated by ``__'' (double underscore).
}

\value{
A matrix with one row per web and one column per link (i.e. per combination of lower and higher species). All-zero columns may often be included.}

%\note{
%}

%\references{ 
% }

\author{ Jochen Fründ}


\examples{
array2linkmx(webs2array(Safariland, vazquenc))
}
\keyword{ package }


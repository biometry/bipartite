\name{listModuleInformation}

\alias{listModuleInformation}

\title{listModuleInformation}

\description{
This function takes an object of class "moduleWeb" and returns information about the names of the nodes of which the computed modules exist.
}

\usage{
listModuleInformation(moduleWebObject)
}

\arguments{
  \item{moduleWebObject}{Object of class "moduleWeb"}
}

\value{
The value of the function is a list of lists of lists of vectors representing the names of the nodes involved in a certain module.

\item{value[[x]][[y]][[1]]}{vector with the names of the "row nodes" comprised by the z.th module in depth x-1 (sic!)}
\item{value[[x]][[y]][[2]]}{vector with the names of the "column nodes" comprised by the z.th module in depth x-1 (sic!)}
}

\author{Rouven Strauss}


\examples{
\dontrun{
data(small1976)

moduleWebObject = computeModules(small1976);
moduleList = listModuleInformation(moduleWebObject);
}
}
\keyword{Methods and Generic Functions}
\keyword{Clustering}

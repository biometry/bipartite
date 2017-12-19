\encoding{UTF-8}

\name{printoutModuleInformation}

\Rdversion{1.1}

\alias{printoutModuleInformation}

\title{printoutModuleInformation}

\description{
This takes an object of the class "moduleWeb" and prints out the information about the modules stored in this object.
It is a formatted print-out of the information one gets as result of the function listModuleInformation(moduleWebObject).
}

\usage{
printoutModuleInformation(moduleWebObject)
}

\arguments{
  \item{moduleWebObject}{An object of class "moduleWeb".}
}

\value{None. This function is called for its side effects of printing the content of its object in a more acceptable format.}

\author{Rouven Strauss}

\examples{
\dontrun{
data(small1976)
moduleWebObject = computeModules(small1976)
printoutModuleInformation(moduleWebObject)
}
}

\keyword{Methods and Generic Functions}

\keyword{Clustering}

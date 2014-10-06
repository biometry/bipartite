\encoding{UTF-8}
\name{moduleWeb-class}
\docType{class}

\alias{moduleWeb-class}

\title{Class "moduleWeb"}

\description{
This class is the output of an application of the function computeModules(...) to a graph. It consists of the matrix representing the original graph which has been passed to computeModules(...) in order to compute modules, a matrix representing the same graph but permutated according to the identified modules, two vectors indicating the permutation of row and column indices, respectively, and information about the modules themselves.
}
\section{Objects from the class}{
Objects from the class should only be created by using the function computeModules(...).
}
\section{Slots}{
	 \describe{
    \item{\code{likelihood}:}{Contains a number with the likelihood-equivalent of the final proposed module structure. This value is the same value as Q (or M), the modularity as given by Newman or Guimerà & Amaral (2005). }	 
    \item{\code{originalWeb}:}{Object of class \code{"matrix"} representing the original bipartite graph in which modules have been computed.}
    \item{\code{moduleWeb}:}{Object of class \code{"matrix"} representing the original bipartite graph but reordered such that plotting modules is possible.}
    \item{\code{orderA}:}{Object of class \code{"vector"} representing the permutation of the rows of the original graph.}
    \item{\code{orderB}:}{Object of class \code{"vector"} representing the permutation of the columns of the original graph.}
    \item{\code{modules}:}{Object of class \code{"matrix"} containing for each module the information about its depth and involved nodes.}
  }
}
\section{Methods}{
Objects of this class are used in following functions:

listModuleInformation(moduleWebObject)

printoutModuleInformation(moduleWebObject)

plotModuleWeb(moduleWebObject, plotModules=TRUE, rank=FALSE, weighted=TRUE, displayAlabels=TRUE, displayBlabels=TRUE, labsize=1, plotsize=12, xlabel="", ylabel="", square.border="white", fromDepth=0, upToDepth=-1)
}
\author{Rouven Strauss}
\examples{
showClass("moduleWeb")
}
\keyword{classes}
\keyword{modules}
\keyword{moduleWeb}
\keyword{modularity}
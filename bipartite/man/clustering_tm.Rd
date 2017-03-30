\name{clustering_tm}
\alias{clustering_tm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Redefined clustering coefficient for two-mode networks }
\description{ This function calculates the two-mode clustering coefficient as proposed by Opsahl (2010). \cr
%Note: If you are having problems with this function (i.e., run out of memory or it being slow for simulations), there is a quicker and much more memory efficient c++ function. However, this function is not fully integrated in R, and requires a few extra steps. Send me an email to get the source-code and Windows-compiled files.  
}
\usage{clustering_tm(net, subsample=1, seed=NULL)}
\arguments{
  \item{net}{A binary or weighted two-mode edgelist }
  \item{subsample}{Whether a only a subset of 4-paths should we used when calculating the measure. This is particularly useful when running out of memory analysing large networks. If it is set to 1, all the 4-paths are analysed. If it set to a value below one, this is roughly the proportion of 4-paths that will be analysed. If it is set to an interger greater than 1, this number of ties that form the first part of a 4-path that will be analysed. Note: The c++ functions are better as they analyse the full network. }
  \item{seed}{If a subset of 4-paths is analysed, by setting this parameter, the results are reproducable.}
}
\value{Returns the outcome of the equation presented in the paper}

\references{Opsahl, T. 2010. Triadic closure in two-mode networks: Redefining the global and local clustering coefficients. \emph{arXiv},1006.0887  }

\author{ Tore Opsahl; http://toreopsahl.com }

\note{ version 1.0.0, taken, with permission, from package tnet  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
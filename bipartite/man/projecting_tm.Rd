\name{projecting_tm}

\alias{projecting_tm}

\title{ Projecting binary and weighted two-mode networks onto weighted one-mode networks. }
\description{ This function is the implemtation of the procedure outlined on \cr
https://toreopsahl.com/2009/05/01/projecting-two-mode-networks-onto-weighted-one-mode-networks/ }
\usage{projecting_tm(net, method = "sum")}

\arguments{
  \item{net}{A two-mode edgelist}
  \item{method}{ The method-switch control the method used to calculate the weights.\cr
  binary sets all weights to 1\cr
  sum sets the weights to the number of cooccurances\cr
  Newman bases the weights on Newman's (2001) method of discounting for the size of collaborations. }
}

\value{Returns a one-mode network}
\references{Opsahl. T. 2009. Projecting two-mode networks onto weighted one-mode networks. Available at: https://toreopsahl.com/2009/05/01/projecting-two-mode-networks-onto-weighted-one-mode-networks/ }

\author{ Tore Opsahl; https://toreopsahl.com }

\note{ version 1.0.0, taken, with permission, from package tnet }

\examples{
#please download and investigate tnet for examples!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
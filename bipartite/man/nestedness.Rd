\encoding{UTF-8}
\name{nestedness}
\alias{nestedness}

\title{Calculates nestedness temperature of presence/absence matrices}

\description{Deprecated: Calculates matrix temperature using the binmatnest approach of Miguel Rodríguez-Gironés}

\usage{
nestedness(m, null.models = TRUE, n.nulls = 100, popsize = 30, 
n.ind = 7, n.gen = 2000, binmatnestout=FALSE)
}

\arguments{
  \item{m}{ \code{m} is the matrix object for which the temperature is calculated. \code{m} will be converted to a binary matrix as temperature is only based on binary data}
  \item{null.models}{Logical; shall the three different null models to check for significance of the matrix temperature be calculated? The null models procedure is 
quite time consuming and therefore we added this switch. Defaults to \code{null.models}=TRUE.}
  \item{n.nulls}{How many null models should be calculated. Defaults to \code{n.nulls=100}.}
  \item{popsize}{For the genetic algorithm some parameters have to be initialised. First is \code{popsize}, default is 30}
  \item{n.ind}{Second is number of individuals picked for the next generation. Default of \code{n.ind} is 7.}
  \item{n.gen}{Third is the number of generations until the genetic algorithm stops. Default of \code{n.gen} is 2000.}
  \item{binmatnestout}{if set to TRUE a file "binmat.out" is saved in the current working directory , which stores the original binmatnest output}
}


\details{
This is a note that the original function provided here, and detailed below, has been retired. \pkg{vegan}'s \code{\link[vegan]{nestedtemp}} is more stable and replaces this call whereever we use it in \pkg{bipartite}. Indeed, calling \pkg{bipartite}'s \code{nestedness} now simply call's \pkg{vegan}'s \code{\link[vegan]{nestedtemp}}. The function will eventually disappear, but for now this fix is provided as legacy support for dependent packages.

All arguments except first are ignored.

There are several implementations of nestedness-calculators, most noticeably NTC (nestedness temperature calculator), BINMATNEST and aninhado (check Wikipedia's entry on the subject: \url{https://en.wikipedia.org/wiki/Nestedness}). While we used BINMATNEST, this does not disqualify any of the others. Miguel was simply the first we contacted and he was readily willing to share his code. 

We used BINMATNEST by calling a tweaked version of the C++ program binmatnest. In principle nestedness temperature is calculated by using a line of perfect order (using a genetic algorithm) to determine the reordering of rows and columns that leads to minimum matrix temperature of given size and fills. The deviation from this minimun temperature is the matrix temperature. In addition nestedness uses different null models to check for statistical significance of the matrix temperature. For details on what BINMATNEST does different, and better, than the original NTC see reference below. 

Notice also that the original software BINMATNEST is available as a stand-alone application. Check out Miguel's homepage: \url{http://www.eeza.csic.es/eeza/personales/rgirones.aspx}
% or download directly: \url{http://www.eeza.csic.es/eeza/personales/rgirones/File/BINMATNEST3.zip}.
}


\references{
Rodríguez-Gironés M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935
}



\keyword{ package}
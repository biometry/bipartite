\encoding{UTF-8}
\name{nestedness}
\alias{nestedness}

\title{Calculates nestedness temperature of presence/absence matrices}

\description{Calculates matrix temperature using the binmatnest approach of Miguel Rodríguez-Gironés}

\details{
This is merely a note that the original function provided here, and detailed below, has been retired. \pkg{vegan}'s \code{\link[vegan]{nestedtemp}} is more stable and replaces this call whereever we use it in \pkg{bipartite}.

There are several implementations of nestedness-calculators, most noticeably NTC (nestedness temperature calculator), BINMATNEST and aninhado (check Wikipedia's entry on the subject: \url{https://en.wikipedia.org/wiki/Nestedness}). While we used BINMATNEST, this does not disqualify any of the others. Miguel was simply the first we contacted and he was readily willing to share his code. 

We used BINMATNEST by calling a tweaked version of the C++ program binmatnest. In principle nestedness temperature is calculated by using a line of perfect order (using a genetic algorithm) to determine the reordering of rows and columns that leads to minimum matrix temperature of given size and fills. The deviation from this minimun temperature is the matrix temperature. In addition nestedness uses different null models to check for statistical significance of the matrix temperature. For details on what BINMATNEST does different, and better, than the original NTC see reference below. 

Notice also that the original software BINMATNEST is available as a stand-alone application. Check out Miguel's homepage: \url{http://www.eeza.csic.es/eeza/personales/rgirones.aspx}
% or download directly: \url{http://www.eeza.csic.es/eeza/personales/rgirones/File/BINMATNEST3.zip}.
}


\references{
Rodríguez-Gironés M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935
}



\keyword{ package}
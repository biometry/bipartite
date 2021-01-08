\encoding{UTF-8}
\name{nestedness}
\alias{nestedness}

\title{Calculates nestedness temperature of presence/absence matrices}

\description{Calculates matrix temperature using the binmatnest programm of Miguel Rodríguez-Gironés}

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
There are several implementations of nestedness-calculators, most noticeably NTC (nestedness temperature calculator), BINMATNEST and aninhado (check Wikipedia's entry on the subject: \url{https://en.wikipedia.org/wiki/Nestedness}). While we here use BINMATNEST, this does not disqualify any of the others. Miguel was simply the first we contacted and he was readily willing to share his code. 

Here we use BINMATNEST by calling a tweaked version of the C++ program binmatnest. In principle nestedness temperature is calculated by using a line of perfect order (using a genetic algorithm) to determine the reordering of rows and columns that leads to minimum matrix temperature of given size and fills. The deviation from this minimun temperature is the matrix temperature. In addition nestedness uses different null models to check for statistical significance of the matrix temperature. For details on what BINMATNEST does different, and better, than the original NTC see reference below. 

Notice also that the original software BINMATNEST is available as a stand-alone application, too. Check out Miguel's homepage: \url{http://www.eeza.csic.es/eeza/personales/rgirones.aspx}
% or download directly: \url{http://www.eeza.csic.es/eeza/personales/rgirones/File/BINMATNEST3.zip}.
}

\value{
  returns a list of matrix descriptors, such as
  \item{temperature }{the matrix temperature}
  \item{parameters of genetic algorithms}{Parameters used for the genetic algorithm}
  \item{nullmodels}{switch if null models have been calculated, 1 for yes, 0 for no}
  \item{p, mean, var}{probability, mean temperature and variance of temperature for the three different null models}
  \item{packing order}{the packing order of the most packed matrix (minimum temperature of a perfectly nested matrix using given size and fills.}
}

\references{
Rodríguez-Gironés M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935
}

\author{Bernd Gruber, based on C++-code by Miguel Rodriguez-Girones.}

\note{ Make sure matrix \code{m} is valid, as error proofing in the C++ function does not fully work, and therefore it is possible that R may crash when using 
strange types of matrices, such as matrices with only one entry.

Also, while this function returns exactly the same values as binmatnest, there seems to be the occasional glitch in the software, in that in inexplicable 
circumstances the matrix is inverted before calculating nestedness. We are currently trying to find the cause of this problem. Sorry!

In the meantime, use \code{nestedtemp} in \pkg{vegan} (as we do in \code{networklevel})!
}

\examples{
	\dontrun{
	data(vazarr)
	nestedness(vazarr) # null models are calculated
	# no null models, much faster for bigger matrices:
	nestedness(vazarr, null.models=FALSE) 
	nestedness(vazarr, n.nulls=300, n.gen=300)
	}
}

\keyword{ package}
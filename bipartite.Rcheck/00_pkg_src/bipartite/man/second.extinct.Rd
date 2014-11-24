\name{second.extinct}

\alias{second.extinct}

\title{Secondary extinctions in bipartite networks}

\description{
Calculates the consequences of removing a species from a bipartite network. With plotting and slope-estimation functionality.
}

\usage{
second.extinct(web, participant = "higher", method = "abun", nrep = 10,
  details = FALSE, ext.row=NULL, ext.col=NULL)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{participant}{\code{high} (default) or \code{low} or \code{both}, depending if you want to exterminate higher trophic level species, lower trophic level species, or randomly choose species of both levels.}
  \item{method}{Random deletion of a species (\option{random}); according to its abundance, with least abundant going extinct first (\option{abundance}; default) or \option{"degree"} to build a sequence from the best-to-least connected species. This is the most extreme case, where the most generalist species goes extinct first (see Memmott et al. 1998). (Note that \option{method="abundance"} does not work with \option{participant="both"}; in this case a random sequence of species from both trophic levels is generated.) \option{external} will use the externally provided vector to determine extinction sequence.}
  \item{nrep}{Number of replicates of extermination sequence. Will not be considered for method abundance.}
  \item{details}{Logical; returns details, i.e. for each replicate the sequence of secondary extinctions. If set to \code{FALSE} (default), replicated runs will be averaged.  Using anything but \code{FALSE} will not yield the correct object for further analysis. This switch is really only for understanding what \code{second.extinct} does and for possible debugging.}
  \item{ext.row}{Optional vector giving the sequence in which lower-level species are to be deleted.}
  \item{ext.col}{Optional vector giving the sequence in which higher-level species are to be deleted.}
}
\details{
The procedure of this function is simple. For example imagine the web to represent a pollination web, in which pollinators die one by one. Set all entries of a  column to zero, see how may rows are now also all-zero (i.e. species that are now not pollinated any more), and count these as secondary extinctions of the primary exterminated pollinator. 

Internally, each extermination is achieved by a call to \code{\link{extinction}}, followed by a call to \code{\link{empty}}, which counts the number of all-zero columns and rows.

Although written for pollination webs (hence the nomenclature of pollinator and plant), it can be similarly applied to other types of bipartite networks. It is called by \code{\link{networklevel}}.
%
%\option{details=FALSE} will not work with \option{participant="both"}, since sequences cannot be averaged, due to different lengths.
}

\value{
Function returns an object of class \dQuote{bipartite} to ensure proper working of function \code{plot.bipartite}.

If \option{details=FALSE}, the returned object contains a matrix with columns representing the number of species going extinct from one step to the next, averaged across all runs.

If \option{details=FALSE}, the returned object will be a list of matrices containing the number of species going extinct at each step.

The objects attribute \dQuote{exterminated} gives the name of the trophic level (pollinator or plant).
}

%\note{You should not combine \option{method="random"} with \option{participant="both"}! The reason is that extinction sequences are of different length in this case and slope.bipartite cannot handle that. Thus, while second.extinct will do the computations alright, you cannot use \code{slope.bipartite} to compute (and plot) slopes from them. Please use \option{method="random"} for each level separately when followed by \code{slope.bipartite}!}

\references{
Memmott, J., Waser, N. M. and Price, M. V. (2004) Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611.
 }

\author{ Carsten F. Dormann }

\note{
Note: The length of an extinction sequence is obviously given by the number of species in the selected trophic level. When setting \option{participant="both"}, lengths will differ for each replicate run, since it is unpredictable in which sequence species go extinct, and hence how many secondary extinctions will pre-empt further primary extinctions.
}

\seealso{ \code{\link{networklevel}} calls \code{second.extinct}; \code{\link{extinction}}
and \code{\link{empty}} are internal helper functions, and \code{\link{slope.bipartite}}
calculates extinction slopes from the output of \code{second.extinct}.}

\examples{
\dontrun{
data(Safariland)
(ex <- second.extinct(Safariland, participant="lower", method="random", nrep=50, 
	details=TRUE))
(ex <- second.extinct(Safariland, participant="lower", method="random", nrep=50, 
	details=FALSE))
}
}

\keyword{package}


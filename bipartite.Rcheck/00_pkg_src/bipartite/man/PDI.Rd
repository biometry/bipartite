\encoding{UTF-8}
\name{PDI}
\alias{PDI}

\title{Paired Differences Index}
\description{ Computes the Paired Differences Index }

\usage{
PDI(web, normalise=TRUE, log=FALSE)
}

\arguments{
  \item{web}{A bipartite interaction web, i.e.~a matrix with higher (cols) and lower (rows) trophic levels.}
  %JF: some changes here:
  \item{normalise}{Logical; divides for each species by the maximum of interactions. Thereby species can be compared among each other and values range between 0 and 1. Defaults to TRUE, which differs from Poisot et al. (2011a). Note that Tim Poisot also recommends the normalised computation as default.}
  \item{log}{logical; since number of interactions is often highly skewed, the log yields a more even spread of PDI-values across species. Defaults to \code{FALSE}.}
}

\details{
There are many ways to skin a cat. This is a more recent addition proposed by Poisot et al. (2011a) and used for mutualistic network analysis by Poisot et al. (2011b). This function can be used alone or through \code{\link{specieslevel}}. 

% JF: some changes here
If P is interaction strength (typically interaction frequency or proportion depending on argument \code{normalized}), then PDI for this species is computed as:

 sum(P1 - Pi) /(H-1), 
 
where P1 is the highest number of interactions in a link, while Pi are the remaining values. H is the number of potential interactors (e.g. plant species if the target species is a pollinator). 
}

\note{
\enumerate{
% JF: changes: (a) as Schoener doesn't use resource range, just generalization, the cite should refer to Poisot2012; (b) corrected that it is not exactly the proportion of unused resources [e.g. 100% unused would mean no interactions, but here 1 means one interaction] 
	\item When a binary web is fed to this function (e.g. \code{PDI(web>0)}) this function returns Poisot et al.'s (2012) \dQuote{resource range}. Resource range has a value of 0 when all resources are used, and a value of 1 when only one resource is used. It is thus more an \dQuote{unused resource range}.
% Actually, PDI builds on this idea (see Poisot et al. 2012 for further details).

%JF: changed completely
  \item This index was originally proposed for performance estimates with standardized resource frequency. It aims to characterize the performance decay as an organism moves away from its optimal resource. As such, it is strongly influenced by the performance on the optimal resource (by its absolute value with option \code{normalised=FALSE} and by its relative value in relation to the sum of all performances with \code{normalised=TRUE}). It is less clear what this index means for resource use data with variation in resource availability. In typical plant-flower visitor webs, PDI values may be close to 1 simply due to highly skewed plant abundance distributions.

%JF: moved here, and changed
	\item When interaction strength is estimated by a count variable (interaction frequency), PDI and many other indices are problematic for species with only one observation (singletons). In this case P1 is 1 (all interactions are on one plant species), all Pi are 0 and hence PDI is 1 - independent of the species' specialization. For a singleton we cannot estimate specificity, only discrimination (i.e. whether it happens to visit a common or rare plant species), as is done by d' (implemented in \code{\link{dfun}}).

}
}

\value{
% JF: changed here to perfect generalist
Returns a vector with PDI values between 0 (perfect generalist) and 1 (perfect specialist).
}

\references{
Dormann, C.F. (2011) How to be a specialist? Quantifying specialisation in pollination networks. \emph{Network Biology} \bold{1}, 1--20

Poisot, T., Lepennetier, G., Martinez, E., Ramsayer, J., and Hochberg, M.E. (2011a) Resource availability affects the structure of a natural bacteria-bacteriophage community. \emph{Biology Letters} \bold{7}, 201--204

Poisot, T., Bever, J.D., Nemri, A., Thrall, P.H., and Hochberg, M.E. (2011b) A conceptual framework for the evolution of ecological specialisation. \emph{Ecology Letters} \bold{14}, 841--851

Poisot, T., E. Canard, N. Mouquet, and M. E. Hochberg. 2012. A comparative study of ecological specialization estimators. \emph{Methods in Ecology and Evolution} \bold{3}, 537--544. doi: 10.1111/j.2041-210X.2011.00174.x.
}

\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }

\seealso{ See also \code{\link{specieslevel}}. }

\examples{
data(Safariland)
PDI(Safariland) # for pollinators
PDI(t(Safariland), log=TRUE) # for plants
}
\keyword{ package}


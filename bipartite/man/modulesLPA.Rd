\encoding{UTF-8}

\name{DIRT_LPA_wb_plus}

\alias{LPA_wb_plus}
\alias{convert2moduleWeb}
\alias{DIRT_LPA_wb_plus}
%\alias{BarbersMatrix}
%\alias{WEIGHTEDMODULARITY}
%\alias{WEIGHTEDMODULARITY2}
%\alias{TRACE}
%\alias{DIVISION}
%\alias{StageOne_LPAwbdash}
%\alias{StageTwo_LPAwbdash}
%\alias{LOCALMAXIMISATION}

\title{Functions "LPA_wb_plus" and "DIRT_LPA_wb_plus"}

\description{
Use \code{\link{computeModules}} to call this function! This function takes a bipartite weighted graph and computes modules by applying Newman's modularity measure in a bipartite weighted version to it. To do so, it uses Stephen J Beckett's DIRTLPAwb+ (or LPAwb+) algorithm, which builds on Liu & Murata's approach. In contrast to the tedious MCMC-swapping QuanBiMo algorithm, this algorithm works by aggregating modules until no further improvement of modularity can be achieved.}

\usage{
DIRT_LPA_wb_plus(MATRIX, mini=4, reps=10)
LPA_wb_plus(MATRIX, initialmoduleguess=NA)
convert2moduleWeb(MATRIX, MODINFO)
}

\arguments{
  \item{MATRIX}{
\code{MATRIX} is the matrix representing the weighted bipartite graph (as an example, see e.g. web \code{\link{small1976}} in this package).
}
 \item{mini}{Minimal number of modules the algorithm should start with; defaults to 4. See explanation of \option{initialmoduleguess} to understand why not starting with one module per species interaction makes sense.}
 \item{reps}{Number of trials to run for each setting of \option{mini}; defaults to 10 but my benefit from higher values.}
 \item{initialmoduleguess}{Optional vector with labels for the modules of the longer web dimension. The \option{initialmoduleguess} argument in \code{LPA_wb_plus} is used as an initial guess of the number of modules the network contains. The code then randomly assigns this many labels across the nodes of one of the types. The classical algorithm (i.e. \code{LPA_wb_plus}) would use an initialmoduleguess of the maximum number of modules a network can contain (i.e. assigning a different label to each node initially). \code{DIRT_LPA_wb_plus} exploits this. By initially assigning fewer than the maximum number of modules in the first instance nodes that may not have been placed together by the classical algorithm are placed together - creating different start points from which to attempt to maximise modularity.Defaults to \code{NA}, i.e. as many modules as there are species in the smaller group.}
 \item{MODINFO}{Object returned by \code{modulesLPA}.}
}

\value{
\code{LPA_wb_plus} computes the modules. \code{DIRT_LPA_wb_plus} is a wrapper calling \code{LPA_wb_plus} repeatedly to avoid getting stuck in some local minimum. Both return a simple list of row and column labels for the modules, as well as the modularity value. Using \code{convert2moduleWeb} turns this into the richer moduleWeb-class object produced by \code{computeModules}. For this object, the plotting function \code{\link{plotModuleWeb}} and summary functions \code{\link{listModuleInformation}} and \code{\link{printoutModuleInformation}} are available.
}

\references{ 
Beckett, S.J. 2016 Improved community detection in weighted bipartite networks. \emph{Royal Society open science} \bold{3}, 140536.

Liu X. & Murata T. 2010. An Efficient Algorithm for Optimizing Bipartite Modularity in Bipartite Networks. \emph{Journal of Advanced Computational Intelligence and Intelligent Informatics (JACIII)} \bold{14} 408--415.

Newman M.E.J. 2004. \emph{Physical Review E} \bold{70} 056131
}

\author{Stephen J Beckett (\url{https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS}), lifted, with consent of the author, by Carsten F. Dormann to bipartite}

\seealso{
\code{\link{computeModules}}; see also class "moduleWeb", \code{\link{listModuleInformation}}, \code{\link{printoutModuleInformation}}
}

\examples{
	\dontrun{
		(res <- DIRT_LPA_wb_plus(small1976))
		mod <- convert2moduleWeb(small1976, res)
		plotModuleWeb(mod)
	}
}
\keyword{Methods and Generic Functions}
\keyword{Clustering}
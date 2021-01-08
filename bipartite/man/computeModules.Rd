\encoding{UTF-8}

\name{computeModules}

\alias{computeModules}
%\alias{cM}
\alias{readModuleData}
%\alias{deleteModuleData}
\alias{metaComputeModules}

\title{computeModules}

\description{
This function takes a bipartite weighted graph and computes modules by applying Newman's modularity measure in a bipartite weighted version to it. 
\code{metaComputeModules} re-runs the algorithm several times, returning the most modular result, to stabilise modularity computation.
}

\usage{
computeModules(web, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
	steps = 1000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
	
metaComputeModules(moduleObject, N=5, method="Beckett", ...)
}

\arguments{
  \item{web}{
\code{web} is the matrix representing the weighted bipartite graph (as an example, see e.g. web \code{\link{small1976}} in this package). This matrix can be binary (i.e. consist only of 0s and 1s), in which case the output will be Newman's (2006) modularity.
}
 \item{method}{Choice between the algorithm(s) provided by Stephen Beckett (2016) or Dormann & Strauss (2016) (\option{method="DormannStrauss"}). Defaults to the much faster and in the majority of cases better algorithm of Beckett. (Note the optional argument \option{forceLPA} to use his slightly inferior but even faster pure LPA algorithm.)}
  \item{deep}{
If \code{deep} is set to FALSE (default), a flat clustering is computed, otherwise submodules are identified recursively within modules. Works only with \option{method="DormannnStrauss"}.
}
  \item{deleteOriginalFiles}{
If \code{deleteOriginalFiles} is set to TRUE (default), the files mentioned above in the description are deleted from the hard drive disk, otherwise not. Applies only to \option{method="DormannnStrauss"}.
}
  \item{steps}{
\code{steps} is the number of steps after which the computation of modules stops if no better division into modules than the current one can be found. Applies only to \option{method="DormannnStrauss"}.
}
 \item{tolerance}{How small should the difference between MCMC-swap results be? At some point computer precision fluctuations make the algorithm fail to converge, which is why we choose a (very low) defaults of 1E-10. Applies only to \option{method="DormannnStrauss"}.}
 \item{experimental}{Logical; using an undescribed and untested version for which no detail is available? (We suggest: not yet.) Applies only to \option{method="DormannnStrauss"}.}
 \item{moduleObject}{Output from running \code{computeModules}.}
 \item{forceLPA}{Logical; should the even faster pure LPA-algorithm of Beckett be used? DIRT-LPA, the default, is less likely to get trapped in a local minimum, but is slightly slower. Defaults to FALSE.  Applies only to \option{method="Beckett"}.}
 \item{N}{Number of replicate runs; defaults to 5. Not really required for \option{method="Beckett"}, which starts in different places anyway.}
 \item{...}{Arguments passed on to \code{computeModules}, which is called internally.}
}

\value{
An object of class "moduleWeb" containing information about the computed modules. For details, please refer to the corresponding documentation file. 
}

\references{ 
Beckett, S.J. 2016 Improved community detection in weighted bipartite networks. \emph{Royal Society open science} \bold{3}, 140536.

Dormann, C. F., and R. Strauß. 2014. Detecting modules in quantitative bipartite networks: the QuanBiMo algorithm. \emph{Methods in Ecology & Evolution} \bold{5} 90--98 (and \bold{arXiv [q-bio.QM]} 1304.3218.)

Liu X. & Murata T. 2010. An Efficient Algorithm for Optimizing Bipartite Modularity in Bipartite Networks. \emph{Journal of Advanced Computational Intelligence and Intelligent Informatics (JACIII)} \bold{14} 408--415.

Newman M.E.J. 2004. \emph{Physical Review E} \bold{70} 056131

Newman, M.E.J. 2006. Modularity and community structure in networks. \emph{Proceedings of the National Academy of Sciences of the United States of America}, \bold{103}, 8577–-8582.

}

\author{Rouven Strauss, with fixes by Carsten Dormann and Tobias Hegemann; modified to accommodate Beckett's algorithm by Carsten Dormann}

\note{For perfectly compartmentalised networks the algorithm may throw an error message. Please add a little bit of noise (e.g. uniform between 0 and 1 or so) or a small constant (1E-5 or so) and it will work again. 

When using the method \option{DormannStrauss}, files are written onto the hard drive during the computation. These files are by default deleted after the computation terminates, unless it breaks. Details of the modularity algorithm can be found in Dormann & Strauß (2013).
}

\seealso{
See also class "moduleWeb", \code{plotModuleWeb}, \code{listModuleInformation}, \code{printoutModuleInformation}, \code{\link{DIRT_LPA_wb_plus}}.
}

\examples{
	\dontrun{
		data(small1976)
		(res <- computeModules(small1976)) 
		plotModuleWeb(res)
		
		# slow:
		res2 <- metaComputeModules(small1976, method="DormannStrauss")
		res2
	}
}
\keyword{Methods and Generic Functions}
\keyword{Clustering}

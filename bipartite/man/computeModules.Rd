\encoding{UTF-8}

\name{computeModules}

\alias{computeModules}
\alias{cM}
\alias{readModuleData}
\alias{deleteModuleData}

\title{Function "computeModules"}

\description{
This function takes a bipartite weighted graph and computes modules by applying M. E. J. Newman's modularity measure in a bipartite weighted version to it. During the computation files are written onto the hard drive disk. These files are by default deleted after the computation terminates. Details of the modularity algorithm can be found in Dormann \& Strauß (2013).
}

\usage{
computeModules(web, deep = FALSE, deleteOriginalFiles = TRUE, 
	steps = 1000000, tolerance = 1e-10, experimental = FALSE)
}

\arguments{
  \item{web}{
\code{web} is the matrix representing the weighted bipartite graph (as an example, see e.g. web "small1976" in this package).
}
  \item{deep}{
If \code{deep} is set to FALSE (default), a flat clustering is computed, otherwise submodules are identified recursively within modules.
}
  \item{deleteOriginalFiles}{
If \code{deleteOriginalFiles} is set to TRUE (default), the files mentioned above in the description are deleted from the hard drive disk, otherwise not.
}
  \item{steps}{
\code{steps} is the number of steps after which the computation of modules stops if no better division into modules than the current one can be found.
}
 \item{tolerance}{How small should the difference between MCMC-swap results be? At some point computer precision fluctuations make the algorithm fail to converge, which is why we choose a (very low) defaults of 1E-10.}
 \item{experimental}{Logical; using an undescribed and untested version for which no detail is available? (We suggest: not yet.)}
}

\value{
An object of class "moduleWeb" containing information about the computed modules. For details, please refer to the corresponding documentation file.
}

\references{ 
Dormann, C. F., and R. Strauß. 2013. Detecting modules in quantitative bipartite networks: the QuaBiMo algorithm. \bold{arXiv [q-bio.QM]} 1304.3218.

Newman M.E.J. 2004. \emph{Physical Review E} \bold{70} 056131
}

\author{Rouven Strauss, with fixes by Carsten Dormann and Tobias Hegemann}

\note{For perfectly compartmentalised networks the algorithm may throw an error message. Please add a little bit of noise (e.g. uniform between 0 and 1 or so) and it will work again. }

\seealso{
See also class "moduleWeb", function "listModuleInformation", function "printoutModuleInformation".
}

\examples{
	\dontrun{
		data(small1976)
		(res <- computeModules(small1976)) # takes several minutes!
		plotModuleWeb(res)
	}
}
\keyword{Methods and Generic Functions}
\keyword{Clustering}

\encoding{UTF-8}

\name{null.distr}

\alias{null.distr}

\title{ Null model based on fitted marginal distribution }

\description{
Given a network, this function fits a distribution to the marginal totals and then draws randomly from this distribution to yield a new network}

\usage{null.distr(N, web, distr="lognormal")}

\arguments{
  \item{N}{Number of null model web to be generated.}
  \item{web}{A (bipartite) network matrix.}
  \item{distr}{The name of the distribution to be fitted to the marginal totals. Currently only the lognormal (\option{"lognormal" or "log-normal"}) and the negative binomial (\option{"negative binomial" or "negbin"}) are supported. Defaults to lognormal.}
}

\details{
  This package provides several functions to generate null models for the observed data (see \code{\link{nullmodel}}). However, this function deviates from any of these in that it does not hold the marginal totals as constant. Rather, it sees them as a random draw of some underlying distribution. This distribution is fitted (in a very ad hoc manner) to the data.  The inspiration for this function comes from Bl¸thgen et al. (2008).
   
 In the next step, the same number of species is drawn from this distribution randomly, and a new matrix is generated as the cross-product of these new vectors. The matrix is then standardised to sum to 1. It now serves as probability of drawing an interaction for any of its cells. 
 
 As many interactions as were observed are drawn (given the above probabilities) and hence a new, null matrix is generated.
 
 In case of the negative binomial fit (and random draw), some 0-values will result from the random draws. As a consequence, the dimension of the matrix can be dramatically lower than the observed. To avoid this, I simply add 1 to each marginal total value. Again, this is very ad hoc and not statistically justified. In fact, values already large should not receive an additional observation (as shown by Dewdney 1998 in a very different context). 
  
 NOTE 1: The fitted distribution is not supposed to represent the true distribution behind the abundances, but merely one way to have new marginal totals. In fact, in many cases the marginal totals aren't lognormal (or negative binomial), but much more skewed than that!  
 
 NOTE 2: Although the dimensions of the new, null web CAN be the same as that of the original, quite often they will be lower. This is because some species have a very low probability of being observed, and only in webs with many observations they will be. Stochasticity may render some species unobserved. \bold{The consequences are potentially large!} As species are lost, relative linkage density goes up automatically, which affects virtually every network index!
}

\value{A list of N new matrices with the same number of interactions, but possibly different dimensions.}

\author{Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>}

\seealso{ \code{\link{nullmodel}}}

\references{
Blüthgen, N., Fründ, J., Vázquez, D. P. and Menzel, F. 2008. What do interaction network metrics tell us about specialisation and biological traits? \emph{Ecology}, \bold{89}, 3387--3399.

Dewdney, A.K. 1998. A general theory of the sampling process with applications to the ``veil line''. \emph{Theoretical Population Biology}, \bold{54}, 294--302.
}

\examples{
\dontrun{
data(Safariland)
null.distr(N=2, Safariland)
null.distr(N=2, Safariland, distr="negbin")

round(networklevel(Safariland, "info"), 3)
sapply(null.distr(N=5, Safariland), function(x) networklevel(x, index="info"))
# highly connected
sapply(null.distr(N=5, Safariland, distr="negbin"), function(x) networklevel(x, 
	index="info")[3])
# similarly highly connected
}
}

\keyword{package}

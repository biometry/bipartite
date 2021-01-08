\encoding{UTF-8}
\name{nullmodel}
\alias{nullmodel}

\title{Generates null models for network analysis}
\description{
A wrapper function for convenient generation of null models for quantitative and binary networks
}
\usage{
nullmodel(web, N=1000, method="r2d", ...)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{N}{number of null models to be generated; defaults to 1000 (more might be better, less probably not).}
  \item{method}{Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web", 3/"vaznull", 4/"shuffle.web", 5/"mgen"; allows for partial match of names; methods 1 to 4 works for quantitative webs, 4 and 5 for binary.}
  \item{...}{arguments to be passed to the function generating the specific null models, see there for options.}
}

\details{
 ADVICE: Look at the same-named function in \pkg{vegan}, as well as the long list of potential null models described in \code{commsim} in that package. It offers a richer and more standardised implementation of null models than this (earlier) function. In particular the method \option{shuffle.web} is potentially confusing, as it calls \pkg{bipartite}'s \code{shuffle.web} for quantitative networks, but \pkg{vegan}'s \option{quasiswab} algorithm for binary. Since \pkg{vegan} also now offers the argument \option{greedyqswap} for quantitative networks, please have a look at \pkg{vegan}'s \code{\link{nullmodel}} function.
 
This is only a wrapper function to facilitate and standardise the generation of null models.
  
These null models assume that integers represent frequencies that are ``individually'' counted, not decimal numbers. Multiplication by 1000 (say) and rounding does NOT necessarily make your value frequencies satisfy this assumption. Null models for ``continuously quantitative'' webs still have to be developed!

%  \option{mgen} is used when options 1, 2 or 3 are applied to a binary network.
A warning is returned when all entries in a quantitative network are 0 or 1 (which suggests a binary network).
}

\value{
  Returns a list of \code{N} null model-generated networks. Species names are (obviously) dropped.
}

\note{
When a quantitative network contains only 1s (as may happen when sampling intensity is low), the quantitative null model will be extremely similar (often identical) to the observed network. This is no error. It is reflecting the fact that this network contains little (no) information beyond the abundances.
}


\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\seealso{ For the functions generating the null model network: \code{\link{shuffle.web}}, \code{\link{swap.web}}, \code{\link{vaznull}}, \code{\link{mgen}}, \code{vegan::simulate} and \code{r2dtable}
}

\examples{
\dontrun{
	data(Safariland)
	nullmodel(Safariland, N=2, method=1)
	nullmodel(Safariland>0, N=2, method=4)
	# analysis example:
	obs <- unlist(networklevel(Safariland, index="weighted nestedness"))
	nulls <- nullmodel(Safariland, N=100, method=1)
	null <- unlist(sapply(nulls, networklevel, index="weighted nestedness")) #takes a while ...
	
	plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))), 
		main="comparison of observed with null model Patefield")
	abline(v=obs, col="red", lwd=2)    
	
	praw <- sum(null>obs) / length(null)
	ifelse(praw > 0.5, 1-praw, praw)    # P-value
	
	# comparison of null model 4 and 5 for binary:
	nulls4 <- nullmodel(Safariland>0, N=100, method=4)
	nulls5 <- nullmodel(Safariland>0, N=100, method=5)
	null4 <- unlist(sapply(nulls4, networklevel, index="weighted nestedness"))
	null5 <- unlist(sapply(nulls5, networklevel, index="weighted nestedness"))
	
	
	plot(density(null4), xlim=range(c(null4, null5)), lwd=2, 
		main="comparison of null models")
	lines(density(null5), col="red", lwd=2)
	legend("topright", c("shuffle", "mgen"), col=c("black", "red"), lwd=c(2,2), 
		bty="n", cex=1.5)
	abline(v=networklevel(Safariland>0, index="weighted nestedness"))
	}
}

\keyword{ package }



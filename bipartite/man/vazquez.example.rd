\encoding{UTF-8}

\name{vazquez.example}
\alias{vazquez.example}
\alias{confint}
\alias{intasymm}
\alias{intereven}
\alias{mlik}
\alias{netstats}
\alias{plotmat}
\alias{quant2bin}
\alias{sortmatr}
\alias{sortmatrext}

\title{ Examples for some analyses }
\description{
Describes how to use bipartite to calculate the statistics presented in Vázquez et al. (2009). Some of these functions are available in bipartite or other packages, and this help page will show how to use them in line with the publication.
}
\details{
  The functions used are:
  \describe{
    \item{confint:}{ Is the same as \code{quantile}}
    \item{intasymm:}{ Can be extracted using \code{\link{specieslevel}}}
    \item{intereven:}{ Is similar to interaction evenness in \code{\link{networklevel}}, but only for a specific option}
 %   \item{mgen:}{Random web based on the number of links; it usually looses ranks (i.e. not all species will still be represented in this random web), and it is based only on the binary web (i.e. all quantitative information is lost); see \code{r2dtable}, \code{\link{vaznull}} and \code{\link{shuffle.web}} for quantitative alternatives.}
    \item{mlik:}{A specific call to \code{dmultinom} and the calculation of the AIC; the number of parameters entering the AIC-calculation is not obvious; this depends on the constraints used by the null model. In the case of \code{r2dtable}, column and row totals are constrained, i.e. ncol+nrow parameters must be given. In the case of \code{\link{swap.web}}, connectance is also constrained, but how many parameters does that imply? One? In \code{\link{shuffle.web}}, we constrain the dimensionality and connectance, i.e. 3 (?) parameters. Vázquez et al. (2009) argue that they constrain only 2 parameters when producing the probability matrix given as pweb in the example below. We tend to disagree: the marginal probabilities of all columns and rows are given, hence k = (ncol(web) + nrow(web)). To our knowledge, there is no mathematical/statistical treatise of this problem.}
    \item{netstats:}{A wrapper calling the other functions, in that sense similar to \code{\link{networklevel}}, but also calling some output from \code{\link{specieslevel}}.}
    \item{plotmat:}{Now part of \code{\link{visweb}}, using the right options.}
    \item{quant2bin:}{A dedicated function to do a simple thing: \code{(web>0)*1}.}
    \item{sortmatr:}{newly defined function: \code{\link{sortweb}}}
    \item{sortmatrext:}{sort matrix by some given sequence; also part of \code{\link{sortweb}}}
  }
  
  In the example below, we use the \pkg{bipartite}/standard R functions whenever possible. 

}

\references{
Vázquez, P.D., Chacoff, N.,P. and  Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology} \bold{90}, 2039--2046.
}

\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> based on code and ideas of Diego Vázquez, Natacha P. Chacoff and Luciano Cagnolo}

\seealso{ See also \code{\link{networklevel}}. }

\examples{
\dontrun{
	data(Safariland)
	
	# confint:
	N100 <- sapply(swap.web(100, Safariland), networklevel, index="nestedness")
	quantile(unlist(N100), c(0.025, 0.975))
	# intasymm: extract values for the asymmetry of interactions and the 
	# dependency matrix for pollinators:
	specieslevel(Safariland)$"higher trophic level"$"interaction push/pull"
	specieslevel(Safariland)$"higher trophic level"$"dependence"
	# for plants:
	specieslevel(Safariland)$"lower trophic level"$"interaction push/pull"
	specieslevel(Safariland)$"lower trophic level"$"dependence"
	
	#intereven
	networklevel(Safariland, index="interaction evenness", intereven="sum")[2]
	# or, as we recommend (see help on networklevel):
	networklevel(Safariland, index="interaction evenness", intereven="prod")[2]
		
	# mlik:
	# calculates the log-likelihood of observing a network, given a probability  
	# matrix of the same size (pweb):
	dmultinom(Safariland>0, prob=pweb, log=TRUE)
	# AIC (the number of parameters is given by how many constraints are put onto the 
	# null model; here, we constrain 9 rows and 27 columns, i.e. sum(dim(binweb))):
	-2*dmultinom(Safariland>0, prob=pweb, log=TRUE) + 2*(sum(dim(binweb)))
	
	# netstats:
	networklevel(Safariland, 
	  index=c("connectance", "interaction evenness", "nestedness", "ISA"))
	mean(specieslevel(Safariland)$"higher trophic level"$"interaction push/pull")
	mean(specieslevel(Safariland)$"lower trophic level"$"interaction push/pull")
	
	#plotmat:
	visweb(t(unname(Safariland)), circles=TRUE, boxes=FALSE)
	
	#sortmatr/sortmatrext:
	sortweb(Safariland, sort.order="inc") #rares species first
	plotweb(sortweb(Safariland, sort.order="dec"), method="normal")
	plotweb(sortweb(web=Safariland, sort.order="seq", 
	  sequence=list(seq.higher=sample(colnames(Safariland)), 
	  seq.lower=sample(rownames(Safariland)))), 
	  method="normal")
	}
}

\keyword{ package}
%	# mgen:
%	binweb <- Safariland>0 #throw away the information on the number of visits
%	# make a matrix with probabilities for each link, based on column and row totals:
%	pweb <- outer(rowSums(binweb)/sum(binweb), colSums(binweb)/sum(binweb), FUN="*")
%	# make a new, emtpy matrix:
%	rbinweb <- matrix(0, nrow=nrow(binweb), ncol=ncol(binweb))
%	# put the links into random places, with probability as given by the observed data:
%	rbinweb[sample(1:prod(dim(binweb)), size=sum(binweb), prob=pweb)] <- 1
%	# this is the new, random realisation given the observed marginal link sums:
%	rbinweb
%	# for this null-web any of the networklevel indices can be calculated


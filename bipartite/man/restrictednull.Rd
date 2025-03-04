\encoding{UTF-8}

\name{restrictednull}
\alias{restrictednull}
\alias{PosteriorProb}

\title{The restricted null model accounting for modules with nested entries}

\description{
Null model for "compound network structure", i.e. modules exhibiting nestedness, roughly maintaining marginal totals and connectance
}
\usage{
restrictednull(web, 
              Prior.Pij = "degreeprob", 
              conditional.level = "modules", 
              N=10, 
              print.null = FALSE, 
              allow.degeneration = FALSE, 
              return.nonrm.species = FALSE, 
              connectance = TRUE, 
              byarea = FALSE, 
              R.partitions = NULL, 
              C.partitions = NULL)
              
PosteriorProb(web, 
              R.partitions = NULL, 
              C.partitions = NULL, 
              Prior.Pij = "degreeprob", 
              conditional.level="modules")
}

\arguments{
	\item{web}{Matrix with observation interactions.}
	\item{Prior.Pij}{Method for computing "a priori" probabilities of interaction between species i and j. Can be defined as:
	  \itemize{
	    \item \option{"degreeprob"}: the default; probability of interaction proportional to overall species degrees;
	    \item \option{"equiprobable"}: probability of interaction identical to all species;
	    \item \option{"degreeprob.byarea"}: probability of interaction proportional to species degrees in each matrix area - see "areas" in \option{conditional.level} for a definition of matrix areas. 
	  }
	}
	\item{conditional.level}{Level to which conditional probability of interaction among species i and j will be conditioned. Can be defined as:
	  \itemize{
	    \item \option{"modules"}: the default; conditional probabilities differing between areas within and outside modules; \option{R.partitions} and \option{C.partitions} must be provided, e.g. based on modularity analyis; see example.
	    \item \option{"matrix"}: conditional probabilities identical in all matrix areas;
      \item \option{"areas"}: a different set of conditional probabilities for each matrix area. A matrix area is a submatrix M[AB] of M formed by all rows of module A and all columns of module B. If A = B, then M[AB] is a module area, otherwise M[AB] is the area between two modules. Therefore, when \option{conditional.level = "areas"}, each area has its own conditional probabilities of interaction.
	  }
	 }
  \item{N}{Integer number of null matrices to be produced.}
  \item{print.null}{Logical: Shall the simulation progress be printed? Defaults to FALSE.}
  \item{allow.degeneration}{Logical. If null matrices are allowed to degenerate, i.e. be of lower dimension that the observed. Defaults to FALSE. If TRUE, interactions are drawn without ensuring that all rows and columns must have at least one interaction.}
  \item{return.nonrm.species}{Logical. Shall the index of non-removed rows and columns should be returned in the output? Defaults to FALSE. Note: if TRUE, each null model will be a list of three, rather than a single matrix!}
  \item{connectance}{Logical. Shall connectance of the null matrices be either exactly (TRUE) or approximately (FALSE) the same as the original matrix. Defaults to TRUE.}
  \item{byarea}{Logical. Shall interactions be drawn independently for each matrix area (i.e. in each submatrix M[ij] of M formed by all rows of module i and all columns of module j)? Defaults to FALSE.}
  \item{R.partitions}{Vector of integers. Partitioning of rows, used only if \option{byarea = TRUE}.}
  \item{C.partitions}{Vector of integers. Partitioning of columns, used only if \option{byarea = TRUE}.}
}


\details{
This is the restricted null model used in Felix et al. (2017), Pinheiro (2019), Pinheiro et al. (2019), Mello et al. (2019) and Queiroz et al. (2020). It was derived from the \code{\link{vaznull}} model. The synthesis presented in this function and null model, reported in a series of studies, was based on the ideas first proposed by Lewinsohn et al. (2006) and followed up by Mello et al. (2009), Flores et al. (2013) and Pinheiro et al (2016).

Our restricted null model was designed for testing for a compound topology, i.e. a modular network structure with internally nested modules. It allows comparing observed and expected values of nestedness between species of the same module (NODFsm), and between species of different modules (NODFdm).
}

\seealso{The function \code{\link{nest.smdm}} computes NODFsm and NODFdm for such networks, while \code{\link{sortmatrix}} and \code{\link{plotmatrix}} facilitate drawing matrices in a way that helps visualizing a compound topology.}

\author{Gabriel Felix, Rafael Pinheiro, and Marco Mello from the Ecological Synthesis Lab (SintECO) in São Paulo. The code was taken, with permission and encouragement, from Gabriel's github repository (\url{https://github.com/gabrielmfelix/Restricted-Null-Model}) and modified by Carsten F. Dormann for conformity with \pkg{bipartite}'s naming style. The helper function \code{PosteriorProb}, which computes the probabilities of an interaction, given the modular structure, is not exported and hence would have to be called as \code{bipartite:::PosteriorProb}. See Gabriel's github page for further details.}

\references{
Bezerra, E. L. S., I. C. Machado, and M. A. R. Mello. 2009. Pollination networks of oil-flowers: a tiny world within the smallest of all worlds. \emph{J. Anim. Ecol.} \bold{78}: 1096--1101. \url{https://pubmed.ncbi.nlm.nih.gov/19515098/}.

Felix, G. M., R. B. P. Pinheiro, R. Poulin, B. R. Krasnov, and M. A. R. Mello. 2017. The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. \emph{bioRxiv} 236687. \doi{10.1101/236687}

Flores, C. O., S. Valverde, and J. S. Weitz. 2013. Multi-scale structure and geographic drivers of cross-infection within marine bacteria and phages. \emph{ISME J.} \bold{7}: 520-–532. \doi{10.1038/ismej.2012.135}.

Lewinsohn, T. M., P. Inácio Prado, P. Jordano, J. Bascompte, and J. M. Olesen. 2006. Structure in plant-animal interaction assemblages. \emph{Oikos} \bold{113}: 174–-184. \doi{10.1111/j.0030-1299.2006.14583.x}

Mello, M. A. R., G. M. Felix, R. B. P. Pinheiro, R. L. Muylaert, C. Geiselman, S. E. Santana, M. Tschapka, N. Lotfi, F. A. Rodrigues, and R. D. Stevens. 2019. Insights into the assembly rules of a continent-wide multilayer network. \emph{Nat. Ecol. Evol.} \bold{3}: 1525--1532. \doi{10.1038/s41559-019-1002-3}

Pinheiro, R. B. P., G. M. F. Félix, A. V Chaves, G. A. Lacorte, F. R. Santos, É. M. Braga, and M. A. R. Mello. 2016. Trade-offs and resource breadth processes as drivers of performance and specificity in a host–parasite system: a new integrative hypothesis. \emph{Int. J. Parasitol.} \bold{46}: 115-–121. \doi{10.1016/j.ijpara.2015.10.002}

Pinheiro, R. B. P., G. M. F. Felix, C. F. Dormann, and M. A. R. Mello. 2019. A new model explaining the origin of different topologies in interaction networks. \emph{Ecology} \bold{100}: e02796. \doi{10.1002/ecy.2796}

}

\examples{
Mod <- computeModules(Safariland)

# Recover the partitions
Part <- module2constraints(Mod)
row.Part <- Part[1:nrow(Safariland)]
col.Part <- Part[(nrow(Safariland)+1):(nrow(Safariland)+ncol(Safariland))]

# Generate randomized networks with the null model of your choice, 
# considering the interaction probabilities calculated before. 
nulls <- restrictednull(web = Safariland, R.partitions = row.Part, C.partitions = col.Part)

# Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, nest.smdm, constraints = Part, weighted = TRUE, decreasing = "abund")
(WNODA.null <- unlist(null[1,])) # WNODArow
(WNODAsm.null <- unlist(null[2,])) # WNODAcol
(WNODAdm.null <- unlist(null[3,])) # WNODAmatrix
# observed values:
nest.smdm(Safariland, weighted = TRUE, decreasing = "abund")
}
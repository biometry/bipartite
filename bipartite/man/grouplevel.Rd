\encoding{UTF-8}
\name{grouplevel}
\alias{grouplevel}
\alias{one.grouplevel}

\title{Analysis of bipartite webs at the level of each of the two levels (groups) of the network}
\description{
	Calculates a variety of indices and values for each group of a bipartite network (\code{one.grouplevel} is the actual function to do the computations and is not intended to be called by the user)
}
\usage{
grouplevel(web, index="ALLBUTDD", level="both", weighted=TRUE, empty.web=TRUE, 
dist="horn", CCfun=mean, logbase="e", normalise=TRUE,  extinctmethod="r", 
nrep=100, fcdist="euclidean", fcweighted=TRUE)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of host.}
  \item{index}{One or more of the following (exact match only!):
   \itemize{
    \item \option{number of species} in the respective trophic level,
    \item \option{mean number of links},
    \item \option{mean number of shared partners},  
    \item \option{cluster coefficient}, 
    \item \option{weighted cluster coefficient},  
    \item \option{togetherness}, 
    \item \option{C score}, 
    \item \option{V ratio}, 
    \item \option{discrepancy}, 
    \item \option{degree distribution}, 
    \item \option{extinction slope}, 
    \item \option{robustness},  
    \item \option{niche overlap},  
    \item \option{generality},
    \item \option{vulnerability},
    \item \option{partner diversity},
    % \item \option{effective partners}, % but still allowed for backward compatibility, returning gen/vul
    \item \option{fc} (or alternatively \option{functional complementarity}), 
    \item \option{ALL} calculates all the above indices (returning a list (since degree distribution fits are a table within this list and cannot be simplified)), 
    \item \option{ALLBUTDD} (default) calculates all indices except degree distribution fits. This latter has the advantage that the output can be simplified to a vector; 
    }
    }
  \item{level}{For which of the two groups (\dQuote{levels}) should these indices be computed? Options are \option{lower}, \option{higher} and \option{both} (default). For index \option{robustness}, the indicated level is the one that will be exterminated in simulations (see \code{\link{second.extinct}}). Although \code{grouplevel} can be employed on its own, it typically will be called through \code{\link{networklevel}}.}
%! JFedit here
\item{weighted}{logical; for those indices which are simply averaged across species to yield the group-level index (e.g. \option{niche overlap}), should this averaging take into account the number of observations for a species? Defaults to \code{TRUE}. The logic behind this default is that we have more faith in an index value when it is based on many, rather than few, observations. Bersier et al. (2002) proposed this weighting and it is also commonly used to compute vulnerability and generality (e.g. in Tylianakis et al. 2006).}
  % In fact, vulnerability/generality are only weighted mean number of effective partners for the lower/higher level, respectively. Thus, with default settings, the same value will be returned twice for these two indices.
  \item{empty.web}{Shall the empty columns and rows be deleted? Defaults to TRUE.}
  \item{dist}{Distance metric to be used to calculate niche overlap (calling \code{vegan::vegdist}); defaults to Horn's index, which is the recommendation of Krebs (1989); for other options see \code{vegdist} in \pkg{vegan}.}
  \item{CCfun}{Method to use when calculating the clustering coefficient. Originally proposed as mean of cluster coefficients for each node. Defaults to \option{median}, because cluster coefficients are strongly skewed. }
%! JFedit here
\item{logbase}{Shall various indices (partner diversity, generality/vulnerability) be calculated to the base of e (default) or 2? Log2 is the proposal for generality and vulnerability by Bersier et al. (2002), while Shannon uses ln. The choice of the base will not affect the results qualitatively, at most by a scaling factor.}
  \item{normalise}{Logical; shall the C-score and togetherness metrics be normalised to a range of 0 to 1? Defaults to \code{TRUE}.}
  \item{extinctmethod}{Specifies how species are removed from matrix: \option{random} or \option{abundance} (partial matching), where abundance removes species in the order of increasing abundance (i.e. rarest first); see Memmott et al. (2004).}
  \item{nrep}{Number of replicates for the extinction sequence analysis.}
  \item{fcweighted}{Logical; when computing "functional complementarity" sensu function \code{\link{fc}}, should the weights of the matrix be used. Defaults to TRUE, but original paper (Devoto et al. 2012) is based on FALSE.}
  \item{fcdist}{Distance measure to be used to compute functional complementarity through \code{\link{fc}}; any measure accepted by \code{\link{dist}} is acceptable.}
}

\details{
This function implements a variety of the many (and still procreating) indices describing network topography at the group level.

%! JFedit here
Note that Bersier et al. (2002) have three levels of values for some of their indices: qualitative (i.e. based on binary networks), quantitative (based on networks with information on the number of interactions observed for each link), and weighted-quantitative (where each species is given a weight according to the number of interactions it has). At present, we implement a mixture of qualitative, quantitative and weighted-quantitative indices and offer the option \option{weighted} to compute the weighted-quantitative version of some of them (\option{mean number of links}, \option{mean number of shared partners}, \option{cluster coefficient}, \option{partner diversity}, \option{generality / vulnerability}). For all others, the mechanics behind the index do not allow a weighted mean to be computed (e.g. the distance-matrix between all species combinations used to compute \option{niche.overlap}).

All indices in this function work with real as well as integer values. 

Extinction slope works on a repeated random sequence of species extinctions (within one trophic level), and calculates the number of secondary extinctions (in the other level). These values are then averaged (over the \option{nrep} runs) and plotted against the number of species exterminated. The proportion still recent (on the y-axis) regressed against the proportion exterminated (on the x-axis) is hence standardised to values between 0 and 1 each. Through this plot, a hyperbolic regression is fitted, and the slope of this regression line is returned as an index of extinction sensitivity. The larger the slope, the later the extinction takes its toll on the other trophic level, and hence the higher the redundancy in the trophic level under consideration. Using \option{plot.it=F} also returns the graphs (set history to recording in the plotting window). Changing the \option{extinctionmethod} to \dQuote{abundance} will always result in the same sequence (by increasing abundance) and hence does not require replication.

Most indices are straightforward, one-line formulae; some, such as betweenness, also require a re-arranging of the matrix; and one, secondary extinction slope, internally requires iterative runs, making the function relatively slow. If you are not interested in the secondary extinction slopes, simply set \option{nrep=1} to make it much faster.


}

\value{
The suffixes LL and HL refer to lower and higher level, respectively. If values for both levels are requested, those for the higher level are given first, followed immediately by those for the lower level.

Depending on the selected indices, some or all of the below (returned as vector if \dQuote{degree distribution} was not requested, otherwise as list):
  \item{mean number of species}{\emph{sic}, possibly weighted (if \option{weighted=TRUE}; the weighted mean is not something typically reported, but it seems a very plausible way to embrace the uncertainty introduced by species with very few interactions).}
  \item{mean number of links}{\emph{sic} (sum of links for each species, averaged over all species in that level), possibly weighted (if \option{weighted=TRUE}).}
  \item{mean number of shared partners}{Based on the distance matrix between species, counting the number of species in the other level that both interact with; based on Roberts & Stone (1990) and Stone & Roberts (1992), i.e. for pollinators will yield mean number of plants shared by any two pollinators (Cannot be weighted.)}
  \item{cluster coefficient}{The cluster coefficient for a level is the (weighted) average cluster coefficients of its members. The cluster coefficient for each species, in turn, is simply the number of realised links divided by the number of possible links. Introduced by Watts & Strogatz (1998) and described in Wikipedia under \url{http://en.wikipedia.org/w/index.php?title=Clustering_coefficient}. If you want to use Tore Opsahl's adaptation to two-modes, please see the next index, based on his function \code{clustering_tm} in \pkg{tnet}. To my knowledge, so far every study has used the \dQuote{wrong} one, i.e. the one presented here as \option{cluster coefficient}.}
  \item{weighted cluster coefficient}{When asking for \dQuote{weighted cluster coefficient}, this version will automatically use interactions as weights unless the data are binary. The computation is based on \code{clustering_tm} in \pkg{tnet}. See there (and more on Tore Opsahl's webpages) for help.}
  \item{niche overlap}{Mean similarity in interaction pattern between species of that level, calculated by default as Horn's index (use \option{dist} to change this.). Values near 0 indicate no common use of niches, 1 indicates perfect niche overlap.}
  %(In how far it makes sense for hosts of predators to call their commonness in enemies \dQuote{niche overlap} is a different issue. There are people calling predators \dQuote{negative resources} (couldn't be asked to look up the reference). I would probably rather call it similarity in predators in this case.)
  \item{togetherness}{Mean number of co-occupancies across all species-combinations; the whole matrix is scanned for submatrices of the form (0,0,1,1), representing perfect matches of co-presences and co-absences. These are counted for each pairwise species combination, and averaged (without weighting). Since the number of species differs between the levels, the same number of co-occupancies will lead to different togetherness-values for the two levels. Based on Stone & Roberts (1992).}
  \item{C score}{(Normalised) mean number of checkerboard combinations across all species of the level. Values close to 1 indicate that there is evidence for disaggregation, e.g. through competition. Value close to 0 indicate aggregation of species (i.e. no repelling forces between species). Since the number of species differs between the levels, the same number of checkerboard patterns will lead to different C-scores for the two levels. See Stone and Roberts (1990) for details.}
  \item{V ratio}{Variance-ratio of species numbers to interaction numbers within species of a level. Values larger than 1 indicate positive aggregation, values between 0 and 1 indicate disaggregation of species. See Schluter (1984) for details.}
  \item{discrepancy}{Discrepancy as proposed by Brualdi & Sanderson (1999); see \code{\link{discrepancy}} for details.}
  \item{degree distribution}{Coefficients and fits for three different functions to a level's degree distributions: exponential, power law and truncated power law. See \code{\link{degreedistr}} for details and references.}
  \item{extinction slope}{Slope of the secondary extinction sequence in one level, following extermination of species in the other level; \code{extinction slope.HL} refers to the robustness of the higher level to extinctions in the lower level (and vice versa); see \code{\link{slope.bipartite}} and \code{\link{second.extinct}} for details.}
  \item{robustness}{Calculates the area below the \dQuote{secondary extinction} curve; \code{robustness.HL} refers to the robustness of the higher level to extinctions in the lower level (and vice versa); see \code{\link{robustness}} for details. Corresponds to \dQuote{extinction slope}.}
  \item{functional complementarity}{\dQuote{Functional complementarity} for a given level. This measure of niche complementarity (as described by Devoto et al. 2012), is computed as the total branch length of a \dQuote{functional dendrogram} based on qualitative differences of interactions of one level with the other. Thus, the \dQuote{functional} aspect of functional complementarity refers to the function of sharing interactions. Should be highly correlated with niche overlap, only binary.}
\item{partner diversity}{(Weighted) mean Shannon diversity of the number of interactions for the species of that level. Choose \option{logbase=2} to change to a log2-based version.}
%#! JFedit: effpartners and gen/vul combined
%\item{effective partners}{(Weighted) mean of exp(Shannon diversity). Simply the Jost (2006)-recommended version of diversity, which happens to be identical to generality (higher level) or vulnerability (lower level), as defined in Bersier et al. (2002).}
 \item{generality/vulnerability}{(Weighted) mean effective number of LL species per HL species (generality; HL species per LL species for vulnerability), weighted by their marginal totals (row sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical to exp(\dQuote{partner diversity}, i.e., simply the Jost (2006)-recommended version of diversity.}
}

\note{
If your web returns and NA for some of the indices, this can be because the index cannot be computed. For example, if the web is full (i.e. no 0-cells), extinction slopes cannot be fitted (singularity of gradient). Check if you can expect the index to be computable! If it is, and \code{grouplevel} doesn't do it, let me know.

Some indices require rather long computation times on large webs. If you want to increase the speed by omitting some indices, here a rough guide: Ask only for the indices you are interested in! Otherwise, here is the sequence of most time-consuming indices: 
    \enumerate{
	 \item For somewhat larger networks (i.e. more than 2 dozen species per level), \option{weighted cluster coefficient} is \emph{very} time consuming (an exhaustive search for 4-loops in the one-mode projection of the network). Omitting it can dramatically boost speed.
	 \item Typically, the slowest function is related to extinction slopes and robustness. Excluding \emph{both} makes the function faster. 
	 \item Degree distributions are somewhat time consuming. 
   }
}

\references{ 

Bascompte, J., Jordano, P. and Olesen, J. M. 2006. Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. 2002. Quantitative descriptors of food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N. 2010. Why network analysis is often disconnected from community ecology: A critique and an ecologist's guide. \emph{Basic and Applied Ecology} \bold{11}, 185--195

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen N. 2007 Specialization, constraints and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Devoto M., Bailey S., Craze P., and Memmott J. (2012) Understanding and planning ecological restoration of plant-pollinator networks. \emph{Ecology Letters} \bold{15}, 319--328. http://dx.doi.org/10.1111/j.1461-0248.2012.01740.x

Dormann, C.F., Fründ, J., Blüthgen, N., and Gruber, B. (2009) Indices, graphs and null models: analysing bipartite ecological networks. \emph{The Open Ecology Journal} \bold{2}, 7--24.

Dunne, J. A., R. J. Williams, and N. D. Martinez. 2002 Food-web structure and network theory: the role of connectance and size. \emph{Proceedings of the National Academy of Science USA} \bold{99}, 12917--12922

Gotelli, N. J., and G. R. Graves. 1996 Null Models in Ecology. Smithsonian Institution Press, Washington D.C.

Jost, L. 2006. Entropy and diversity. \emph{Oikos} \bold{113}, 363-–375. doi: 10.1111/j.2006.0030-1299.14714.x

Krebs, C. J. 1989 \emph{Ecological Methodology}. Harper Collins, New York.

Memmott, J., Waser, N. M. and Price M. V. 2004 Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Roberts, A. and Stone, L. 1990 Island-sharing by archipelago species. \emph{Oecologia} \bold{83}, 560--567

Schluter, D. (1984) A variance test for detecting species associations, with some example applications. \emph{Ecology} \bold{65}, 998-1005.

Stone, L. and Roberts, A. (1990) The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79.

Stone, L. and Roberts, A. (1992) Competitive exclusion, or species aggregation? An aid in deciding. \emph{Oecologia} \bold{91}, 419--424

Tylianakis, J. M., Tscharntke, T. and Lewis, O.T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Watts, D. J. and Strogatz, S. (1998) Collective dynamics of \sQuote{small-world} networks. \emph{Nature} \bold{393}, 440--442
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\seealso{ This function can (and typically will) be called, with all its arguments, by \code{\link{networklevel}}. Several indices have their own function as implementation: \code{\link{second.extinct}}, \code{\link{degreedistr}}, \code{\link{C.score}} and \code{\link{V.ratio}}
}

\examples{
\dontrun{
data(Safariland)
grouplevel(Safariland)
grouplevel(Safariland, level="lower", weighted=FALSE) #excludes degree distribution fits
}
}

\keyword{ package }
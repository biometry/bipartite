\encoding{UTF-8}
\name{specieslevel}

\alias{specieslevel}

\title{Calculate various indices for network properties at the species level}

\description{
Apart from the properties of the entire web, also its participants can be described
specifically. Various simple numbers and indices are calculated and returned.
}

\usage{
specieslevel(web, index="ALLBUTD", level="both", logbase=exp(1), low.abun=NULL, 
	high.abun=NULL, PDI.normalise=TRUE, PSI.beta=c(1,0), nested.method="NODF", 
	nested.normalised=TRUE, nested.weighted=TRUE, empty.web=TRUE)
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher
  trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of host.}
  \item{index}{Vector of indices to be calculated for each trophic level of the web;
            options are: 
    \itemize{
       		\item \option{degree}, 
		\item \option{ND} for normalised degrees, 
		%\item \option{dependence}, 
		\item \option{species strength} as sum of dependencies for each species, 
		\item \option{nestedrank} as rank in a nested matrix,
		\item \option{interaction push pull} for interaction push/pull (our version of dependence asymmetry: see details),  
	        	\item \option{PDI} for Paired Differences Index, 
	        	\item \option{resource range} for Poisot et al. (2012)'s index of unused resources,
	        	\item \option{species specificity} (or coefficient of variation of interactions),
        		\item \option{PSI} for pollination service index (or pollinator support index, depending on the trophic level), 
	        	\item \option{NS} for node specialisation index, 
		\item \option{betweenness} for betweenness, 
		\item \option{closeness} (both automatically also return their weighted counterparts proposed by Tore Opsahl in package \pkg{tnet}),
        		\item \option{Fisher} for Fisher's alpha index, 
		\item \option{diversity} for Shannon diversity of interactions of that species, 
		\item \option{effective partners} for the effective number of interacting partners,  
		\item\option{proportional generality} a quantitative version of normalised degree,
		\item\option{proportional similarity} specialisation measured as similarity between use and availability,
		\item \option{d} for Blüthgen's d', 
		\item \option{ALL} for all the aforementioned. 
	        	\item \option{ALLBUTD} (default) excludes only the dependence matrix and leads to the output being simplified to a matrix per trophic level.}
            }
  \item{level}{For which level(s) should the indices be computed? Options are \option{both} (default), \option{lower} and \option{higher}. Output for one trophic level will be returned as a data.frame (unless dependencies are also reported), while for two levels a list of two matrices (higher and lower) will be returned.}
  \item{logbase}{numeric; gives the base to which the various diversity indices (partner diversity, effective partners) shallb be calculated, typcially exp(1) (default) or 2. Log2 is the proposal for generality and vulnerability by Bersier et al. (2002), while Shannon typically uses the natural logarithm, ln. The choice of the base will not affect the results qualitatively, only by a scaling factor.}
  \item{low.abun}{Optional vector of independent abundances of the \emph{lower} trophic level to be passed on to \code{\link{dfun}} and used as `resource availability' in indices proportional similarity and proportional generality.}
  \item{high.abun}{Optional vector of independent abundances of the \emph{higher} trophic level to be passed on to \code{\link{dfun}} and used as `resource availability' in indices proportional similarity and proportional generality.}
  \item{PDI.normalise}{Logical; should the normalised PDI be computed? Defaults to TRUE. See function \code{\link{PDI}} for details and reasoning.}
  \item{PSI.beta}{A length 2 vector of parameter determining the effect of floral constancy and pollen mixing on the proportion of same-species pollen delivered to a plant. Defaults to \code{c(1, 0)}. See section details for details.}
  \item{nested.method}{One of \option{NODF}, \option{binmatnest}, \option{wine} or \option{sort}. See \code{\link{nestedrank}} for details.}
  \item{nested.normalised}{Logical; defaulting to TRUE. Divides the nestedrank-1 by the number of species -1, thereby ranging it between 0 (most generalist) and 1 (most specialised). Species sequences from different networks are now comparable.}
  \item{nested.weighted}{Logical; should the ``nestedrank'' be computed based on weighted network (default) or binary version (FALSE)?}
  \item{empty.web}{Shall the empty columns and rows be deleted? Defaults to TRUE. option FALSE not yet fully tested}
}

\details{
This function implements a variety of the many (and still procreating) indices describing species properties. Some are embarrassingly simple (such as number of interacting species for each species). Others are variations on Shannon's diversity index applied to within species.

Indices based on graph theory (such as NDI, closeness, betweenness) require the data to form a connected graph. When the network is compartmented (as would be seen when plotting it using \code{\link{plotweb}}), these indices will be computed for the each compartment. However, single-link compartments (only one partner in each trophic level) will not form a proper graph and hence the indices will have a value of NA.

Most indices are straightforward, one-line formulae; some, such as d', also require a re-arranging of the matrix. We (Dormann, Blüthgen, Gruber) came up with a new one, called \dQuote{Pollination Service Index} or psi, for which a few more details seem appropriate.

%PSI calculates the average contribution per visit for each pollinator species (which in itself depends 
%on the specialisation and abundance of the bees, as well as the abundance of the plant species). It is hence a measure
%of the \sQuote{usefulness} of each pollinator for the observed plant community. The value is ranged
%between 0 (pollinator largely useless for this plant community) and 1 (perfect specialisation by abundant pollinators).
%The usefulness of this value lies in the fact that a secondary extinctions will increase as the psi-value of a
%extinct pollinator increases: psi=1 implies a full dependence of at least one species on this pollinator, while a
%psi of 0 indicates a pollinators redundancy.

\bold{Pollination Service Index (PSI)}

This index estimates the importance of a pollinator for all plant species. PSI is comprised of three calculation steps: firstly, we calculate, for each pollinator species, the proportion to which it visits each plant species (or, phrased anthropomorphically, the number to the question: which proportion of my visits are to dandelion?). Secondly, we calculate the proportion to which a plant is visited by each bee species (Which proportion of my pollinators are red mason bees?). Multiplying, these two proportions gives the portion of own pollen for each plant species (because this depends both on a pollinators specialisation (step 1) and the plant's specific receptiveness (step 2). Finally, we sum the proportions own pollen delivered across all plant species. This value is the PSI-value. At its maximum (which depends on the value of beta) it shows that all pollen is delivered to one plant species that completely depends on the monolectic pollinator. At its minimum, 0, it indicates that a pollinator is irrelevant to all plant species. Note that PSI can assume values from 0 to 1 (for beta=1) for species of any frequency: a bee been found only once on a plant species visited by no-one else receives a PSI of 1, even if in total 14 million visits were recorded.

(This is all very complicated. So here is another attempt (by Jochen) to explain the PSI:
For PSI, importance of a pairwise interaction (for the plant) is calculated as: 'dependence'_i_on_j * per.visit.efficiency_i_visitedby_j, where per.visit.efficiency_i_visitedby_j = (average proportion visits to i by j in all visits by j)^beta.

It assumes that the order of plant species visited is random (no mixing, no constancy). To account for that not being true, beta could be adjusted. However, this really waits for good empirical tests.)

We envisage a penalty for the fact that a pollinator has to make two (more or less successive) visits to the same plant species: the first to take the pollen up, the second to pollinate the next. Thus, using beta=2 as an exponent in step 1 would simulate that a pollinator deposits all pollen at every visit. In a sense, beta=2 represents a complete turnover of pollen on the pollinator from one visit to the next; only the pollen of the last-visited species is transferred. That is certainly a very strong penalisation. At present we set the exponent to beta=1, because the step of controlling for \dQuote{pollen purity} is already a major improvement. It assumes, implicitly, that pollen is perfectly mixed on the pollinator and hence pollen deposited directly proportional to frequency of visits to the different plants Also, we have no idea to which extent pollen gets mixed and/or lost during foraging flights, and the true exponent remains elusive. For a value of beta=0, PSI simplifies (and is equal) to species strength.

For the perspective of the plant's effect on pollinators (then PSI = pollinator support index), this index makes less sense. Here we would rather use beta=0, because %any plant's pollen is as good as any others for the well-being of the pollinator. This assumption is certainly simplistic, but 
pollen value is not related to number of visits, so we cannot compute it from the network. %(Using beta=1 and the inverse of the matrix, this would give the proportion of self-pollen from a given community. This could serve as an index to compare pollination service from different webs. Interesting!)
Similarly, for other networks, such as host-parasitoids, beta=0 seems plausible, since for the host it does not matter, whether a parasitoid has visited another species before or not. In this case (beta=0), PSI is simply equal to species strength. Not just pollen turnover/carryover on the pollinator is important and influences beta, but all these considerations depend on the assumption how the proportion of conspecific pollen affects pollination (assuming many visits per flower visitation sequence). (a) If only presence of any conspecific pollen on bee is sufficient for pollination, carryover (how long pollen from one visits remains on bee) matters, beta is anywhere between 0 (infinite carryover) and 2 (one-step carryover).
(b) If the proportion of conspecific pollen on bee determines pollination success (linear relationship), carryover does not matter, the proportion can be assumed to be in an equilibrium, and beta=1.

Our choice of defaults (c(1,0)) will yield species strength for plants, and PSI for pollinators, assuming, for the latter, that pollen mixes perfectly.
}
\value{
For both the \dQuote{higher trophic level} and the \dQuote{lower trophic level} a list with the following components:

%  \item{number of species}{Sic.}
  \item{species degree}{Sum of links per species.}
  \item{normalised degree}{As degree, but scaled by the number of possible partners; see \code{\link{ND}}.}
 % \item{dependence}{Web entry divided by column totals (for higher trophic level) or row totals (for lower trophic level). Expresses the proportion of interactions with each species, and sums to 1 for each.}
  \item{species strength}{Sum of dependencies of each species (used, e.g., in Bascompte et al. 2006). It aims at quantifying a species' relevance across all its partners. The alternative version of Barrat et al. (2004; also used by Poisot et al. 2012) as the sum of interactions of a species seems too trivial a measure, reflecting abundance rather than anything else. Do not take this to be the much-discussed \dQuote{interaction strength} in food web papers, which focusses on pairwise interactions (reviewed in Berlow et al. 2004)!}
  \item{interaction push pull}{Direction of interaction asymmetry based on dependencies: positive values indicate that a species affects the species of the other level it interacts with stronger then they affect it (\dQuote{pusher}); negative values indicate that a species is, on average, on the receiving end of the stick (\dQuote{being pulled}); formula based on Vázquez et al. (2007), but push/pull is our own nomenclature. Values are highly correlated with species strengths (see below), but standardised to fall between -1 (being pulled) and 1 (pushing). Compared to \dQuote{strength}, this index quantifies the net balance, rather than the average effect.}
  \item{nestedrank}{Quantifies generalism by the rank of a species in a network matrix re-arranged for maximal nestedness (Alarcon et al. 2008). A low rank (e.g. 1, 2) indicates high generality, while high ranks (up to the number of species in that level) indicate specialism or rarity.}
  \item{PDI}{Paired Differences Index as proposed by Poisot et al. (2011a,b), by default using a normalised version (\option{PDI.normalise=TRUE}); ranges between 0 (generalist) and 1 (specialist); see \code{\link{PDI}} for details and comments.}
  \item{resource range}{Poisot et al.'s (2012) ``resource range'' is a somewhat strange name for something that has a value of 0 when all resources are used, but a value of 1 when only one resource is used. It is, in fact, closer to an ``unused resource range''.  The aforementioned Paired Difference Index is a generalisation of resource range, which is equal to resource range when the web is binary.}
  \item{species specificity}{Coefficient of variation of interactions, normalised to values between 0 and 1, following the idea of Julliard et al. (2006), as proposed by Poisot et al. (2012). Values of 0 indicate low, those of 1 a high variability (and hence suggesting low and high specificity). Since not corrected for number of observations, this index will yield high specificity for singletons, even though there is no information to support this conclusion.} 
  \item{PSI}{Pollination Service Index for the higher trophic level, and the equivalent Pollinator Support Index for the lower trophic level. See Details above for more explanations.} %Expresses the value of, say, a pollinator for all plant species as the sum (across all plant species) of the proportion of pollen deposited per visit. The more specialised the pollinator and the fewer pollinators the plants it pollinates, the higher is its value for the plant community. 
%It is a really cool measure, but it's not easy to get one's head around it. Alternatively, one can see it as the inverse of redundancy of the species for the web: the higher the value, the less redundant.}
  \item{node specialisation index}{Another measure of specialisation, based on the path length between any two higher-trophic level species. Species sharing hosts/prey have an FS-value of 1. See specific function \code{\link{nodespec}} for details, problems and reference.}
  \item{betweenness}{A value describing the centrality of a species in the network by its position on the shortest paths between other nodes; see \code{\link{BC}} and \code{betweenness} in \pkg{sna} and the dedicated section in the vignette (sec. 5.4.1).}
  \item{weighted betweenness}{Computes betweenness (proportion of shortest paths through this species), but based on weighted representation of the network. It calls \code{betweenness_w} from \pkg{tnet} and often differs considerably from its binary counterpart!}
  \item{closeness}{A value describing the centrality of a species in the network by its path lengths to other nodes; see \code{\link{CC}} and \code{closeness} in \pkg{sna}.}
  \item{weighted closeness}{Computes closeness (in one of its varieties), but based on weighted representation of the network. It calls 
\code{closeness_w} from \pkg{tnet} and is usually very similar to its binary counterpart. Note that NAs indicate that these species belong to a different compartment and hence no closeness distance could be calculated.}
  \item{Fisher alpha}{Fisher's alpha diversity for each species (see \code{fisher.alpha} in \pkg{vegan} for details).}
  \item{partner diversity}{Shannon diversity (when using \option{logbase="e"}) or per-species generality/vulnerability (when using 
\option{logbase=2}) of the interactions of each species. See also \code{\link{networklevel}} for the aggregated version of this index (i.e. averaged across all species in a trophic level).}
  \item{effective partners}{\option{logbase} to the power of \dQuote{partner.diversity}: Bersier et al. (2002) interpret this as the effective 
number of partners, if each partner was equally common. Note that \dQuote{partner} is a bit euphemistic when it comes to predator-prey or host-parasitoid networks.}
  \item{proportional generality}{`Effective partners' divided by effective number of resources (`logbase' to the power of `resource diversity'; which is calculated from high.abun/low.abun if provided, and else from marginal totals); this is the quantitative version of proportional resource use or normalised degree (i.e., the number of partner species in relation to the potential number of partner species); note that this index can be larger than 1, e.g. when a species selects for a balanced diet.}
  \item{proportional similarity}{Specialization measured as dissimilarity between resource use and availability (estimated from high.abun/low.abun if provided, else from marginal totals); proposed by Feinsinger et al. (1981).}
  \item{d}{Specialisation of each species based on its discrimination from random selection of partners. More specifically, it returns d', which is calculated based on the  raw \code{d}, \code{dmin} and \code{dmax} for each species (see \code{\link{dfun}}). See Blüthgen et al. (2006) for details.}
}

\references{
% JF edit: I just changed the reference for resource range to Poisot et al. 2012 and deleted the Schoener cite (who doesn't mention anything like this index); Schoener 1989 would rather be a citation for generality(generalization) and vulnerability in their classic (binary) versions (before Bersier, vulnerability was defined this way among foodwebbers, compare the help for networklevel)
% JF edit: also changed (under "value") the description of PDI and resource range

Alarcon, R., Waser, N.M. and Ollerton, J. 2008. Year-to-year variation in the topology of a plant-pollinator interaction network. \emph{Oikos} \bold{117}, 1796--1807

 Barrat, A., M. Barthélemy, R. Pastor-Satorras, and A. Vespignani. 2004. The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences of the USA} \bold{101}, 3747–-3752. doi: 10.1073/pnas.0400087101.
 
Bascompte, J., Jordano, P. and Olesen, J. M. (2006) Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Berlow, E. L., A. M. Neutel, J. E. Cohen, P. C. de Ruiter, B. Ebenman, M. Emmerson, J. W. Fox, V. A. A. Jansen, J. I. Jones, G. D. Kokkoris, D. O. Logofet, A. J. McKane, J. M. Montoya & O. Petchey (2004) Interaction strengths in food webs: issues and opportunities. \emph{Journal of Animal Ecology} \bold{73}, 585-–598

Blüthgen, N., Menzel, F. and Blüthgen, N. (2006) Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 9

Dormann, C.F. (2011) How to be a specialist? Quantifying specialisation in pollination networks. \emph{Network Biology} \bold{1}, 1--20

Feinsinger, P., Spears, E.E. and Poole,R. W. (1981) A simple measure of niche breadth. \emph{Ecology} \bold{62}, 27--32.

Julliard, R., Clavel, J., Devictor, V., Jiguet, F. and Couvet, D. (2006) Spatial segregation of specialists and generalists in bird communities. \emph{Ecology Letters} \bold{9}, 1237-–1244

Martín Gonzáles, A.M., Dalsgaard, B. and Olesen, J.M. (2010) Centrality measures and the importance of generalist species in pollination networks. \emph{Ecological Complexity}, \bold{7}, 36--43

Opsahl, T. & Panzarasa, P. (2009). Clustering in weighted networks. Social Networks, 31, 155--163

Poisot, T., Lepennetier, G., Martinez, E., Ramsayer, J., and Hochberg, M.E. (2011a) Resource availability affects the structure of a natural bacteria-bacteriophage community. \emph{Biology Letters} \bold{7}, 201--204

Poisot, T., Bever, J.D., Nemri, A., Thrall, P.H., and Hochberg, M.E. (2011b) A conceptual framework for the evolution of ecological specialisation. \emph{Ecology Letters} \bold{14}, 841--851

Poisot, T., E. Canard, N. Mouquet, and M. E. Hochberg (2012) A comparative study of ecological specialization estimators. \emph{Methods in Ecology and Evolution} \bold{3}, 537-–544. doi: 10.1111/j.2041-210X.2011.00174.x.

Vázquez, D. P., Melian, C. J., Williams, N. M., Blüthgen N., Krasnov B. R. and Poulin, R. (2007) Species abundance and asymmetric interaction strength in ecological networks. \emph{Oikos} \bold{116}, 1120--1127
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} & Jochen Fründ}

\note{
A comparison of specialisation indices is provided in Dormann (2011); the PDI is missing (since it was published later).

Dependencies are still an open field of debate. Dependencies are calculated as the value in a matrix divided by the \code{rowSums} (for the lower trophic level) or the \code{colSums} (for the higher trophic level). As such, any pollinator observed only once will receive a dependency-value of 1, indicating perfect dependence on this plant species. That may or may not be true. In any case it is based on a sample size of 1, that is why the dependency asymmetry (which is based on the dependencies for both trophic levels) has come under criticism and may be rather sensitive to singletons.

We here provide the code to calculate the strength of a species (i.e. sum of its dependencies), based on the current proposal by Bascompte et al. (2006). It may be a good idea to remove all singletons from the web before calculating this index, to investigate whether it is indeed driven by those scarce observations.

The maximum value of the uncorrected d (i.e. the maximal potential specialization defining d'=1) is not a trivial issue. We treat it here in the same way as given in the BMC Ecology paper, but please have a look at the (raw) code for further comments.
}

\seealso{ \code{\link{networklevel}} for some further comments; \code{\link{dfun}}, \code{\link{nodespec}}, which are called by this function }

\examples{
data(Safariland)
\dontrun{
specieslevel(Safariland)}
specieslevel(Safariland, index="ALLBUTD")[[2]]
}


\keyword{ package}

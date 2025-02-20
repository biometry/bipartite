 \encoding{UTF-8}
\name{networklevel}
\alias{networklevel}


\title{Analysis of bipartite webs at the level of the entire network}
\description{
	Calculates a variety of indices and values for a bipartite network
}
\usage{
networklevel(web, index="ALLBUTDD", level="both", weighted=TRUE, 
   ISAmethod="Bluethgen",  SAmethod = "Bluethgen", extinctmethod = "r", 
   nrep = 100, CCfun=median, dist="horn", normalise=TRUE, empty.web=TRUE, 
   logbase="e", intereven="prod", H2_integer=TRUE, fcweighted=TRUE, 
   fcdist="euclidean", effective=FALSE, legacy=FALSE, mass.action.norm=FALSE)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). 
Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{index}{One or more of the following (exact match only!). First the group of \dQuote{pure} network indices, then those computed for each level.
   \itemize{
      \item \option{connectance},
      \item \option{web asymmetry},
      \item \option{links per species}, 
      \item \option{number of compartments}, 
      \item \option{compartment diversity}, 
      \item \option{cluster coefficient}, which will compute both the network-wide cluster coefficient as well as those for each level,
      \item \option{nestedness} (where ties may induce variation when changing the sequence of species in \option{web}), 
      \item \option{NODF},
      \item \option{weighted nestedness} (where ties may induce variation when changing the sequence of species in \option{web}), 
      \item \option{weighted NODF}, 
      \item \option{spectral radius},
      \item \option{ISA} (or alternatively \option{interaction strength asymmetry} or 
      	\cr \option{dependence asymmetry}), 
      \item \option{SA} (or alternatively \option{specialisation asymmetry}),
      \item \option{linkage density}, 
      \item \option{weighted connectance},
      \item \option{Fisher alpha}, 
 %   \item \option{mean interaction diversity}, 
      \item \option{interaction evenness}, 
      \item \option{Alatalo interaction evenness}, 
      \item \option{Shannon diversity}, 
      \item \option{H2};  
   } and/or those invoked through \code{\link{grouplevel}}:
   \itemize{
     \item \option{number of species} in the respective trophic level,
     \item \option{mean number of links},
     \item \option{mean number of shared partners},  
     \item \option{weighted cluster coefficient},  
     \item \option{degree distribution}, 
     \item \option{togetherness}, 
     \item \option{C score}, 
     \item \option{V ratio}, 
     \item \option{discrepancy} (where ties may induce variation when changing the sequence of species in \option{web}), 
     \item \option{extinction slope} (where ties may induce variation when changing the sequence of species in \option{web}), 
     \item \option{robustness} (where ties may induce variation when changing the sequence of species in \option{web}),  
      \item \option{niche overlap},  
      \item \option{generality},
      \item \option{vulnerability},
      %\item \option{effective partners}, 
      \item \option{fc} (or alternatively \option{functional complementarity}).
    } Furthermore, there are some groups of indices that can be called:
    \itemize{
      \item \option{ALL} calculates all indices (returning a list (since degree distribution fits are a table within this list and cannot be simplified)), 
     \item \option{ALLBUTDD} (default) calculates all indices except degree distribution fits. This latter has the advantage that the output can be simplified to a vector; 
     \item \option{info} returns more general information on the network;
     \item \option{binary} returns a best-of selection of indices based on a binary network;
     \item \option{quantitative} returns a best-of selection of indices based on quantitative networks;
     \item \option{topology} returns indices more abstractly describing network properties. Also CHECK details below!
     }
   }
  \item{level}{For which level should the level-specific indices be computed: \option{both} (default), \option{lower} or \option{higher}?}
  \item{weighted}{Logical; should the weighted average be computed for indices that are averaged across species (at the group level)? Defaults to TRUE.}
  \item{ISAmethod}{Method to use for calculating interaction strength (= dependence) asymmetry; original by \option{Bascompte} is yielding artefact results based only on the asymmetry of the web (as shown by example in Blüthgen et al. 2007 analytically in Blüthgen 2010) and should hence be avoided; \option{Bluethgen} (default) excludes singletons and corrects for low number of interactions (range -1 to 1).}
  \item{SAmethod}{How to aggregate d'-based specialisation values: mean of log-transformed dependencies (\option{log}) or Blüthgen's marginal totals-weighted mean (default); see Blüthgen et al. (2007).}
  \item{extinctmethod}{Specifies how species are removed from matrix: \option{random}, \option{degree} or \option{abundance} (partial matching). See \code{\link{second.extinct}} for details an option to predefine the sequence externally; idea from Memmott et al. (2004).}
  \item{nrep}{Number of replicates for the extinction sequence analysis.}
%  \item{plot.it.extinction}{logical; plot the extinction sequence graph; defaults to FALSE. }
%  \item{plot.it.dd}{logical; plot the degree distribution fits?; defaults to FALSE. }
  \item{CCfun}{Method to use when calculating the clustering coefficient. Originally proposed as mean of cluster coefficients for each species. Defaults to \option{median}, because cluster coefficients are strongly skewed. }
  \item{dist}{Distance metric to be used to calculate niche overlap. Any of \pkg{vegan}'s \code{vegdist}-metrics can be used; defaults to Horn's index, which is the recommendation of Krebs (1989). Binary percent niche overlap would be computed with \cr  \option{dist = "jaccard"}.}
  \item{normalise}{Logical; shall the C-score and togetherness metrics be normalised to a range of 0 to 1? Defaults to \code{TRUE}.}
  \item{empty.web}{Shall the empty columns and rows be deleted? Defaults to TRUE.}
  \item{logbase}{Shall the various diversity indices (linkage density, partner diversity, generality/vulnerability, interaction evenness) be calculated to the base of e (default) or 2? Log2 is the proposal for generality and vulnerability by Bersier et al. (2002), while Shannon uses ln. The choice of the base will not affect the results qualitatively, at most by a scaling factor. Note that for all these indices, we follow common practice and define 0 * log(0) = 0.}
  \item{intereven}{Shall all cells of the matrix be used to calculate the interaction evenness \cr
  (\option{intereven = "prod"})? Or, as given by Bersier et al. (2002) and Tylianakis et al. (2007), should only the realised links be used (\option{intereven = "sum"}; default)? Prod and sum refer to using the log of the \bold{product} of matrix dimensions (i.e. all cells) or the log of the \bold{sum} of non-zero cells (i.e. number of links)  as denominator in the evenness formula. See last paragraph of the details-section for views on these two options!}
  \item{H2_integer}{Logical; indicates whether values in web are integers. Passed on to \code{\link{H2fun}}; see there for details.}
  \item{fcweighted}{Logical; when computing ``functional complementarity'' sensu function \code{\link{fc}}, should the weights of the matrix be used. Defaults to TRUE, but original paper (Devoto et al. 2012) is based on FALSE.}
  \item{fcdist}{Distance measure to be used to compute functional complementarity through \code{\link{fc}}; any measure accepted by \code{\link{dist}} is acceptable.}
  \item{effective}{logical; should interaction evenness, Alatalo evenness, H2' and diversity be expressed as ``effective'' diversity? For Jost (2010, page 210) this seems to be the better decomposition of an index into a diversity and an evenness component. On a downside, H2' is then no longer ranged in [0,1], and ``diversity'' becomes the same as ``vulnerability'' and ``generality''. For these reasons, and backward compatibility, the default is \code{FALSE}.}
  \item{legacy}{Logical; should the old (pre-2.00) version of \code{networklevel} be used? To be backward compatible, the old \code{networklevel}-function is still available (\code{.networklevel}) and can be called by setting \option{legacy=TRUE}. This is only for the transition period until all papers in the making have been published (or binned). Index names and sometimes unclear focal level were downsides of the old implementation, which is now remedied. Thus, the use of \option{legacy=TRUE} and the direct call of \code{.networklevel} are strongly discouraged!}
  \item{mass.action.norm}{Logical; shall the spectral radius be normalised by mass action? Defaults to FALSE. See \code{\link{spectral.radius}} for details.}
}

\details{
For explanations of any of the indices computed for a level (i.e. those with HL and/or LL suffix), please see \code{\link{grouplevel}} for details.

This function implements a variety of the many (and still procreating) indices describing network topography. Some are embarrassingly simple and mere descriptors of a network's outer appearance (such as number of species in each trophic level or the number of links (= non-zero cells) in the web). Others are variations on Shannon's diversity index applied to within column or within rows. Only extinction slope is newly implemented here, and hence described in a bit more detail.

Currently, you \emph{cannot} get the \emph{qualitative} version of quantitative indices such as vulnerability!% by simply calling \code{networklevel} on \dQuote{binarised} data: \code{networklevel(Safariland>0, index="vulnerability")}. (Why you would want to do that, however, is currently beyond me.)%(Thanks to Jason Tylianakis for proposing to put this clarification into the help, and for Yvonne Fabian for correcting the clarification!)

Integers or continuous values - what are the quantities in quantitative webs? Some web metrics expect in their typical formulation that the entries in the web-matrix are integers - e.g. H2' is defined relative to minimum and maximum based on marginal totals. Blüthgen et al. (2006) use an algorithm assuming values can only be integers. If your quantities are not constrained to be integers, multiplication and rounding may or may not give consistent results, depending on rounding errors and the factor applied. Multiplication with high numbers such as 10 000 seems to be OK. For H2' a simplified calculation applicable to continuous numbers is available (by declaring option \option{H2_integer=FALSE} in \code{\link{H2fun}}). Note that values of H2' based on integers are not directly comparable to H2' based on continuous values (for sparse webs, H2'_continuous is much higher than H2'_integer). We tentatively think that other indices are hardly affected by non-integer values or by multiplication and rounding. Please let us know your experience.
}

\value{
The suffixes LL and HL refer to lower and higher level, respectively

Depending on the selected indices, some or all of the below (returned as vector if \dQuote{degree distribution} was not requested, otherwise as list):
	\item{connectance}{Realised proportion of possible links (Dunne et al. 2002): sum of links divided by number of cells in the matrix (= number of higher times number of lower trophic level species). This is the \emph{standardised number of species combinations} often used in co-occurrence analyses (Gotelli & Graves 1996)}
	\item{web asymmetry}{Balance between numbers in the two levels: positive values indicate more higher-trophic level species, negative more lower-trophic level species; implemented as \code{(ncol(web)-nrow(web))/sum(dim(web))}; web asymmetry is a null model for what one might expect in dependence asymmetry: see Blüthgen et al. (2007).}
	\item{links per species}{Mean number of links per species (qualitative): sum of links divided by number of species.}
	\item{number of compartments}{Compartments are sub-sets of the web which are not connected (through either higher or lower trophic level) to another compartment. Mathematically, they are Jordan blocks, but this implementation is rule-based (and fast). They are also nicely visualised in the \code{\link{visweb}} function.}
	\item{compartment diversity}{Shannon's diversity of compartment sizes (size = number of species from both levels); see Tylianakis et al. (2007).}
	\item{cluster coefficient}{The cluster coefficient for a network is the average cluster coefficients of its members, i.e. simply the number of realised links devided by the number of possible links. Introduced by Watts & Strogatz (1998) and described in Wikipedia under \url{https://en.wikipedia.org/w/index.php?title=Clustering_coefficient}. The cluster coefficient can be computed both for the entire network, as well as for each level (for the latter indicated by suffix HL or LL).}
	\item{nestedness}{Nestedness temperature of the matrix (0 means cold, i.e. high nestedness, 100 means hot, i.e. chaos). \code{networklevel} calls \code{\link[vegan]{nestedtemp}}! If you are interested in the different null models, please use the function \code{\link{nested}} or \code{\link[vegan]{nestedtemp}} directly.}
	\item{NODF}{Another index for nestedness, calling \code{\link[vegan]{nestednodf}}. High values indicate nestedness. According to the analysis of Almeida-Neto et al. (2008, 2010), NODF is more consistent and \dQuote{better} than usual measures of nestedness.}
	\item{weighted nestedness}{A nestedness version that considers interaction frequencies (and is hence weighted), proposed by Galeano et al. (2007) and implemented in \code{\link{wine}}. It ranges between 1 (perfect nestedness) and 0 (perfect chaos). Note that this is the OPPOSITE interpretation of nestedness temperature!}
	\item{weighted NODF}{Another quantitative (=weighted) index for nestedness, building on NODF (see \code{\link[vegan]{nestednodf}}). High values indicate nestedness. According to the analysis of Almeida-Neto et al. (2008, 2010), NODF is more consistent and \dQuote{better} than usual measures of nestedness.}
	\item{spectral radius}{Another nestedness index, proposed by Staniczenko et al. (2013), possibly correcting for "mass action". See \code{\link{spectral.radius}} and/or \code{\link{massaction}} for details.}
	\item{interaction strength asymmetry}{(selected using \option{index = "ISA"}). Explaining dependence asymmetry is also a measure of specialisation, across both trophic levels. Proposed by Bascompte et al. (2006) and critised and alterations proposed by Blüthgen et al. (2007). The latter also show that dependence asymmetry can be almost entirely explained by web asymmetry (see above). Positive values (only possible of \option{ISAmethod = "Bluethgen"}) indicate higher dependence in the higher trophic level. See function \code{\link{specieslevel}} and its index \option{interaction push/pull}, which quantifies the balance of affecting and being effected by other species. Similarly, index \option{strength} quantifies the average effect of each species on all its partners.}
	\item{specialisation asymmetry}{(selected by using \option{index="SA"}). Asymmetry (higher vs. lower trophic level) of specialisation now based on d' (see \code{\link{dfun}}), which is insensitive to the dimensions of the web. Again, two options of calculation are available: the one proposed by Blüthgen et al. (2007), where they weight the specialisation value for each species by its abundance (\option{SAmethod="Bluethgen"}) or where d'-values are log-transformed (arguing that d'-values are indeed log-normally distributed: \option{SAmethod="log"}). Since the mean d-value for the lower trophic level is subtracted from that of the higher, positive values indicate a higher specialisation of the higher trophic level.}
    \item{linkage density}{Marginal totals-weighted diversity of interactions per species (quantitative). Actually, this is computed as the average of vulnerability and generality (Bersier et al. 2002). Does not respond to setting \code{weighted=FALSE}.}
    \item{weighted connectance}{Linkage density divided by number of species in the network (Bersier et al. 2002). This will respond to whether non-interacting species (e.g. unparasitised hosts) are included or not!}
    \item{Fisher's alpha}{An alternative measure of interaction diversity (using \code{\link[vegan]{fisherfit}}).}
%         \item{interaction diversity LTL}{Mean diversity of predators per prey. LTL refers to Lower Trophic Level. LTL is analog to Shannon's diversity index for interactions (see below), but calculated for each species in the LTL, then averaged. Or, in other words, it is the log-version of vulnerability, but averaged without weighting.}
%         \item{interaction diversity HTL}{Mean diversity of prey per predator. HTL refers to Higher Trophic Level.}
    \item{interaction evenness}{Shannon's evenness for the web entries. Note that the two options are rather different. By definition, IE = H/Hmax; H = -sum(p.i.mat*log(p.i.mat)), where p.i.mat = matrix/sum(entries in matrix). This means, when calculating H, do we treat all possible links as species, and the interactions (cell values) as measure of their abundance? By definition, Hmax = ln(N). The key question is: What is the right value for N? Since we treat the matrix cells as species, it is (clearly?) the number of matrix cells, i.e. number of higher trophic level species x number of lower trophic level species. We think this logic justifies our default \code{"prod"}. However, others argue in favour of N=number of links. Please see note for our discussion on this point. Note that \option{effective=TRUE} will turn this into Jost (2019)'s ``effective evenness''.}
%What else?  Were we to use the interpretation of Tylianakis et al. (2007), then Hmax = ln(sum of entries in matrix). This means, we equate ``number of links'' = ``number of species''. That means, each link is a potential species. What should that mean? Why should that measure ``interactions evenness''? Why should we move from a view of ``cells are species'' when calculating H to a view of ``interactions are species'' when calculating N? To say the least, it doesn't seem consistent. 
    \item{Alatalo interaction evenness}{A different measure for web entry evenness, as proposed by Müller et al. (1999). }
    \item{Shannon diversity}{Shannon's diversity of interactions (i.e. network entries).}
    \item{H2}{H2' is a network-level measure of specialisation. It ranges between 0 (no specialisation) and 1 (complete specialisation). More specifically, H2' is a measure of discrimination, i.e. calculated in comparison to no specialisation (see \code{\link{H2fun}} for details. To avoid confusion of keys (apostrophe vs. accent), we call H2' only H2 here.}
    \item{others now to come:}{all other indices are returned as output from \code{\link{grouplevel}}. Please see there for details, we here only provide minimal listing.}
    \item{number of species}{}
    \item{mean number of shared partners}{in this level}
    \item{cluster coefficient}{for this level (same for both levels if \option{weighted=FALSE}).}
    \item{weighted cluster coefficient}{}
    \item{niche overlap}{Mean similarity in interaction pattern between species of the same level, calculated by default as Horn's index (\option{dist="horn"}).}
    \item{togetherness}{Mean number of co-occupancies across all species combinations.}
    \item{C score}{Mean (normalised) number of checkerboard combinations across all species.}
    \item{V ratio}{Variance-ratio of species numbers to individual numbers within species for that level.}
    \item{discrepancy}{Discrepancy as proposed by Brualdi & Sanderson (1999); see also \code{\link{discrepancy}} for details. However, \code{networklevel} actually calls \bold{vegan}'s \code{nesteddisc} to better handle ties.}
    \item{degree distribution}{See \code{\link{degreedistr}} for details and references.}
    \item{extinction slope}{Slope of the secondary extinction sequence in that level, following extermination of species in the other level.}
    \item{robustness}{Area below the \dQuote{secondary extinction} curve; see \code{\link{robustness}} for details. Corresponds to \dQuote{extinction slope}.}
    \item{functional complementarity}{for a given level.}
    %#! JFedit: changed from functional diversity
    \item{partner diversity}{(Weighted) mean Shannon diversity of the number of interactions for the species of that level. Choose \option{logbase=2} to change to a log2-based version.}
    %#! JFedit: effpartners and gen/vul combined
    %\item{effective partners}{(Weighted) mean of exp(Shannon diversity). Simply the Jost (2006)-recommended version of diversity, which happens to be identical to generality (higher level) or vulnerability (lower level), as defined in Bersier et al. (2002).}
    \item{generality/vulnerability}{(Weighted) mean effective number of LL species per HL species (generality; HL species per LL species for vulnerability), weighted by their marginal totals (row sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical to exp(\dQuote{partner diversity}, i.e., simply the Jost (2006)-recommended version of diversity.}
}

\note{
\bold{Sum or Prod: How to calculate interaction evenness?} I shall first put down my argument for ``prod'' and then Jason Tylianakis' arguments for ``sum''. 

Carsten: ``I do not want to defend a position I cannot hold against the flood of qualified criticism, and shall be happy to change the default to option ``sum'' (i.e. Jason's proposal). Nevertheless, I shall make a very brief attempt to defend my (and Nico's point of view). Imagine a completely different situation: I have ``counted'' birds in a landscape. From a more meticulous colleague I know that there are 27 bird species breeding at the moment, but on that two mornings that I went out, I could only hear 15. Now I want to calculate the Shannon diversity (and evenness) of birds in that landscape. The ``normal'' (in the sense of established) approach to use the data from my 15 species. But hold on: I KNOW there are more species out there. I don't know how many (i.e. there may be more than the 27 my colleague has found), but there are at least 27. If I only use the data from my 15 species, I will get a higher evenness value than when I also include the 12 zeros. My conclusion would be: I don't want to overestimate evenness only because I couldn't look long enough, thus I use all 27 values.''

Jason: ``I would disagree because what you ``know'' is based on your meticulous colleague's `sampling', which will also have its limits. If all you wanted was to know the total number of species there (assuming none have gone extinct), then what you propose is fine. However, the problem comes when you want to compare sites, and then sampling effort should be standardised. In most cases we know we don't have a full representation of the diversity (or food web) of an area, but we know for a given spatial or temporal sampling scale that one site differs from another in certain ways, and to me that is the most important. Anyway, it is all a question of scale and the precise question being asked. So what about making it an option in bipartite that you can either choose to divide by the realised links (give our 2007 paper as a ref, so people know it's comparable to that) or divide by the number of potential links, if that's the question people want to ask?''
There you go: it's your choice!

\bold{NA values:}
All error and warning messages are (or at least should be) suppressed! If your web returns and NA for some of the indices, this can be because the index cannot be computed. For example, if the web is full (i.e. no 0-cells), extinction slopes cannot be fitted (singularity of gradient). Check if you can expect the index to be computable! If it is, and networklevel doesn't do it, let me know.

\bold{Reducing computation time:}
Some indices require rather long computation times on large webs. If you want to increase the speed by omitting some indices, here a rough guide: Ask only for the indices you are interested in! Otherwise, here is the sequence of most time-consuming indices: 
    \enumerate{
	 \item The slowest function is related to extinction slopes and robustness. Excluding both makes the function faster. 
	 \item \option{weighted cluster coefficient} is also very time consuming (an exhaustive search for 4-loops in the one-mode projection of the network). Omitting it can dramatically boost speed.
	 \item Degree distributions are somewhat time consuming. 
	 \item Fisher's alpha is computed iteratively and hence time consuming. 
	 \item Nestedness and weighted nestedness are not the fastest of routines. 
	 \item Number (and diversity) of compartments calls a recursive and hence relatively slow algorithm. 
	 \item H2 and specialisation asymmetry require an iterative, heuristic search algorithm. Finally, excluding discrepancy can also moderately decrease computation time.
   }
   
\bold{Does the species sequence in the data matter?}
Obviously, it shouldn't, and for most indices it doesn't. However, particularly indices based on binary representation of \option{web} will have ties, where several species have the number of links. In this case, it does matter how the matrix is sorted before simulating extinctions or computing discrepancy. There is no (known) foolproof way to get this sequence "right" (see also "Details" of help for \code{\link[vegan]{nesteddisc}}). Re-running the same code with a shuffled network may thus yield (slightly) different values for (weighted) nestedness, togetherness, discrepancy, extinction slopes and robustness. (Thanks to Valentin Stefan for making us explicitly addressing this issue!)
}

\references{ 
Almeida-Neto, M., Loyola, R.D., Ulrich, W., Guimaraes, P., Guimaraes, Jr., P.R. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}, 1227--1239

Almeida-Neto, M. & Ulrich, W. (2011) A straightforward computational approach for measuring nestedness using quantitative matrices. \emph{Environmental Modelling & Software} \bold{26}, 173--178

Bascompte, J., Jordano, P. and Olesen, J. M. 2006 Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. (2002) Quantitative descriptors of food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N. (2010) Why network analysis is often disconnected from community ecology: A critique and an ecologist's guide. \emph{Basic and Applied Ecology} \bold{11}, 185--195

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen N. 2007 Specialization, constraints and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Burgos, E., H. Ceva, R.P.J. Perazzo, M. Devoto, D. Medan, M. Zimmermann, and A. Maria Delbue (2007) Why nestedness in mutualistic networks? \emph{Journal of Theoretical Biology} \bold{249}, 307--313

Corso G, de Araújo AIL, de Almeida AM (2008) A new nestedness estimator in community networks. \emph{arXiv} 0803.0007v1 [physics.bio-ph]

Devoto M., Bailey S., Craze P., and Memmott J. (2012) Understanding and planning ecological restoration of plant-pollinator networks. \emph{Ecology Letters} \bold{15}, 319--328. http://dx.doi.org/10.1111/j.1461-0248.2012.01740.x

Dormann, C.F., Fründ, J., Blüthgen, N., and Gruber, B. (2009) Indices, graphs and null models: analysing bipartite ecological networks. \emph{The Open Ecology Journal} \bold{2}, 7--24.

Dunne, J. A., R. J. Williams, and N. D. Martinez 2002. Food-web structure and network theory: the role of connectance and size. \emph{Proceedings of the National Academy of Science USA} \bold{99}, 12917--12922

Galeano, J., Pastor, J.M. and Iriondo, J.M. 2008. Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Gotelli, N. J., and G. R. Graves. 1996 Null Models in Ecology. Smithsonian Institution Press, Washington D.C.

Jost, L. 2010. The relation between evenness and diversity. \emph{Diversity} \bold{2}, 207--232. https://doi.org/10.3390/d2020207 

Krebs, C. J. 1989. Ecological Methodology. Harper Collins, New York.

Memmott, J., Waser, N. M. and Price M. V. 2004. Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999. The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Roberts, A. and Stone, L. 1990 Island-sharing by archipelago species. \emph{Oecologia} \bold{83}, 560--567

Rodríguez-Girónes M.A., and Santamaría L. 2006. A new algorithm to calculate the nestedness temperature of presence-absence matrices. \emph{Journal of Biogeography} \bold{33}, 924--935

Schluter, D. 1984. A variance test for detecting species associations, with some example applications. \emph{Ecology} \bold{65}, 998-1005.

Staniczenko, P., Kopp, J. & Allesina, S. (2013) The ghost of nestedness in ecological networks. \emph{Nat. Commun.} \bold{4}, 1391. \doi{10.1038/ncomms2422}

Stone, L. and Roberts, A. 1990. The checkerboard score and species distributions. \emph{Oecologia} \bold{85}, 74--79.

Stone, L. and Roberts, A. 1992. Competitive exclusion, or species aggregation? An aid in deciding. \emph{Oecologia} \bold{91}, 419--424

Tylianakis, J. M., Tscharntke, T. and Lewis, O.T. 2007. Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Watts, D. J. and Strogatz, S. 1998. Collective dynamics of \sQuote{small-world} networks. \emph{Nature} \bold{393}, 440--442
}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\seealso{ Some functions are implemented separately: \code{\link{H2fun}}, \code{\link{second.extinct}}, \code{\link{degreedistr}}, \code{\link{C.score}} and \code{\link{V.ratio}}
}

\examples{
\dontrun{
data(Safariland)
networklevel(Safariland)
networklevel(Safariland, index="ALLBUTDD") #excludes degree distribution fits
}
}

\keyword{ package }
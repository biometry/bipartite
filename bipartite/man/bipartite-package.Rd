\encoding{UTF-8}
\name{bipartite-package}

\alias{bipartite-package}
\alias{bipartite}

\docType{package}

\title{
Analysis of bipartite ecological webs
}

\description{
Bipartite provides functions to visualise webs and calculate a series of indices commonly used to describe pattern in (ecological) networks, a.k.a. webs. It focusses on webs consisting of only two levels, e.g. pollinator-visitation or predator-prey webs. Visualisation is important to get an idea of what we are actually looking at, while the indices summarise different aspects of the webs topology.}

\details{
%We only had three types of bipartite webs in mind when writing this package: seed-disperser, plant-pollinator and host-parasitoid systems. In how far it makes sense to use these functionalities for other systems (or indeed for these systems) lies in the hands of the user. Please refer to the literature cited for details on the theory behind the indices.
%
Input for most analyses is an interaction matrix of m species from one group (\dQuote{higher}) with n species from another group (\dQuote{lower}), i.e. a n x m matrix, where higher level species are in columns, lower level species in rows. Column and row names can be provided. This is fundamentally different from \dQuote{one-mode} webs, which are organised as k x k matrix, i.e. one group of species only, in which each species can interact with each other. Such a format is incompatible with the functions we provide here. (Note, however, that functions \code{\link{as.one.mode}} and \code{\link{web2edges}} are convenience functions to morph bipartite networks into one-mode webs. Furthermore, some indices build on one-mode networks and are called from bipartite.)

Before you start with the network, you have to get the data into the right shape. The function \code{\link{frame2webs}} aims to facilitate this process. Arranging a web, e.g. by size, is supported by \code{\link{sortweb}}.

The typical first step is to \bold{visualise} the network. Two functions are on offer here: one (\code{\link{visweb}}) simply plots the matrix in colours depicting the strength of an interaction and options for re-arranging columns and rows (e.g. to identify compartments or nesting). The other function (\code{\link{plotweb}}) plots the actual web with participants (as two rows of rectangles) connected by lines (proportional to interaction strength). Both can be customised by many options.

The second step is to \bold{calculate indices} describing network topography. There are \bold{three} different levels this can be achieved at: the entire web (using function \code{\link{networklevel}}), at the level of each group (also using function \code{\link{networklevel}}) or the individual species (using function \code{\link{specieslevel}}). Most other functions in the package are helpers, although some can be called on their own and return the respective result (\code{\link{dfun}}, \code{\link{H2fun}} and \code{\link{second.extinct}} with \code{\link{slope.bipartite}}).

The third step is to \bold{compare results to null models}. Many interaction matrices are very incomplete snapshots of the true underlying network (e.g. a one-week sampling of a pollination network on a patch of 4 x 4 meters). As a consequence, many species were rarely observed, many are singletons (only one recording). To make analyses comparable across networks with different sampling intensity and number of species per group, we need a common yardstick. We suggest that users should use a null model, i.e. an algorithm that randomises entries while constraining some web properties (such as dimensions, marginal totals or connectance). The function \code{\link{nullmodel}} provides a few such null models, but this is a wide field of research and we make no recommendations (actually, we do: see Dormann et al. 2009 and Dormann 2011, both shipping in the doc-folder of this package). You can also simulate networks using \code{\link{genweb}} or \code{\link{null.distr}}.

Finally, bipartite comes with 23 quantitative pollination network data sets taken from the NCEAS interaction webs data base (use \code{data(package="bipartite")} to show their names) and it has a few miscellaneous functions looking at some special features of bipartite networks (such as modularity: \code{\link{computeModules}} or apparent competition: \code{\link{PAC}}).

\bold{Speed}: The code of bipartite is almost exclusively written in R. You can increase the speed a bit (by 30 to 50 \%, depending on the functions you use) by compiling functions on-the-fly. To do so, you need to load the \pkg{compiler} package and type: \code{enableJIT(3)}. The first time you call a function, it will be compiled to bytecode (just-in-time: jit), which takes a few seconds, but the second call will be substantially faster than without compilation. In the few tests we have run, this improvement was NOT substantial (i.e. a few tens of percent), indicating, I guess, that our R code wasn't too bad. See \pkg{compiler}-help files or \url{http://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/} for details.

See help pages for details and examples.

\tabular{ll}{
Package: \tab bipartite\cr
Type: \tab Package\cr
Version: \tab 2.06\cr
Date: \tab 2016-04-18\cr
License: \tab GPL \cr
}
}

\section{versionlog}{
Please see help page \code{\link{versionlog}} for all changes and updates prior to version 2.00. This page will only list most recent changes. 

\itemize{
  \item 2.07 (release date: 18-Apr-2016)
   \describe{ % begin describe 2.07
   \item{Some explanation addded to \code{\link{czvalues}},}{ where a z-value of NA is returned if a species is alone (in its trophic level) in a module. This is due to the way z-values are computed, and not a bug.}
   \item{Function \code{\link{nestedcontribution}}}{was not exported in the namespace. Fixed. Thanks to various people reporting this.}	
   } % end describe 2.07

  \item 2.06 (release date: 29-Sep-2015)
   \describe{ % begin describe 2.06
   \item{Bug fix in \code{\link{C.score}},}{ which did not compute the maximum number of possible checkerboards correctly, and hence let the normalised C-score to be incorrect. Now it uses a brute-force approach, which works fine but takes its time.}	
   \item{Function \code{\link{nestedcontribution}}}{was not exported (i.e. not listed in the namespace file). Fixed. Thanks to Wesley Dátillo for reporting.}
   \item{Help page of \code{\link{specieslevel}}}{now correctly described a species' degree as sum of its links. Thanks to Bernhard Hoiß for the correction!}
   \item{C++-warnings addressed:}{outcommented some unused variables in dendro.h and removed some fprintf-warnings in bmn5.cc}
   \item{Little bug fix in \code{\link{vaznull}}:}{Threw an error when matrix was without 0s. Thanks to Thais Zanata for reporting.}
   } % end describe 2.06
   
  \item 2.05 (release date: 24-Nov-2014)
   \describe{ % begin describe 2.05
   \item{New function \code{\link{nestedcontribution}}}{which computes the contribution of each species to the overall nestedness, based on Bascompte et al. 2003 and as used by Saavedra et al. 2011. Many thanks to Daniel Stouffer for contributing this function!}	
  \item{New function \code{\link{mgen}}:}{this function is based on R-code written by Diego Vázquez (many thanks for sending the code), with a bit of brushing up options by CFD. The function takes a probability matrix generated by whatever mechanism and builds a null model from it. This is a niffty little idea, making null modelling concerned with generating ideas on what makes an interaction probable and leaving the step of producing and integer-network of simulated interactions to this function.}
   \item{minor fixes in \code{\link{networklevel}}}{``weighted connectance'' was only returned when ``linkage density'' was in ``index'' call; now also available on its own. Also sligthly revised the help file.}
   \item{\code{\link{nested}} with option \option{weighted NODF} called the unsorted version of this function,}{while calling the same index in \code{\link{networklevel}} called the sorted. This is not nice (although not strictly wrong). Now both call the sorted version and users have to directly invoke \code{nestednodf} for the unsorted option. Many thanks to Julian Resasco for reporting!}
  \item{Changes to the help page of \code{\link{vaznull}}:}{I (CFD) misread the original paper introducing this null model and hence assumed that\code{\link{vaznull}} would constrain marginal totals \bold{and} connectance. However, this was not intended in Diego Vázquez original implementation and never stated anywhere (except in the help pages of this function here in bipartite). Hence, the help pages were changed to now reflect both intention and actual action of this function. This also means that currently only one null model with constrained marginal totals \bold{and} connectance is available: \code{\link{swap.web}}. Many thanks to Diego for clearing this up!}
   \item{Some example code had to be re-written}{to adapt to the upcoming/new structure of \pkg{vegan}, which got rid of function \code{commsimulator} (replaced by \code{simulate}). Many thanks to Jari Oksanen for informing me about this!}
   \item{Added an error message}{to function \code{\link{second.extinct}} for the case that a user wants to provide an extinction sequence for both trophic levels. There is no obvious way to simulate this across the two groups, and hence it is not implemented. Also added error messages for non-matching vector/web dimensions and alike.}
   } % end describe 2.05

  \item 2.04 (release date: 25-Mar-2014)
   \describe{ % begin describe 2.04
  \item{R-C++-communication bug fixed in \code{\link{computeModules}}:}{This bug has been a constant thorn in my side. Somehow the C-code behind \code{computeModules} could only be called once. On second call, it returned an error because somehow it kept some old files in memory. So far, I used a work-around (unloading and re-loading the dynamic library), which only worked on Windows and Mac. I still don't fully understand it, but thanks to Tobias Hegemann (whom I paid for being more competent than myself) we now have a function running bug-free on all platforms. (Deep sigh of relief.)}
  \item{The call of index ``functional complementarity'' through \code{\link{networklevel}} did not work.}{Fixed this legacy issue, which was due to a confusion created by the index' earlier name of ``functional diversity''.}
  \item{Help page to \code{\link{specieslevel}} gave incomplete name for one index:}{Should be \option{interaction push pull}; also the function itself had the ``push pull''-bit messed up. Thanks to Natacha Chacoff for reporting!}
  \item{Sequence of indices differed between lower and higher level. (Fixed.)}{Both should be the same and should fit the description in the help file. Thanks to Jimmy O'Donnell for reporting!}
   } % end describe 2.04

  \item 2.03 (release date: 15-Jan-2014)
   \describe{ % begin describe 2.03
  \item{Some ghost text led to conflicts with the updated package checking.}{Ghost text deleted. Thanks to Brian Ripley of the R-Team and CRAN for not only reporting the issue but also pointing to its solution!}
  \item{Option \option{empty.web} added to \code{\link{specieslevel}}:}{Similar to the argument in \code{\link{networklevel}}; non-interacting species from the network were always excluded so far; new option \option{FALSE} not fully tested yet.}
  \item{Minor bug fix in \code{\link{specieslevel}}:}{ ``pollination support index'' returned ``PSI''; ``PDI'' now referenced correctly as ``paired differences index''.}
  \item{Simplification in \code{\link{grouplevel}} and correspondingly in \code{\link{networklevel}}:}{Previously, \code{index="generality"} or \code{"vulnerability"} was identical to \code{"effective partners"} with option \code{weighted=TRUE}, but different for \code{weighted=FALSE} (to which only \code{"effective partners"} responded). We reduced this to one index called "generality" or "vulnerability" (depending on the focal group), but which will now give the non-weighted mean if option \code{weighted=FALSE}. It can still be called by "effective partners" for backward compatibility.}
  \item{Function \code{\link{grouplevel}} used \code{fd} wrongly!}{Instead of returning the value for rows, it returned the functional diversity for columns (and vice versa). We also used the opportunity to rename the index to its correct name: ``functional complementarity'' and the function to \code{\link{fc}}. Help pages for \code{\link{fc}} and \code{\link{grouplevel}} were adapted accordingly. Thanks to Mariano Devoto for pointing out this mistake!}
	\item{New index ``weighted connectance'' in function \code{\link{networklevel}}:}{This index is simply computed as linkage density divided by number of species in the network. Note that using \option{empty.web=TRUE} will affect this value (which is intended). Thanks to Becky Morris for suggesting to add this index here.}
	\item{Help page for function \code{\link{PDI}} corrected.}{Thanks to Timothy Poisot for reporting some issues in the help page.}
   } % end describe 2.03
  \item 2.02 (release date: 30-Sep-2013)
   \describe{ % begin describe 2.02
   	\item{Glitch fixed in \code{\link{grouplevel}} (thus also affecting \code{networklevel}).}{Networks with only one species in one of the two levels resulted in errors, rather than simply return NA for C-score and secondary extinction computation. Thanks to whoever it was for reporting (at the INTECOL workshop).}
	\item{Minor bug fixes in  \code{\link{specieslevel}}:}{Gave error messages for closeness and betweenness if the network had no shortest path. Now returns a warning and NAs instead. Reported: JF.}
	\item{Minor bux fix in  \code{\link{networklevel}}:}{Failed to work when an index was listed twice in the function call. Reported: JF.}
	\item{New function \code{\link{r2dexternal}}:}{This function is a null model algorithm like Patefields (\code{r2dtable}, but it excepts externally measured abundances to compute the null model-expectation. Experimental.}
	\item{Memory leak in \code{\link{computeModules}} fixed.}{Because some object was not deleted, memory consumption of this function shot through the roof (with time). Since R has a somewhat weird way of handling memory, I think that also subsequent operations were slower (because the dynamically expanded memory is not being shrunken again, which is a problem if you use the hard drive as RAM). Thanks to Florian Hartig for investing the time to fix it!}
   } % end describe 2.02
  \item 2.01 (release date: 28-Jun-2013)
  This release features smoothing of various glitches that were introduced when we cleaned up the code for version 2.00.
   \describe{ % begin describe 2.01
   	\item{New index for \code{\link{specieslevel}}:}{Computes the nestedness rank (as proposed by Alarcon et al. 2008). Can also be employed directly using the \bold{new function} \code{\link{nestedrank}} with options for weighting for number of interactions per link, normalising the rank and different method to compute the nestedness-arranged matrix.}
  	\item{Polishing \code{\link{specieslevel}}:}{Now returns an error message if the index is not recognised, instead of an empty list.}
	\item{Function \code{\link{plotweb}}}{received an option to plot additional individuals of a species in different ways. For a host-parasitoid network, some hosts are not parasitised. This data vector can now be interpreted in two ways, making the plotting function a bit more flexible.} %Thanks to Jochen Fründ for implementing.
	\item{Function \code{\link{degreedistr}}}{can now be invoked for each level separately. Also arguments can be passed to the plotting options.}
	\item{New data set \code{\link{junker2013}}:}{a nice and large pollination network. Thanks to Robert Junker for providing this data set!}
	\item{Fixed computation of secondary extinction slopes for both levels simultaneously}{for random extinction sequences. This was so far not possible, because the function did not combine extinction sequences of different lengths. This was simply an oversight, reported by Richard Lance. (Thanks!)}
   } % end describe 2.01
  
  \item 2.00 (release date: 15-Mar-2013)
 A new version number usually indicates substantial changes. In this case, we have re-named and re-grouped some of the output of \code{\link{networklevel}} and \code{\link{specieslevel}} for greater consistency and transparency. Beware! Running the same functions now (2.00 and up) will yield different results to <2.00 (because the same values are now in a different sequence).
  
We also started carefully renaming indices and re-writing help files. The main reason is that we started this work thinking of pollination networks. Over time, however, other types of ecological networks came into focus, and now also non-ecological networks are on the table. Thus, we started (and shall continue) referring to lower and higher levels, rather than plant and pollinators, hosts and predators or even trophic levels. Thus, in our emerging nomenclature the two levels are referred to as \dQuote{groups} (their members remain \dQuote{species} interacting with their \dQuote{partners} in the other group). 
 
Please read (or at least skim) the help pages before using a function of version 2.00 for the first time. 
  
In function \code{\link{specieslevel}} indices can now be computed for levels separately (or together). Few user-visible changes, but complete re-structuring under the hood. Option \option{species number} was moved to \code{\link{grouplevel}} as \option{number of species}.
  
In the new function \code{\link{grouplevel}} we collected all indices that can be computed for each of the two groups (i.e. trophic or other levels). Indices can be computed for each group separately or for both simultaneously. All group-level indices are also accessible through \code{\link{networklevel}}!

In the new function \code{\link{linklevel}} we collected all indices that can be computed for each cell of the bipartite matrix. Currently, there are few such indices, however.

In function \code{\link{networklevel}} we dropped the plotting options. Users wanting to plot degree distributions or extinction slopes are encouraged to use the functions \code{\link{degreedistr}} and \code{\link{slope.bipartite}}, respectively.

Furthermore, due to licensing issues, we copy-pasted several functions from the package \bold{tnet}, created and maintained by Tore Opsahl, to bipartite. We have so far called these functions from tnet, but only recently did R start to enforce license compatibility, which caused this step (bipartite being GPL and tnet being CC by-NC 3.0). We are really very grateful to Tore for allowing us to include the following functions: \code{\link{as.tnet}}, \code{\link{betweenness_w}}, \code{\link{closeness_w}}, \code{\link{clustering_tm}}, \code{\link{distance_w}}, \code{\link{symmetrise_w}}, \code{\link{tnet_igraph}}.
 
  Here a more detailed list of changes:
  \describe{
  	\item{\code{\link{networklevel}}}{
		\itemize{
			\item Function call and output now more consistent in naming and sequence. When higher and lower level indices are given (e.g. extinction slopes, number of shared partners), the first will always be the one referring to the property of the \emph{lower} level. From a pollinator network perspective, the first value in such a pair describes a plant-level index, the second a pollinator-level index.
			\item Indices \option{mean interaction diversity} dropped from \code{\link{networklevel}}. We found no reference to this metric and saw little use for it. It is very similar to vulnerability/generality and can easily be computed from the output of \code{\link{specieslevel}} as \code{mean(specieslevel(web, index="diversity"))}.
			\item  Now also accepts non-integer values as input. The argument \option{H2_integer} will then automatically be set to FALSE. Will return NA for those indices that cannot be computed (e.g. Fisher's alpha). As a knock-on effect, \code{\link{H2fun}} had to be slightly adapted to round to machine precision when searching for H2min. (A somewhat technical detail, but making \code{H2fun} getting caught sometimes.)
		} % end itemize
	}
	\item{New function \code{\link{grouplevel}}}{in which we collected indices that can be computed for each of the two groups (i.e. trophic or other levels). Indices can be computed for each group separately or for both simultaneously. All group-level indices are also accessible through \code{\link{networklevel}}!}
	\item{New function \code{\link{linklevel}}}{in which we collect indices that can be computed for each cell of the bipartite matrix.}
%	\item{Index \option{ISA} = \option{interaction strength asymmetry} = \option{dependence asymmetry}} dropped from \code{\link{networklevel}}{We found no study that constructively used this metric and saw little use for it. It can easily be computed from the output of \code{\link{specieslevel}} in two lines of code: \code{out <- specieslevel(web, index="dependence"); mean(abs(out[[1]][[1]]-out[[2]][[1]])/pmax(out[[1]][[1]], out[[2]][[1]]), na.rm=TRUE)}. }
	\item{New option to \code{\link{PDI}}:}{\option{normalise=FALSE} offers the option of using the index as originally proposed, although we prefer to use TRUE and made this the default. }
	\item{Corrected network \code{\link{bezerra2009}}.}{Network was actually the transpose of the correct network and hence wrongly had plant species as columns.}
	\item{New function \code{\link{endpoint}}}{computes end-point degrees following Barrat et al. (2004); one of the indices computed at \code{\link{linklevel}}.}
	\item{New function \code{\link{frame2webs}}}{helps organising data into one or more webs.}
	\item{New function \code{\link{webs2array}}}{helps organising webs into one array.}
	\item{Function \code{\link{specieslevel}}}{gained two new indices (thanks to Jochen Fründ): \option{proportional} \option{similarity} and \option{proportional generality}. See help page of that function for details.}
	\item{New function \code{\link{npartite}}}{Experimental function to analyse more-than-2-level networks.}
	\item{\code{\link{visweb}}}{now obeys the label size to make sure labels are always in the plotting area. Thanks to Zachary Grinspan %(no, he's not a character from Harry Potter, but he does have a sense of humour) 
	for drawing our attention to this issue.}
	\item{Little bug fix in \code{\link{second.extinct}}}{Function failed for argument \option{participant="both"} because I filled the extinction sequence with the wrong number of 0s (to achieve always the same dimensionality of results in repeated runs). Thanks to Carine Emer for reporting!}
	\item{\code{\link{specieslevel}}}{failed to work for non-matrix data (i.e. \code{data.frames}). It now coerces \code{data.frames} to \code{matrix} as a first step and hence should work also on \code{data.frame}s. Thanks to Marina Wolowski for drawing our attention to this problem.}
	\item{Minor bug fix in \code{\link{dfun}}:}{When external abundances were provided with a 0 in it, \code{dfun} could throw up \code{Inf}-values. Reported by Indrani Singh and fixed by Jochen Fründ.}
	\item{Settings for functions called by \code{\link{nested}}}{are now enshrined in stone. The initial reason was to set only the default for one function (\code{\link{nestedness}}) to a faster setting (\option{null.models=FALSE}), but then I decided to restrict all settings to the defaults of the functions called (except for this one option).}
	\item{Bug fix for the rarely used function \code{\link{null.t.test}}:}{Did not work if only one index was given.}

} % end of describe 2.00
} % end of versionlog's itemize

} % end of section versionlog

\author{
Carsten F. Dormann, Jochen Fründ and Bernd Gruber, with additional code from many others (referred to in the respective help file), noticeably from Tore Opsahl's tnet package.

Maintainer: Carsten Dormann \email{carsten.dormann@biom.uni-freiburg.de}

}
\references{
Alarcon, R., Waser, N.M. and Ollerton, J. 2008. Year-to-year variation in the topology of a plant-pollinator interaction network. \emph{Oikos} \bold{117}, 1796--1807

Almeida-Neto, M. and Ulrich, W. (2011) A straightforward computational approach for measuring nestedness using quantitative matrices. \emph{Environmental Modelling & Software}, \bold{26}, 173--178

Bascompte, J., Jordano, P. and Olesen, J. M. (2006) Asymmetric coevolutionary networks facilitate biodiversity maintenance. \emph{Science} \bold{312}, 431--433

Bersier, L. F., Banasek-Richter, C. and Cattin, M. F. (2002) Quantitative descriptors of food-web matrices. \emph{Ecology} \bold{83}, 2394--2407

Blüthgen, N., Menzel, F. and Blüthgen, N. (2006) Measuring specialization in species interaction networks. \emph{BMC Ecology} \bold{6}, 12

Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B. and Blüthgen, N. (2007) Specialization, constraints, and conflicting interests in mutualistic networks. \emph{Current Biology} \bold{17}, 1--6

Corso G., de Araújo A.I.L. and de Almeida A.M. (2008) A new nestedness estimator in community networks. \emph{arXiv}, 0803.0007v1 [physics.bio-ph]

Dalsgaard, B., A. M. Martín González, J. M. Olesen, A. Timmermann, L. H. Andersen, and J. Ollerton. (2008) Pollination networks and functional specialization: a test using Lesser Antillean plant-hummingbird assemblages. \emph{Oikos} \bold{117}, 789--793

Devoto M., Bailey S., Craze P. & Memmott J. (2012) Understanding and planning ecological restoration of plant-pollinator networks.
\emph{Ecology Letters} \bold{15}, 319--328

Dormann, C.F., Fründ, J., Blüthgen, N., and Gruber, B. (2009) Indices, graphs and null models: analysing bipartite ecological networks. \emph{The Open Ecology Journal} \bold{2}, 7--24

Dormann, C.F. (2011) How to be a specialist? Quantifying specialisation in pollination networks. \emph{Network Biology} \bold{1}, 1--20

Galeano J., Pastor J.M. and Iriondo J.M. (2008) Weighted-Interaction Nestedness Estimator (WINE): A new estimator to calculate over frequency matrices. \emph{arXiv} 0808.3397v1 [physics.bio-ph]

Martín Gonzáles, A.M., Dalsgaard, B. and Olesen, J.M. (2009) Centrality measures and the importance of generalist species in pollination networks. \emph{Ecological Complexity}, \bold{7}, 36--43

Memmott, J., Waser, N. M. and Price, M. V. (2004) Tolerance of pollination networks to species extinctions. \emph{Proceedings of the Royal Society B} \bold{271}, 2605--2611

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. (2004) Experimental evidence for apparent competition in a tropical forest food web. \emph{Nature} \bold{428}, 310--313

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. (2005) Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449--462.

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. (1999) The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

Poisot, T., Lepennetier, G., Martinez, E., Ramsayer, J., and Hochberg, M.E. (2011a) Resource availability affects the structure of a natural bacteria-bacteriophage community. \emph{Biology Letters} \bold{7}, 201--204

Poisot, T., Bever, J.D., Nemri, A., Thrall, P.H., and Hochberg, M.E. (2011b) A conceptual framework for the evolution of ecological specialisation. \emph{Ecology Letters} \bold{14}, 841--851

Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205

Vázquez, D. P. and Aizen, M. A. (2004) Asymmetric specialization: A pervasive feature of plant-pollinator interactions. \emph{Ecology} \bold{85}, 1251--1257

Vázquez, D.P., Chacoff, N.,P. and  Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology} \bold{90}, 2039--2046.
}

\keyword{ package }

\examples{
\dontrun{
data(Safariland)
plotweb(Safariland)
visweb(Safariland)
networklevel(Safariland)
specieslevel(Safariland)
}
}
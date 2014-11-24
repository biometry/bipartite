\encoding{UTF-8}
\name{versionlog}
\alias{versionlog}
\alias{versionlog}

\title{Log of bipartite versions and changes}
\description{Log tracking changes, bug fixes and new functions in \pkg{bipartite}}

\usage{
    versionlog()
}

\section{versionlog}{
 \itemize{
  \item 1.18 (release date: 06-Sep-2012)
  \describe{
	\item{New function and index \code{\link{fd}}:}{\code{\link{fd}} computes the functional diversity for the rows of a web as a measure of ecological niche complementarity (Devoto et al. 2012). Function written and provided by Mariano Devoto. \code{\link{fd}} is also included in the computation of \code{\link{networklevel}}.}  
	\item{New function \code{\link{czvalues}}}{to compute c- and z-values of the higher trophic level for modular networks. Requires a successful identification of modules through \code{\link{computeModules}}. These can then be used to identify connector species in a network.}
	\item{Speed-up of \code{\link{plotweb}}:}{Thanks to Dirk Raetzel, whose smart improvement made this function orders of magnitude faster! (Amazing how a total of 7 lines of code in a total of over 300 for this function can make such a vast improvement!)}
	\item{Speed-up of \code{\link{compart}}:}{Thanks to Juan M. Barreneche, whose smart improvement made this function orders of magnitude faster while using dramatically less memory.}
	\item{New option for simulating extinctions in \code{\link{second.extinct}} (and necessarily \code{\link{extinction}}):}{Use \option{method="external"} to impose a specific sequence of exterminations onto the function. Provide this sequence as an index vector (e.g. \code{c(4,2,3,1)}) through the argument \option{ext.row} or \option{ext.col} (for lower and higher species level, respectively). Requested/proposed by Matt Koski. (Implementation is far from elegant and will be cleaned up in a future version.)}
	\item{Modification of \code{\link{visweb}}:}{The user is now able to use user-defined colours (def.col) for interaction ranks (\option{square="d"}). Thereby the number of interactions can define the colour of the cell in the plot. Makes visweb-plot prettier.}
	\item{Important bug fix in \code{\link{robustness}} and \code{\link{second.extinct}}:}{\code{\link{robustness}} sometimes interpreted wrongly which trophic level was subject to simulated extinction. That was easy to fix. En route, however, I also discovered unexpected behaviour in the \code{\link{second.extinct}} function: it gobbled up the last species to be removed from the web. I do not think that extinction slope estimates are affected much, but robustness results may well be exhibiting a bias towards 0.5! Thanks to Silvia Santamaria for reporting this!}
	\item{Little bug fix in functions \code{\link{BC}} and \code{\link{CC}}:}{Functions always used unweighted networks for computing betweenness and closeness. The default is now to use weighted information, an option which can be switched off. Thanks to Michael Pocock for reporting and suggesting a fix! (Notice the help file for these functions and the large effect that the choice of settings can have on the results!)}
	\item{Little bug fix in function \code{\link{computeModules}}:}{Failed to work for networks with empty columns or rows. These will now be deleted (with a warning) before the analysis. Thanks to Julien Renoult for reporting.}
	\item{Little bug fix in function \code{\link{plotweb}}:}{Failed to work for networks with only one column or row. Thanks to Zhijiang Wang for reporting.}
	\item{Minor bug fix in \code{\link{networklevel}}:}{Caused an error if a selection of indices was provided.}
	\item{Corrections in the help of \code{\link{networklevel}}:}{Web asymmetry was explained the wrong way around (thanks to Bernhard Hoiß for reporting!) and there was some confusion in the help for whether quantitative or weighted-quantitative indices were implemented (Yvonne Fabian for correcting the clarification!).}
  }

  \item 1.17 (release date: 20-AUG-2011)
       \describe{
             \item{New specialisation index at the species level in function \code{\link{PDI}}:}{Computes the Paired Differences Index proposed by Poisot et al. (2011a,b). Also automatically now included in \code{\link{specieslevel}}.}
  	\item{Bug fix in \code{\link{specieslevel}}:}{Computation of betweenness was causing an error when the network was containing a single-link compartment. Thanks to Marco Mello for spotting!}
	\item{Another PDF is included (Dormann 2011):}{This paper compares specialisation indices, based on the function \code{\link{specieslevel}}. All source code to reproduce these results is included. It can serve as a guide to the \code{\link{specieslevel}} function and its interpretation.}
      }

  \item 1.16 (release date: 29-May-2011)
     \describe{
  	\item{Bug fix and extended documentation in \code{\link{dfun}}:}{1. Sometimes d.min was larger than d.raw yielding negative d'-values. Since the search for d.min is heuristic, it can sometimes fail to yield lower values than d.raw. In these (rare) cases, d.min is now set to d.raw. 2. Documentation for a specific (but not uncommon or unreasonable) constellation added: When independent abundances were provided, the empty rows/columns are purposefully \emph{not} removed from the web (because they now still contain information). Logically (and as implemented), this leads to d-values for these species of NA. This makes sense: the pollinator, say, has never been observed on any of the flowers, so how can we quantify its specialisation? While reasonable, it was undocumented behaviour. Thanks to Ana María Martín Gonzáles for drawing our attention to it!} 
	\item{Speed improvement in \code{\link{shuffle.web}}:}{Thanks to Paul Rabie, a faster implementation of the \code{shuffle.web} null model is available. Choosing the option \option{legacy=FALSE} will deliver the null models about 3 times faster.}
     }

  \item 1.15 (release date: 05-Mar-2011)
  \describe{
	\item{Bug fix in \code{\link{swap.web}}:}{Occasionally, an interaction too many was swapped, leading to a connectance slightly too low compared to the original matrix. Spotted and reported by Sven Hanoteaux. (For ``normal'' webs, this bug has luckily caused only very small discrepancies between the realised and the intended null web. We thus regard this error to be (luckily) of virtually no influence on the results reported in our Open Ecology Journal paper. Please also note that the help of \code{swap.web} recommends using \code{\link{vaznull}}.)}
	\item{Bug fix in \code{\link{shuffle.web}}:}{Some interactions were lost when there were more rows than columns. Spotted and reported by Sven Hanoteaux. (Please note that while a shameful mistake, it luckily does not affect the results reported in our Open Ecology Journal paper, where all webs analysed were asymmetric ``in the right way''.)}
  	\item{Minor bug fix in \code{\link{networklevel}}:}{No output was provided when only index ``Fisher alpha'' or ``mean interaction diversity'' was selected. Thanks to Sven Hanoteaux for reporting and fixing!}
  }	
	
  \item 1.14 (release date: 21-Dez-2010)
  \describe{
  	\item{New function \code{\link{plotPAC}}:}{A plotting function to go with the \code{\link{PAC}}-function, providing a nice circularly arranged bubble plot, which has been copied (in style) from Morris et al. (2005). Various colour options are available.}
	\item{New index ``weighted NODF'' included in \code{\link{networklevel}} and \code{\link{nested}}:}{Gustavo Carvallo has updated the existing function \code{nestednodf} in \pkg{vegan} to compute a weighted (quantitative) version as proposed by Ulrich & Almeida Neto (2011 Env. Mod. Soft). This is now linked into \pkg{bipartite}. }
	%Notice that weighted nestednodf requires vegan_1.18-14, which can be installed from R-forge until it is released on CRAN.
  	\item{New data set:}{bezerra2009 is a``full'' flower-visitation network of oil-collecting bees in Brazil. See help for details.}
	\item{Minor bug fix in \code{\link{specieslevel}}:}{The newly (1.13) added functionality for index ``weighted closeness'' did not work for webs without column and row names (such as null models). Also ``weighted betweenness'' was occasionally playing up.}
	\item{Minor bug fix in \code{\link{specieslevel}}:}{Betweenness and friends cannot be computed for fully compartmented webs (i.e. where each species interacts only with its own partner in the other trophic level). The function then failed instead of returning NA. Thanks to Nadine Sandau for reporting.}
	\item{Minor bug fix in \code{\link{dfun}}:}{When providing independent abundances, the computation of dmin was sometimes incorrect (particularly for small webs) due to a typo. Fixed by Jochen Fründ.}
  	\item{Minor time-saving improvement in \code{\link{degreedistr}}:}{When provided with 3 or less point to compute the degree-distribution fit from now immediately returns NA-filled matrices. Saves the time of going through all starting values and still failing.}
	\item{Change of name of function argument in \code{\link{web2edges}}:}{The former \option{quantitative} did not really reflect what it does, so it was changed to \option{weight.column}.}
	\item{Minor improvements to \code{\link{networklevel}}:}{Robustness and weighted nestedness were wrapped in \code{try}-functions (for very small networks). }
    \item{New functionality to find modules in networks.}{This is still experimental and proper references etc. will be added after testing! Please do not use for production yet (at least not without consultation). Thanks to Rouven Strauß for developing this.} 
}


  \item 1.13 (release date: 29-Sep-2010)
  \describe{
    \item{Option added to \code{\link{H2fun}}}{to compute H2' for non-integer web. Thanks to Jochen Fründ for providing the inspiration and the code. This option is also now available in \code{\link{networklevel}}.}
    \item{Changes to \code{\link{nullmodel}}:}{Upon using the modification entered in version 1.12, it turns out this is no improvement but rather a step in the wrong direction. So this function is back to its pre-1.12 working with a little bug fixed in swap.web (the helper function `upswap' got an option too many) that caused the alteration in the 1.12 version. Thanks to Jochen Fründ for insisting on rolling this back.}
    \item{Bug fix in \code{\link{networklevel}}:}{Alatalo interaction evenness could not be called separately due to a change of name within the function. Thanks to Rachel Gibson for reporting!}
   \item{New index added to \code{\link{networklevel}}:}{In addition to the binary cluster coefficient, the function now also (and automatically) computes the weighted clustering coefficient introduced by Tore Opsahl in his package \pkg{tnet}.}
    \item{New function \code{\link{web2edges}}:}{Transforms a web-matrix into an edge list, as used in \pkg{tnet} or other software (e.g. Pajek).}
 }

  \item 1.12 (release date: 21-May-2010)
  \describe{   
    \item{Bug fix in \code{\link{vaznull}}:}{In very poorly sampled networks, \code{\link{vaznull}} could fail because after the initial filling there were no more interactions to distribute. How an if-statement can sometimes make a difference.}
    \item{Error message in \code{\link{nullmodel}} changed to a warning.}{Up to this release, the function created null models as specified. However, when a supposedly quantitative network was in fact binary, it returned an error. Now it proceeds, using the null model generating algorithm ``mgen'', and returns a warning.}
}	
  
  \item 1.11 (release date: 10-May-2010)
  \describe{   
    \item{Changes to function \code{\link{as.one.mode}}:}{Gains options to do projections of the bipartite (=two-mode) into the one-mode modus required e.g. by \pkg{sna}. In most publications, one focusses on only one trophic level (e.g. the pollinators) and represents the bipartite network as a one-mode, pollinators-only network. This can now be achieved using the various options in the function. See its help for details.}
    \item{Bug fixes in \code{\link{BC}} and \code{\link{CC}}:}{Upon re-reading the paper motivating the inclusion of these functions, I noticed an error in the previous implementation. The key point, and one that is contentious in network theory, is that there is no standard best way to project a two-mode (bipartite) network onto a one-mode network. Up until know, I used the default of \code{\link{as.one.mode}}. However, Martín Gonzáles et al. (2010) use the more common projection (``higher'' and ``lower'', respectively), and this is also now implemented in \code{\link{BC}} and \code{\link{CC}}. The interested reader may want to follow this problem up by reading the pre-print of Tore Opsahl (``Triadic closure in two-mode networks'') to see that this one-mode projection actually distorts the statistics employed on them.}
}
  
  \item 1.10 (release date: 02-May-2010)
  \describe{
    \item{New function \code{\link{null.distr}}:}{Fits a lognormal or negative binomial distribution to the marginal totals of a network and then draws random values from the thus specified distribution. These values represent the marginal totals to be used for constructing a null model network. The function thus serves as a way to produce null models without maintaining exactly the marginal totals, but only their distribution. Please refer to the help for some cautionary remarks on its use! In particular, null model networks are likely to be smaller than the original and there is obviously no guarantee that the marginal totals are actually distributed in the specified way!}
    \item{Polishing the Vázquez pollination networks:}{All species in these 8 networks are now consistently named (rather than numbered).}
    \item{Polish to \code{\link{specieslevel}}:}{Requesting Fisher's alpha sometimes causes convergence problems in the underlying function. This is now caught by a try, returning NA.}
    \item{Rescale in \code{\link{nested}} gave wrong scaling for NODF:}{Fixed.}
  }

  \item 1.09 (release date: 22-Mar-2010)
  \describe{
    \item{Bug fix in \code{\link{degreedistr}}:}{Degree distributions for the higher trophic level were calculated incorrectly. The main problem was that the scaling constant was omitted, leading to ridiculous fits. (Thanks to Roberto Molowny for reporting, correction proposal and discussion!) On improving this, I also ventured to provide multiple starting values and hence the function returns a fit for all three curves more often (but not always).}
  }
  
  \item 1.08 (release date: 16-Mar-2010)
  \describe{
    \item{Bug fix in \code{\link{ND}}:}{Rows and columns were mixed up and a +1 was missing- oh dear. (Thanks to Marco Mello for reporting!)}
  }
  
  \item 1.07 (release date: 18-Feb-2010)
  \describe{
    \item{New function \code{\link{vaznull}}:}{Null model with constrained totals and connectance proposed by Diego Vázquez, similar to \code{\link{swap.web}}, but ``better''. See notes in \code{\link{swap.web}}-help for justification of this statement. We recommend this null model for constraining both marginal totals and connectance.}
    \item{Added a new option to calculate extinction slopes:}{So far, random extinction sequences and the rarest-to-commenest sequence were implemented. Now, \code{\link{second.extinct}} also includes the option \option{method="degree"} to build a sequence from the best-to-least connected species. This is the most extreme case, where the most generalist species goes extinct first (see Memmott et al. 1998). (In response to request by Simone Bazarian.)}
    \item{Adaptations of \code{\link{nullmodel}}:}{Now includes the new vaznull as a method (3). This leads to changes in the sequence of methods! Old code 
may hence have to be adapted accordingly!}
  }
  \item 1.06 (release date: 18-Dec-2009)
  \describe{
    \item{New functions \code{\link{ND}}, \code{\link{BC}} and \code{\link{CC}}:}{Simple functions to calculate normalised degree, betweenness centrality and closeness centrality. These functions and the example allows a reproduction of the type of analysis carried out in Martín González et al. (2009).}
    \item{Changes to \code{\link{specieslevel}}}{, which now calls \code{\link{ND}}, \code{\link{BC}} and \code{\link{CC}}, too.}
    \item{Bug fix in \code{\link{as.one.mode}}:}{now allows also data.frames to be turned into one-mode-style representations. Previously, only matrices could be used.}
    \item{Suppression of errors and warnings in \code{\link{networklevel}}:}{When used on a full network (i.e. one without zeros), some indices in \code{\link{networklevel}} are undefined (e.g. extinction slopes). This led to a long output of warning and error messages, although internally I used the try-function to capture errors. Now, these messages are suppressed.}
  }
  
  \item 1.05 (release date: 05-Dec-2009)
    \describe{
      \item{Help description for \code{\link{PAC}} was convoluted}{ and its suggestions for the interpretation downright wrong. Thanks to Matthew Wainhouse for reporting and simplifying it!}
      \item{Resolved a long-standing issue of a warning message.}{This occurred when detaching the package (\code{detach(package:bipartite)}) and was caused by somehow wrongly using \code{.Last.lib}. Deleting it solved the problem.}
      \item{Fixed bug in \code{\link{specieslevel}}:}{Calling this function with only one index caused it to return empty lists. Simple mistake, but better without it.}		
    }
  
  \item 1.03 (release date: 06-Nov-2009)
    \describe{
      \item{Bug fix in \code{\link{V.ratio}}:}{A small mistake for a human, but a huge bug for mankind. Sorry. (Detail: I forgot to square sigma.i, leading to strange results.) Thanks to Giorgio Mancinelli for reporting!}
      \item{Error message for non-existent indices in \code{\link{networklevel}}:}{So far, \code{\link{networklevel}} returned \code{NULL} when an index was selected that does not exist (e.g. \code{"shannon diversity"} instead of \code{"diversity"}). Now, a helpful (?) error message is returned.}
      \item{Change in defaults to index \code{"interaction evenness"} in \code{\link{networklevel}}:}{After a fruitful discussion with Becky Morris and Jason Tylianakis, I changed the default to \code{"sum"}. I also reproduce some of our communication in the help to this function, under details, to make the ecological assumptions behind either option a bit more transparent. There are good reasons for either option.}
    }
  \item 1.02 (release date: 11-Sep-2009)
    \describe{
      \item{Function \code{\link{plotweb}}:}{Now more trophic networks can be plotted by staggering bipartite networks on top of each other. See multitrophic examples in \code{\link{plotweb}}.} 
      \item{Minor corrections to \code{\link{dfun}}:}{This function did not return \emph{exactly} the values of the website-version. Jochen Fründ corrected this. Please read the help of the function (final paragraph) for details.}
      \item{New function \code{\link{nested}}:}{This convenience function collects the various ways to calculate nestedness of a network in order to facilitate comparison of nestedness analyses. To do so, it heavily borrows from \pkg{vegan}.}
      \item{Bug fix to \code{\link{discrepancy}}:}{Would return a silly value (half of the number of rows) for empty matrices. This had no effect when called by \code{\link{networklevel}}, since the matrix would have been emptied. Thanks to Roberto Molowny for reporting!}
      \item{Bug fix to \code{\link{networklevel}}:}{Due to a missing space, the option \option{ALLBUTDD} did not work properly. Thanks to Etienne Laliberte for reporting!}
    }
  \item 1.00 (release date: 06-Aug-2009)
    \describe{
      \item{Complete overhaul of \code{\link{networklevel}}:}{After a workshop on bipartite networks in ecology, a few more indices were added (Fisher's alpha diversity of interactions, mean interaction diversity, mean number of predators) and the whole output reorganised. It now follows a gradient from less to more interesting (in our view) indices, and from indices for binary to those for quantitative networks. Also, I added several options for which indices to report (\option{index="info", "binary", "quantitative", "topology"}). Most interestingly, perhaps, there is now a quantitative, Shannon-based series of indices. Starting with the \dQuote{mean interaction diversity} (i.e. the Shannon-diversity of interactions of a species, averaged across all species in that trophic level), over \dQuote{Shannon diversity} of interactions, to H2 (i.e. Shannon diversity scaled between max and min possible for this web characteristics).}
      \item{Additions to, and overhaul of \code{\link{specieslevel}}:}{Similar to the above item, some indices were added, and the output simplified when calling the option \option{index="ALLBUTD"} (only one D here!): a list with two matrices is now returned. Fisher's alpha for each species, vulnerability/generality and effective number of species for each target species are also now included. Index sequence has changed.}
      \item{New function \code{\link{PAC}}:}{Calculates the Potential for Apparent Competition following the formula in Müller et al. (1999) and Morris et al. (2005) and the suggestion by Becky Morris and Owen Lewis. More than a theoretical concept, it was experimentally shown to be relevant (Morris et al. 2004).}
      \item{Bug fix in \code{\link{H2fun}} and \code{\link{dfun}}:}{A line of code went missing at some point, so the maximum packing density was not optimal (but still good) in these functions. As a result, reported H2- and d-values were sometimes 0 when they should be only very close to 0.}
      \item{Switch for error reporting in \code{\link{degreedistr}}:}{By default now suppresses error reporting when \code{nls} fails to fit a degree distribution due to too few data points. This leaves the user of \code{\link{networklevel}} somewhat less confused. Warning message now also indicates for which trophic level there were too few data points.}
      \item{Bug fix in \code{\link{C.score}}:}{Failed when the web was very dense or very sparse, because the maximum number of checkerboard patterns was 0 then.}
    }
  \item 0.94 (release date: 01-Aug-2009)
    \describe{
      \item{New function \code{\link{nullmodel}}:}{A convenience wrapper function to generate different types of null models.}
      \item{Small changes to \code{\link{networklevel}}:}{This function returns a list of indices. If we exclude the computation of degree distribution fits, this would be coerced to a vector. We added the option \option{index="ALLBUTDD"} to calculate all indices BUT degree distributions. The output is then returned as vector. This is much more convenient when using networklevel on many data sets (using \code{sapply}).}
    }

  \item 0.93 (release date: 30-Jun-2009)
    \describe{
      \item{bug fix in \code{\link{slope.bipartite}} and \code{\link{robustness}}:}{The function always selected column 3 of the object, instead of 2 for lower and 3 for higher trophic level; thanks to Antonio Rivera for spotting and reporting this error! Notice that this error must have slipped in somewhere around version 0.90/0.91, because I checked and the results reported in the Open Ecology Journal paper are valid! I seem to remember that I (CFD) modified \code{\link{slope.bipartite}} when \code{\link{robustness}} was added: never change a wining team!}
    }
  \item 0.92 (release date: 02-Jun-2009)
    \describe{
      \item{more colour options in \code{\link{visweb}}:}{the arguments \option{box.border} and \option{box.col} now allow a specification of the colour of the boxes and their borders.}
      \item{bug fix \code{\link{empty}}:}{the function returned 0 for a 1x1 matrix. Although not written for such a case, it should still do what it says on the 
tin. Thanks to Mariano Devoto for spotting and reporting!}
      \item{bug fix \code{\link{wine}}:}{returned NA for square matrices.}
      \item{bug fix \code{\link{plot.wine}}:}{gave decimal places for row and column names for very small networks.}
    }

  \item 0.91 (release date: 06-May-2009)
  \describe{
    \item{new function \code{\link{wine}}:}{This function replaces the (now deprecated) function \code{nestedness.corso} in calculating a (weighted) nestedness for bipartite networks. It was developed and implemented by Marcelino de la Cruz, Juan M. Pastor, Javier Galeano and Jose Iriondo. It is also called by \code{\link{networklevel}}. A plotting function is also available, depicting the contributions of each observed link to the web's nestedness. - \code{nestedness.corso} is now removed from the package. The main reason is that it served as an interim solution for \code{\link{wine}}, and the Corso-way of calculating nestedness is just one more of already too many. If you intend to use \code{\link{wine}} on binary data and interpret that as a Corso-equivalent, beware of the following two main differences (thanks to Jose Iriondo for summarising them): \emph{First, the nestedness of Corso et al. varies between 0 and 1, with the highest nestedness is reached at 0 and 1 corresponds to random, whereas in \code{\link{wine}} is just the opposite (this is because the Manhattan distances are calculated with regard to opposite sides of the matrix). Secondly, the \sQuote{d} in the nestedness of Corso is the sums of the \sQuote{dij}s whereas in \code{\link{wine}}, the \sQuote{d} (= \sQuote{win} value of the object produced by \code{\link{wine}}) is the average of the \sQuote{dij}s above 0.} So, we recommend NOT to use wine for calculating Corso's nestedness, but rather download the source code for \code{nestedness.corso} from an older version of \pkg{bipartite}.}
    \item{new function \code{\link{robustness}}:}{A better way to quantify the effect of species loss on the extinctions in the other trophic level; kindly provided by Mariano Devoto. This index is also part of \code{\link{networklevel}}.}
    \item{new data set \code{ollerton2003}:}{Another quantitative pollination network from the NCEAS database (see \code{\link{ollerton2003}}).}
  }

  \item 0.90 (release date: 24-Mar-2009)
  \describe{
    \item{example \code{vazquez.example}:}{We added several new functionalities mainly to be able to analysis data and use network statistics as suggested by 
Vázquez et al. 2009. You can access this example by typing \code{?vazquez.example}.}
    \item{new feature in \code{\link{visweb}}:}{Can now plot different sized circles to represent interactions, as proposed by Vázquez et al. (2009).}
    \item{new feature in \code{\link{networklevel}}:}{We added an option to calculate interaction evenness either based on all possible links or just on realized links.}
    \item{new function \code{\link{sortweb}}:}{Can be used to sort webs in different ways.}
    \item{data set \code{inouye1988}:}{Another pollination network from the NCEAS database (see \code{\link{inouye1988}}).}
    \item{Function \code{\link{compart}}:}{We replaced the CA-based approach to detecting compartments by a comprehensive and recursive approach. The latter is not affected by ties in the data set (i.e. species with the same number of links). In quantitative webs and for the networks included in bipartite, the old function was working fine, but in more recent trials it failed to detect 2-species compartments. In turn, we had to adapt \code{\link{networklevel}} and \code{\link{plotweb}}.}
  }

\item 0.85 (release date: 10-Mar-2009)
  \describe{
	\item{Function \code{\link{plotweb}}}{New feature: text labels can now be printed in different colours. All colours can be passed as vectors and vectors are recycled if not of appropriate length }
	}
\item 0.84 (release date: 25-Feb-2009)
  \describe{
	\item{Function \code{\link{plotweb}}}{New feature: arrow=\dQuote{center.up}, arrow=\dQuote{center.down}, arrow=\dQuote{center.both}: this results in the standard presentation of bipartite networks, where interactions are displayed as centered triangles. }
	}

\item 0.83 (release date: 9-Feb-2009)
	pdf of OEJ-paper added as vignette

\item 0.82 (release date: 27-Jan-2009)
\describe{
	\item{Function \code{\link{plotweb}}}{Two new features: 1. colors of borders of boxes and interactions can be specified via setting bor.col. 2. 
Labels can now be rotated by specifying text.rot=90.}
}

\item 0.81 (release date: 06-Jan-2009)
\describe{
	\item{Function \code{\link{nodespec}} replaces the now obsolete function \code{functspec}.}{The inventor of functional specialisation, 
Bo Dalsgaard, understands the term \sQuote{functional specialisation} to be restricted to the way in which plants are specialised to pollinators. Although I disagree (and think that \sQuote{functional} actually means very little until defined for a given problem), I renamed the \sQuote{functional specialisation index} into \sQuote{node specialisation index}. This is not a particularly clever name, but at least it indicates that the position of nodes in a network is important when it is calculated.}
	\item{Bug in \code{\link{networklevel}}:}{Call to \code{nestedness.corso} overwrote results of \code{\link{nestedness}}.}
	\item{Various changes related to the analysis of very small webs.}{Very small webs should probably not be used at all for testing theories! 
Indices are usually VERY sensitive to the exact number of species, number of observations etc. Still, sometimes we simply want to calculate some index, and then \pkg{bipartite} should handle such small networks, too. So, when going through several dozens of very small webs (sometimes only containing one species in one of the two trophic levels), several functions did not perform correctly (usually to such minor programming issues such as matrices being converted to vectors when \code{[} was used with \option{drop=TRUE}, i.e. the default). Changes affected the functions \code{discrepancy}, \code{nestedness.corso}, \code{shuffle.web}, \code{compart} and \code{empty}. Their output remains identical, only they now also work for small webs.}
	\item{Changes to \code{nestedness.corso} with \option{weighted=TRUE}}{Galeano et al. do not describe how to deal with ties, nor do they make clear if the packed matrix should be sorted by marginal sum of links or marginal sum of interactions! Previously, we used marginal sum of interactions (because it is a weighted index), but now we moved to marginal sum of links, because that is how I interpret their paper after a third (or forth) re-reading. Furthermore, because the real maximum chaos cannot be derived (to my knowledge) algorithmically, we use the 95\% quantile of 500 randomisations as maximum. This will lead to a consistently overestimated nestedness, but it is less sensitive to the number of replicates than the max. (Also, there was an error in the description of the value returned: 0 is nested, 1 is chaos!)}
}

\item 0.8 (release date: 21-Dec-2008)
\describe{
	\item{bug fix in \code{\link{shuffle.web}}}{Function didn't work correctly when the web contained more rows than columns. (Thanks to Anna Traveset for spotting and reporting!)}
	\item{minor modifications in networklevel}{Sometimes additional information is available and species can be included in a network, although they have no interactions with other species. In this case, one might want to use these species in the network, too. A new option (empty.web), allows to keep empty rows and columns, although for some functions an emptied web must be used (e.g. degree distributions).}
	\item{bug fix in \code{\link{H2fun}}}{As in the last \code{\link{H2fun}} bug fix, sometimes H2 became negative.}
	\item{bug (?) fix in \code{\link{networklevel}}}{Shannon diversity is based on the log of interactions. If this value is 0, as it is for most network entries, an NA is produced. As a consequence, Shannon's H (now also given in output) is based not on all interactions, but only those > 0. In consequence, Shannon evenness should also only be SH divided by the number of realised interactions (\code{log(sum(web>0))}). That is now the case.}
}

\item 0.74 (release date: 24-Oct-2008)
\describe{
	\item{functional specialisation (\code{functspec}) bug fix}{Paths were double the true length, hence minimum was 2, rather than 1.}
	\item{\code{\link{H2fun}} bug fix}{Since the search for H2min is heuristic, H2uncorr can sometimes be lower than H2min; in that case, H2fun returned a value greater 1, while it should be one exactly.}
	\item{new function \code{nestedness.corso} }{Calculates (weighted) nestedness according to Corso et al. (2008) and Galeano et al. (2008).}
	\item{new function \code{\link{discrepancy}}}{Calculates discrepancy according to Brualdi & Sanderson (1999), deemed to be best ever measure of nestedness; also gives an example for a binary null model analysis based on \pkg{vegan}'s \code{commsimulator}.}
	\item{other}{Correction of several minor typos on the help pages; removal of "~" in help files; same citation style throughout; new cross references (especially for the nestedness functions); in \code{\link{networklevel}}, nestedness is now calculated using \pkg{vegan}'s \code{nestedtemp} due to matrix inversion problems reported for binmatnest.}
}

\item 0.73 (release date: 1-Sept-2008)
\describe{	
	\item{new feature}{\code{\link{plotweb}}}{Named abundance-vector for each level can be used.}
	\item{new function \code{\link{plotweb2}} (not debugged!)}{For plotting tripartite networks.}
}

\item 0.72 (release date: 12-June-2008)
  \describe{
	\item{new function: functional specialisation with \code{functspec}}{See Dalgaard et al. (2008).}
	\item{new function: interface to \pkg{sna} through \code{\link{as.one.mode}}}{Allows calculation of path lengths, centrality, betweenness and other indices developed for one-mode networks.}
	\item{bug-fix}{Error in \code{\link{plotweb}} when no species labels were given.}
  }
 } % end itemize within versionlog
} % end section versionlog
                        
\note{This function is only invoked for its side effect of opening the help page. I simply didn't know how to do it any other way ...}   
    
\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }
    
\keyword{ package}
    
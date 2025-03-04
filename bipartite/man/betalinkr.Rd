\encoding{UTF-8}
\name{betalinkr}
\alias{betalinkr}
\alias{betalinkr_multi}

\title{
Calculate network dissimilarity (beta diversity) and its components
}

\description{
This function (\code{betalinkr}) is a new implementation of network dissimarility, as proposed by Tim Poisot (originally implemented in the \emph{betalink} package). Following Poisot, dissimilarity (of a pair of networks) is partitioned into the dissimilarity due to difference in species composition (``ST'') and dissimilarity due to rewiring (``OS'', dissimilarity of shared species subweb). Different partitioning approaches and many binary and quantitative indices are available. The recommended method for additive partitioning is \code{commondenom} (originally proposed by Novotny 2009), for which further partitioning into different aspects of species composition differences (\code{partition.st}) and into replacement and richness difference components (\code{partition.rr}) are also available.

\code{betalinkr_multi} is a metafunction that calls \code{betalinkr} for all pairs of networks (passing arguments), and returns a data.frame.
}

\usage{
betalinkr(webarray, index="bray", binary=TRUE, partitioning="commondenom",
          proportions=!binary, function.dist="vegdist", distofempty="zero",
          partition.st=FALSE, partition.rr=FALSE)

betalinkr_multi(webarray, ...)
}
\arguments{
  \item{webarray}{
    Input data, an array with three dimensions, third dimension has length 2 and separates the two webs to compare. First two dimensions are the species (typically first Lower and second Higher). For convenience, if \option{webarray} is a list of two webs, it will be converted to array format first using \code{webs2array} (that function can also be used to convert single web matrices to array format, but then needs to be called explicitly). NOTE: When using \code{betalinkr_multi}, \option{webarray} \bold{must} be an array! Use \code{\link{webs2array}} when providing networks (see last example below).
  }
  \item{index}{ 
    The dissimilarity index, passed to "method" of either vegdist or betadiver (see there for naming). If \option{partitioning="commondenom"} (and thus no external function is used), it has to be one of \option{sorensen} (for Sorensen or Bray-Curtis dissimilarity index) or \option{jaccard} (for Jaccard or Ruzicka); \option{bray} is also allowed and gives the same result as \option{sorensen}.  }
  \item{binary}{
  Should binary data or quantitative data be used (i.e., quantitative or binary versions of dissimilarity indices). If \code{TRUE}, webs will be converted to binary.
  }
  \item{partitioning}{
    How should the components (``ST'' and ``OS'') be calculated? With \option{poisot}, the original approach will be applied, which calculates ST as WN - OS; not recommended for partitioning as it underestimates the contribution of species composition. But with this setting, OS will reflect the uncorrected dissimilarity of shared species subwebs, which may be interesting as such.

    With \option{commondenom} (recommended setting), an alternative partitioning approach is used that should give fair contributions for OS and ST. It keeps a common denominator (of the dissimilarity index) for WN, OS and ST (i.e., the denominator is based on all links of the pair of webs instead of on subsets). Thus, only the numerator of the dissimilarity index is partitioned, which ensures additivity directly by partitioning the sets of links (or link weights).
  }
  \item{proportions}{
    Should data be standardized to proportions before calculating quantitative dissimilarity metrics? Note that this is done here directly before calculating dissimilarities, thus gives a different result than if input data were already standardized. Will be set to TRUE by default if \option{binary=FALSE} and to FALSE otherwise.
  }
  \item{function.dist}{
    Which function to use for calculating dissimilarity? The name of \bold{vegan} package function, either \option{vegdist} (quantitative indices available) or \option{betadiver} (gives compatibility to the 24 numbered indices of Koleff et al. (2003) as used by \pkg{betalink} package, but only binary indices). Note that this argument is ignored if \code{partitioning="commondenom"}, which implements the dissimilarity index calculation directly.
  }
  \item{distofempty}{
    Can be \option{zero} or \option{na}. How should dissimilarity be defined when there are either no links to use for OS (i.e. only links involving species just found in one of the 2 webs) or no links to use for ST (i.e. only ``rewiring links'' present)? \option{zero} is appropriate when interested in the contribution of components OS and ST, whereas \option{na} is appropriate when b_os should be interpreted separately as dissimilarity of shared species subwebs
  }
  \item{partition.st}{
    A secondary partitioning of the species turnover component, following Novotny (2009), which separates ST further into dissimilarity due to absence of resource species (lower, ST.l), absence of consumer species (higher, ST.h) or absence of both (ST.lh)
        
    Only works with \option{partitioning="commondenom"}.
  }
  \item{partition.rr}{
    A secondary partitioning of dissimilarity (WN and OS) into ``true'' turnover components (WN.repl, OS.repl) and a component due to richness difference (difference in network totals; WN.rich, OS.repl). This follows Legendre (2014, and references therein) and was first applied to networks by Noreika et al. (2019). The function implements ``Podani family'' indices according to Legendre's terminology.
    
    Only works with \code{partitioning="commondenom"}.
  }
  \item{...}{Arguments passed on to \code{\link{betalinkr}}.}
}
\details{
The basic idea to calculate dissimilarity (betadiversity) between networks (links instead of species) has been proposed before. Poisot et al. (2012) came up with the idea to separate between rewiring and species turnover as causes of network dissimilarity (but see Novotny 2009). They proposed to calculate rewiring link dissimilarity (``beta_OS'' or simply ``OS'') by focusing on the sub-web containing only species observed in both webs (i.e. excluding links with species unique to one of the webs). Species turnover link dissimilarity (``beta_ST'' or simply ``ST'') is then calculated as total dissimilarity minus rewiring dissimilarity, but this assumes an additivity that is rarely given with dissimilarity indices.
Although dissimilarity values may not be additive, the number of ``unique links'' (i.e. only observed in one of the two webs) can be well partitioned into additive components.
% Here, I propose an alternative approach to partition the two abovementioned components of network dissimilarity.

With option \option{partitioning="poisot"}, many different dissimilarity indices can be used with \code{function.dist} either \code{vegdist} or \code{betadiver}, but no guarantee that all of them can be usefully interpreted (following Legendre 2014, Jaccard-family and Sorensen-family dissimilarity indices are recommended; note that with vegdist you get Sorensen-family with \option{bray}, and with betadiver Jaccard is 15 and Sorensen is 1, although these go under different names there as Koleff et al. (2003) use the names Jaccard/Sorensen for the corresponding similarity metrics).

The alternative approach (with \option{partitioning="commondenom"}) was inspired by Legendre (2014). It avoids using an existing dissimilarity function, but rather implements calculation of dissimilarity indices directly (thus more limited options for these: only binary and quantitative Jaccard-type and Sorensen-type). Here, the same denominator as for the total dissimilarity WN is used also for its components OS and ST, thus ensuring additivity. This method was actually already proposed by Novotny (2009), who further splitted ST (see \option{partition.st}). If your goal is indeed to partition dissimilarity into its (additive) components, use this method!

Note that if you are interested only in \emph{dissimilarity between subwebs (not as an additive component)}, you should use \option{partitioning="poisot"} and look at OS while ignoring ST. This works also for quantitative data (with \code{binary=FALSE}). If this is your goal, you should probably set \option{distofempty="na"}, returning NA/NaN if there are no interactions between shared species.

To generate results identical to Poisot's \code{betalink} function in package \emph{betalink}, use these settings: \option{partitioning="poisot", function.dist="betadiver", distofempty="na" and binary=TRUE}.

Again: the output for OS can differ strongly depending on the chosen options! You have to decide: do you want dissimilarity, which is inherently not additive (use \option{partitioning="poisot"}), or do you want an additive dissimilarity component, which is not a dissimilarity itself (use \option{partitioning="commondenom"})?
}
\value{
A named vector of four (or more) dissimilarities (components), naming follows Poisot et al. 2012.
\item{S}{beta_S, the dissimilarity in species composition}
\item{OS}{beta_OS, the dissimilarity (component) explained by ``rewiring'' among shared species}
\item{WN}{beta_WN, the dissimilarity between the two networks}
\item{ST}{beta_ST, the dissimilarity (component) explained by difference in species community composition}
\item{(others)}{possibly more elements, same names as above, but with secondary partitioning added after ``.'', see \code{partition.st} and \code{partition.rr}}
For \code{betalinkr_multi}, output is a dataframe with one row per pair of webs compared.
}

\references{
Poisot, T., E. Canard, D. Mouillot, N. Mouquet, D. Gravel, and F. Jordan. 2012. The dissimilarity of species interaction networks. \emph{Ecology Letters} \bold{15}, 1353--1361. %\doi{10.1111/ele.12002}

Legendre, P. 2014. Interpreting the replacement and richness difference components of beta diversity. \emph{Global Ecology and Biogeography} \bold{23}, 1324--1334. %\doi{10.1111/geb.12207}

Koleff, P., Gaston, K.J., and J.J. Lennon. 2003. Measuring beta diversity for presence–absence data. \emph{Journal of Animal Ecology} \bold{72}, 367--382.

Novotny, V. 2009. Beta diversity of plant-insect food webs in tropical forests: a conceptual framework. \emph{Insect Conservation and Diversity} \bold{2}, 5--9. %\doi{10.1111/j.1752-4598.2008.00035.x}

Noreika, N., Bartomeus, I., Winsa, M., Bommarco, R., and  E. Öckinger. 2019. Pollinator foraging flexibility mediates rapid plant-pollinator network restoration in semi-natural grasslands.\emph{Scientific Reports}, \bold{9}, 1--11. \doi{10.1038/s41598-019-51912-4}
}

\author{Jochen Fründ}

\note{
  This function allows to use quantitative dissimilarity indices, which are usually recommended. However, for quantitative networks it is far from trivial how to correctly separate (or even define) which part of the dissimilarity is due to rewiring and which due to difference in species composition! Here I use the concept that all variation between the subwebs of shared species can be attributed to rewiring, BUT this will most likely not be correct. Even if all species are shared among two networks, quantitative species dissimilarity may be large (different [relative] abundances), and this will most likely lead to changes in network frequencies (changing link weights, or even missing links) that should not be called rewiring. How to correctly define and measure that is open to discussion, but thus I would still consider the values for beta_OS as overestimates.
  
  This function should also work for one-mode networks.

  Why the name of the function? Short for ``betalink revised''.
  
  Why the names of the output values? These are the indices of beta (for betadiversity); I am keeping the names used by Tim Poisot and guess what they stand for: S (Species), OS (Only Shared species links), WN (Whole Network links), ST (Species Turnover links). Some authors would call all this dissimilarity (or dissimilarity components) and reserve the term betadiversity for something not standardized between 0 and 1.
  
  For partitioning into replacement and richness difference, note that replacement in WN can mean richness difference in OS (if shared species switch from interaction with other shared species in one network to interactions with non-shared species in another network, changing the size of shared-species subweb), so finding OS.rich > WN.rich is not a bug.
  
  Thanks to Carsten Dormann, Benjamin Schwarz, Benoit Gauzens, Nacho Bartomeus and Timothee Poisot for their contributions to this function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[vegan]{vegdist}} and \code{\link[vegan]{betadiver}} for vegan package functions calculating dissimilarity / betadiversity, used here unless \option{partitioning="commondenom"}.
  
  For reshaping web data to the array input format expected here, see \code{\link{webs2array}} and \code{\link{frame2webs}}.
}

\examples{
# two examples that give the same results as would the 
# \code{betalink} function in the package of the same name
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), 
    partitioning="poisot")
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), 
    function.dist="betadiver",index=1, partitioning="poisot")

# same data, with recommended partitioning method plus further partitioning
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)), 
    partitioning="commondenom", partition.st=TRUE)

# another example (no shared links)
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
  lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
  webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))

# more than two webs:
betalinkr_multi(webs2array(Safariland, vazquenc, vazarr), index="jaccard")
}

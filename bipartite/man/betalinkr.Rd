\encoding{UTF-8}
\name{betalinkr}
\alias{betalinkr}
\alias{betalinkr_multi}

\title{
Calculate network dissimilarity (beta diversity) and its components
}

\description{
This function (\code{betalinkr}) is a new implementation of network dissimarility, as proposed by Tim Poisot (originally implemented in the \emph{betalink} package). Following Poisot, dissimilarity (of a pair of networks) is partitioned into the dissimilarity due to difference in species composition ("ST") and dissimilarity due to rewiring ("OS", dissimilarity of shared species subweb). Different partitioning approaches and many binary and quantitative indices are available.

\code{betalinkr_multi} is a metafunction that calls \code{betalinkr} for all pairs of networks (passing arguments), and returns a data.frame.
}

\usage{
betalinkr(webarray, index="bray", binary=TRUE, partitioning="poisot", 
  proportions=FALSE, function.dist="vegdist", distofempty="na")

betalinkr_multi(webarray, ...)
}
% the defaults are currently set to reproduce Poisots results; should rather change to my preferred settings (but then there should be an argument "set.defaults" that can be "poisot")
\arguments{
  \item{webarray}{
    Input data, an array with three dimensions, third dimension has length 2 and separates the two webs to compare. First two dimensions are the species (typically first Lower and second Higher). For convenience, if webarray is a list of two webs, it will be converted to array format first using \code{webs2array} (that function can also be used to convert single web matrices to array format, but then needs to be called explicitly).
  }
  \item{index}{ 
    The dissimilarity index, passed to "method" of either vegdist or betadiver (see there for naming). If \code{partitioning="commondenom"} (and thus no external function is used), it has to be one of "sorensen" (for Sorensen or Bray-Curtis dissimilarity index) or "jaccard" (for Jaccard or Ruzicka).  }
  \item{binary}{
  Should binary data or quantitative data be used (i.e., quantitative or binary versions of dissimilarity indices). If \code{TRUE}, webs will be converted to binary.
  }
  \item{partitioning}{
    How should the components ("ST" and "OS") be calculated? With \code{"poisot"}, the original approach will be applied, which calculates ST as WN - OS; not recommended as this underestimates the contribution of differences in species composition (but with this setting, OS will reflect the uncorrected the dissimilarity of shared species subwebs). Two alternative approaches that give fair contributions for OS and ST are implemented:
    
    With \code{"adjusted"}, ST is calculated equivalently to the calculation of OS, and then OS and ST are adjusted so that they sum to WN. See details below.
    
    With \code{"commondenom"}, an alternative partitioning approach is used, which keeps a common denominator (of the dissimilarity index) for WN, OS and ST (i.e., the denominator is based on all links of the pair of webs instead of on subsets). Thus, only the numerator of the dissimilarity index is partitioned, which ensures additivity directly by partitioning the sets of links (or link weights).
  }
  \item{proportions}{
    Should data be standardized to proportions before calculating quantitative dissimilarity metrics? Note that this is done here directly before calculating dissimilarities, thus gives a different result than if input data are already standardized.
  }
  \item{function.dist}{
    Which function to use for calculating dissimilarity? Either the name of \code{vegan} package function, either \code{"vegdist"} (quantitative indices available) or \code{"betadiver"} (gives compatibility to the 24 numbered indices of Koleff et al. 2003 as used by \emph{betalink} package, but only binary indices). Note that this argument is ignored if \code{partitioning="commondenom"}, which implements the dissimilarity index calculation directly.
  }
  \item{distofempty}{
    Can be \code{"zero"} or \code{"na"}. How should dissimilarity be defined when there are either no links to use for b_os (i.e. only links involving species just found in one of the 2 webs) or no links to use for b_st (i.e. only "rewiring links" present)? "zero" is appropriate when interested in the contribution of components b_os and b_st, whereas "na" is appropriate when b_os should be interpreted separately as dissimilarity of shared species subwebs
  }
  \item{...}{Arguments passed on to \code{\link{betalinkr}}.}
  % another argument "2nd.decomposition" that further calculates B_richness and B_replacement components of each dissimilarity will come; maybe only available with partitioning="commondenom"
}
\details{
The basic idea to calculate dissimilarity (betadiversity) between networks (links instead of species) has been proposed before. Poisot et al. (2012) came up with the idea to separate between rewiring and species turnover as reasons for network dissimilarity. They proposed to calculate rewiring link dissimilarity ("beta_OS" or simply "OS") by focusing on the sub-web containing only species observed in both webs (i.e. excluding links with species unique to one of the webs). Species turnover link dissimilarity ("beta_ST" or simply "ST") is then calculated as total dissimilarity minus rewiring dissimilarity, but this assumes an additivity that is rarely given with dissimilarity indices.

However, the number of "unique links" (i.e. only observed in one of the two webs) can be well partitioned into additive components. Here, I propose two alternative approaches to partition the two abovementioned components of network dissimilarity.

In the first approach (with \code{partitioning="adjusted"}, beta_OS is first calculated following Poisot et al., and beta_ST is then calculated equivalently, but excluding rewiring links instead of links of unique species from the comparison. In a second step, these raw values are adjusted (multiplying by WN/(OS+ST)) in order to fulfill OS+ST=WN. With this partitioning approach, many different dissimilarity indices can be used with \code{function.dist} either \code{vegdist} or \code{betadiver}, but no guarantee that they all of them can be usefully interpreted (following Legendre 2014, Jaccard-family and Sorensen-family dissimilarity indices are recommended; note that with vegdist you get Sorensen-family with "bray", and with betadiver Jaccard is 15 and Sorensen is 1, although these go under different names there as Koleff et al. 2003 uses the names Jaccard/Sorensen for the corresponding similarity metrics).

The second approach (with \code{partitioning="commondenom"}) was inspired by Legendre 2014 and avoids using an existing dissimilarity function, but implements calculation of dissimilarity indices directly (thus more limited options for these). Here, the same denominator as for the total dissimilarity WN is used also for its components OS and ST, thus ensuring additivity.
Note that if you are interested only in dissimilarity between subwebs (not as an additive component), you should use \code{partitioning="poisot"}.
}
\value{
A named vector of four dissimilarities, naming follows Poisot et al. 2012.
\item{S}{beta_S, the dissimilarity in species composition}
\item{OS}{beta_OS, the dissimilarity (component) explained by "rewiring" among shared species}
\item{WN}{beta_WN, the dissimilarity between the two networks}
\item{ST}{beta_ST, the dissimilarity (component) explained by difference in species community composition}
}

\references{
Poisot, T., E. Canard, D. Mouillot, N. Mouquet, D. Gravel, and F. Jordan. 2012. The dissimilarity of species interaction networks. \emph{Ecology Letters} \bold{15}, 1353–-1361. doi: 10.1111/ele.12002

Legendre, P. 2014. Interpreting the replacement and richness difference components of beta diversity. \emph{Global Ecology and Biogeography} \bold{23}, 1324--1334. doi: 10.1111/geb.12207

Koleff, P., Gaston, K.J., and J.J. Lennon. 2003. Measuring beta diversity for presence–absence data. \emph{Journal of Animal Ecology}
\bold{72}, 367--382.
}

\author{Jochen Fründ}

\note{
  This function allows to use quantitative dissimilarity indices, which are usually recommended. However, for quantitative networks it is far from trivial how to correctly separate (or even define) which part of the dissimilarity is due to rewiring and which due to difference in species composition! Here I use the concept that all variation between the subwebs of shared species can be attributed to rewiring, BUT this will most likely not be correct. Even if all species are shared among two networks, quantitative species dissimilarity may be large (different [relative] abundances), and this will most likely lead to changes in network frequencies (changing link weights, or even missing links) that should not be called rewiring. How to correctly define and measure that is open to discussion, but thus I would still consider the values for beta_OS as overestimates, even if they are smaller in the "adjusted" (or "commondenom") than in the "poisot" case.
  
  This function should also work for one-mode networks.

  Why the name of the function? Short for "betalink revised".
  
  Why the names of the output values? These are the indices of beta (for betadiversity); I am keeping the original names and guess what they stand for: S (Species), OS (Only Shared species links), WN (Whole Network links), ST (Species Turnover links).
  
  Thanks to Carsten Dormann, Benjamin Schwarz, Benoit Gauzens, Nacho Bartomeus and Timothee Poisot for their contributions to this function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{vegdist}} and \code{\link{betadiver}} for vegan package functions calculating dissimilarity / betadiversity, used here unless \code{partitioning="commondenom"}.
  
  For reshaping web data to the array input format expected here, see \code{\link{webs2array}} and \code{\link{frame2webs}}.
}

\examples{
# two examples that give the same results as would the 
# \code{betalink} function in the package of the same name
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)))
betalinkr(webs2array(list(Safariland=Safariland, vazarr=vazarr)),function.dist="betadiver",index=1)

# another example (no shared links)
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
  lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
  webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
betalinkr(frame2webs(testdata, type.out="array"), distofempty="zero") 
# avoids NA/NaN output for the components, despite no interactions between shared species

# more than two webs:
betalinkr_multi(webs2array(Safariland, vazquenc, vazarr), index="jaccard")
}

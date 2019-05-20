\encoding{UTF-8}
\name{webs2array}
\alias{webs2array}

\title{Puts two or more webs into one array of webs}

\description{
Function to put several webs into an array, blowing the dimensions up to the union of species
}

\usage{
webs2array(...)
}


\arguments{
  \item{...}{matrices containing webs (see example), separated by comma. Alternatively, a list of webs can be given as single argument.}
}

\details{
Some analyses may require a direct comparison of two webs, e.g. computing their similarity in number of interactions per link (e.g. Poisot et al. 2012). To be able to do that, we first need to blow the single webs up to have the same dimensions, i.e. padding all species not observed in this web with 0s. This function produces a new raw matrix based on the union of species for each trophic level and puts these for all species into one array.

If the function is used on a list of webs that has no names, new names will be created in the format "web1", "web2", etc. 
}

\value{
An array of dimensions (number of species in lower level, number of species in higher level, number of webs).}

\note{
When building the webs from a table of observed interactions, function \code{\link{frame2webs}} allows to build a web-array when using \option{type.out="array"}. That would be preferable (faster and only one step)!

Combining webs into one array makes sense only for webs featuring overlapping communities! 
}

\references{ 
Poisot, T., E. Canard, D. Mouillot, N. Mouquet, D. Gravel, and F. Jordan. 2012. The dissimilarity of species interaction networks. \emph{Ecology Letters} \bold{15}, 1353–-1361. doi: 10.1111/ele.12002
 }

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} and Jochen Fründ}


\examples{
data(Safariland, vazquenc, vazquec)
allin1 <- webs2array(Safariland, vazquenc, vazquec)

# now we can compute distance between two webs:
vegdist(t(cbind(as.vector(allin1[,,1]), as.vector(allin1[,,2]), as.vector(allin1[,,3]))), 
  method="jacc")
  
# alternative input option: weblist
allin1 <- webs2array(list(Safari=Safariland, Llao=vazllao))
}
\keyword{ package }


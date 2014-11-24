\encoding{UTF-8}

\name{sortweb}
\alias{sortweb}

\title{ Function to sort bipartite webs }
\description{ This function sorts bipartite webs by either increasing/decreasing row and column totals or by a given sequence
}

\usage{
sortweb(web, sort.order="dec", sequence=NULL)
}

\arguments{
  \item{web}{ A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level 
species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{sort.order}{ sort.order can be either\cr
            \bold{inc}: sorted by increasing row/column totals\cr
            \bold{dec}: sorted by decreasing row/column totals\cr
            \bold{seq}: sorted by a given order, see sequence for how to specify an order\cr}
  \item{sequence}{list of two with two named vectors: \code{seq.lower} and \code{seq.higher}, which specify the order in which species are 
sorted. To be able to sort by names obviously the given web must be a NAMED matrix, i.e. has column and row names. If you want to order 
an unnamed web, you can either order it by row/column totals are you have to simply specify the sequence as follows (example puts last row/column to the front):
web[c(9, 1:8), c(27, 1:26)]}
}

\value{
  Returns an ordered bipartite web.
}
\references{ Vázquez, P.D., Chacoff, N.,P. and  Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology}in press.  }

\author{ Bernd Gruber }

\seealso{ For plotting and ordered web see  \code{\link{plotweb}}, \option{method="normal"} or  \code{\link{visweb}}, \option{type="none"}.}

\examples{
data(Safariland)
web <- Safariland

sortweb(Safariland, sort.order="dec") 
#rarest species first:
plotweb(sortweb(Safariland, sort.order="inc"), method="normal")   
visweb(sortweb(Safariland,sort.order="inc"), type="diagonal", 
	square="compartment", text="none", frame=TRUE)

# sorted by a given (here random) sequence
sequence <- list(seq.higher=sample(colnames(Safariland)), 
	seq.lower=sample(rownames(Safariland)))
web.ordered <- sortweb(web, sort.order="seq", sequence=sequence)

}
\keyword{ package }

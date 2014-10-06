\encoding{UTF-8}
\name{junker2013}
\alias{junker2013}
\docType{data}
\title{
Flower visitation network
}
\description{
A large (56 plant species by 257 visitor species) network published by Junker et al. (2013). 
}
\usage{data(junker2013)}
\format{
  The format is:
 int [1:56, 1:257] 0 0 0 0 0 0 1 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:56] "Achillea.millefolium" "Alliaria.petiolata" "Bellis.perennis" "Bunias.orientalis" ...
  ..$ : chr [1:257] "Agrypnus.murinus" "Ampedus.pomorum" "Anaspis.frontalis" "Anthaxia.nitidula" ...
}
\details{
I modified some entries in the table: (1) There were two instances of Prunus.sp.1.Kirsche (cherry), which I summed and represented as one. Similarly, there (2) two species named Apidae_sp._1 and Apidae_sp.1 which I merged and (3) the exact same thing for Apidae_sp._2 and Apidae_sp.2. In all cases, only one or two observations were added to the column containing more counts. I do not think that these changes will have any effect on the analyses.
}
\references{
Junker, R. R., Blüthgen, N., Brehm, T., Binkenstein, J., Paulus, J., Schaefer, H. M. and Stang, M. 2013. Specialization on traits as basis for the niche-breadth of flower visitors and as structuring mechanism of ecological networks. \emph{Functional Ecology} \bold{27}, 329–-341
}
\examples{
data(junker2013)
\dontrun{plotweb(junker2013)}
}
\keyword{datasets}

\encoding{UTF-8}
\name{Safariland}

\alias{Safariland}

\docType{data}

\title{ A pollination web from Argentina}
\description{
This pollination web was published by Vázquez and Simberloff (2003). See there for details on how it was measured.
}

\usage{data(Safariland)}


\details{
  The dataset consists of a matrix with 9 rows, representing plant species and 27 columns, representing different pollinators. Values in the matrix are observed flower visitations.

  This dataset is fairly representative of a pollination web: more pollinators than plants, and an awful lot of zeros in the matrix.

  The study was conducted in four grazed and four ungrazed sites in and around Nahuel Huapi National Park and surrounding areas in Rio Negro, Argentina from September 1999  to February 2000. For each site, the plant-pollinator interaction network was described.
  
  The authors recorded their data by counting the number of visits of each flower visitor species to each plant species. Data are presented as an interaction frequency matrix,  in which cells with positive integers indicate the frequency of interaction between  a pair of species, and cells with zeros indicate no interaction.

}

\source{
These data can be downloaded, together with the other datasets, on the NCEAS interactionweb website  \url{http://www.nceas.ucsb.edu/interactionweb/html/datasets.html}. See also there for further details, both on the data and their usage.
}

\references{
Vázquez, D. P. 2002 Interactions among Introduced Ungulates, Plants, and Pollinators: A Field Study in the Temperate Forest of the Southern Andes. Doctoral Dissertation Thesis, University of Tennessee, Knoxville, Tennessee, USA.

Vázquez, D. P., and D. Simberloff. 2002 Ecological specialization and susceptibility to disturbance: conjectures and refutations. \emph{American Naturalist} \bold{159}, 606--623.

Vázquez, D. P., and D. Simberloff. 2003 Changes in interaction biodiversity induced by an introduced ungulate. \emph{Ecology Letters}, \bold{6}, 1077--1083.
}

\examples{
data(Safariland)
plotweb(Safariland)
}

\keyword{datasets}

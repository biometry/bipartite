\name{ollerton2003}
\alias{ollerton2003}
\docType{data}
\title{ollerton2003}
\description{
 A flower visitation network from an upland grassland site in the KwaZulu-Natal reagion, South Africa
}
\usage{data(ollerton2003)}
\format{
  A data frame with 9 observations on the following 56 variables.
}
\details{
The study was conducted in the KwaZulu-Natal region of South Africa. During 3 months of fieldwork the flower visitors and pollinators to a assemblage of nine asclepiads at an upland grassland site were studied Two of the specialized pollination systems that were documented are new to the asclepiads: fruit chafer pollination and pompilid wasp pollination (the latter is almost unique in the angiosperms).

The authors recorded their data by counting the number of individua flower visitors observed and/or caught on each plant species. The total number of individuals observed on each plant species provide a rough estimate of the level of visitation that each species received. Data are presented as an interaction frequency matrix, in which cells with positive integers indicate the frequency of interaction between a pair of species, and cells with zeros indicate no interaction. 
}
\source{
NCEAS data base on interaction webs: \url{https://iwdb.nceas.ucsb.edu/resources.html#plant_pollinator} 
}
\references{
Ollerton, J., S.D. Johnson, L. Cranmer and S. Kellie (2003) The pollination ecology of an assemblage of grassland asclepiads in South Africa. \emph{Annals of Botany} \bold{92} 807-834  
}
\examples{
data(ollerton2003)
plotweb(ollerton2003)
## maybe str(ollerton2003) ; plot(ollerton2003) ...
}
\keyword{datasets}

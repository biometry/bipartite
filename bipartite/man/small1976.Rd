\name{small1976}
\alias{small1976}
\docType{data}
\title{ A flower visitation network from a peat bog in Ottawa, Canada}
\description{
The study took place in the Mer Bleue peat bog of Ottawa, Canada in 1973. The paper is a preliminary evaluation of the pollination relationships of the major entomophilous plant species of the Mer Bleue.
}
\usage{data(small1976)}

\details{
The authors recorded their data by counting the number of individual flower visitors caught on each plant species. The total number of individuals collected on each plant species provide a rough estimate of the level of visitation that each species received. Data are presented as an interaction frequency matrix, in which cells with positive integers indicate the frequency of interaction between a pair of species, and cells with zeros indicate no interaction.
}
\source{\url{https://iwdb.nceas.ucsb.edu/resources.html#plant_pollinator}}
\references{
Small, E. 1976. Insect pollinators of the Mer Bleue peat bog of Ottawa. Canadian Field Naturalist 90:22-28.
}
\examples{
data(small1976)
## maybe str(small1976) ; plot(small1976) ...
}
\keyword{datasets}

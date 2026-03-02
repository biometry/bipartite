\name{inouye1988}
\alias{inouye1988}
\docType{data}
\title{ A pollination network from the Snowy Mountains of New South Wales, Australia}
\description{
 This data set reports a community-level study of the pollination biology of alpine plants in Kosciusko National Park in the Snowy Mountains of south-eastern New South Wales, Australia. The flora and their associated insect pollinators were observed from December 1983 until March 1984.
}
\usage{data(inouye1988)}
\format{
  A data frame with 41 observations on the following 83 variables, with plant species in rows and pollinators in columns.
}
\details{
 The authors recorded their data by counting the number of individual flower visitors caught on each plant species. The total number of individuals collected on each plant species provide a rough estimate of the level of visitation that each species received. Data are presented as an interaction frequency matrix, in which cells with positive integers indicate the frequency of interaction between a pair of species, and cells with zeros indicate no interaction.
}
\note{
Male and female pollinators were summed when moving this data set from NCEAS to bipartite.
}
\source{
  NCEAS data base on interaction webs: \url{https://iwdb.nceas.ucsb.edu/resources.html#plant_pollinator}

}
\references{
  Inouye, D.W. and G.H. Pyke (1988) Pollination biology in the Snowy Mountains of Australia: comparisons with montane Colorado, USA. \emph{Australian Journal of Ecology} \bold{13}: 191--210. 
}
\examples{
data(inouye1988)
plotweb(inouye1988)
}
\keyword{datasets}

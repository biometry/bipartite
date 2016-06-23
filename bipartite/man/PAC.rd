\encoding{utf8}
\name{PAC}
\alias{PAC}

\title{Potential for Apparent Competition}
\description{
Quantifies, for each pair of lower trophic level species, the potential for showing apparent competition with another species, mediated through the higher trophic level.
}
\usage{
PAC(web)
}

\arguments{
  \item{web}{A host-parasitoid network (or alike), where the entries represent the \bold{sum of parasitoids emerging from the interactions between parasitoid and host} (i.e. number of interactions * number of parasitoid individuals emerging from each host). Only if there is only one parasitoid per host this web will be the same as that used in all other calculations in this package!}
}

\details{
Calculates the potential for apparent competition (Holt 1977), following the formula given in Müller et al. (1999) and Morris et al. (2005). See also Morris et al. (2004) for an experimental test.
}

\value{
 Returns a k x k matrix with entries d.ij, where k is the number of species in the lower trophic level and i and j are lower trophic level species. The matrix represents the effect of column species on row species. Diagonal entries are \dQuote{apparent intraspecific competition}.
}

\note{
The idea is that in host-parasitoid networks one host also affects other hosts by the number of parasitoid that hatch from it and are thus added to the pool of parasitoids. An abundant, large host can (involuntarily) contribute many parasitoids to the pool, thus also increasing the parasitoid burden of other hosts. This looks like competition between the two hosts, while in fact it is mediated through the other trophic level. 

Whether this concept can be usefully applied to mutualist networks (such as flower visitation networks, aka pollination webs) is still under-investigated. The study of Carvalheiro et al. (2014) provide an example for constructive use in mutualistic networks.
}

\references{ 
Carvalheiro, L.G., Biesmeijer, J.C., Benadi, G., Fründ, J., Stang, M., Bartomeus, I., Kaiser-Bunbury, C.N., Baude, M., Gomes, S.I.F., Merckx, V., Baldock, K.C.R., Bennett, A.T.D., Boada, R., Bommarco, R., Cartar, R., Chacoff, N., Dänhardt, J., Dicks, L. V., Dormann, C.F., Ekroos, J., Henson, K.S.E., Holzschuh, A., Junker, R.R., Lopezaraiza-Mikel, M., Memmott, J., Montero-Castaño, A., Nelson, I.L., Petanidou, T., Power, E.F., Rundlöf, M., Smith, H.G., Stout, J.C., Temitope, K., Tscharntke, T., Tscheulin, T., Vilà, M. & Kunin, W.E.  2014 The potential for indirect effects between co-flowering plants via shared pollinators depends on resource abundance, accessibility and relatedness. \emph{Ecology Letters} \bold{17}, 1389–1399.

Holt, R. D. 1977 Predation, apparent competition and the structure of prey communities. \emph{Theoretical Population Biology} \bold{12}, 197--229.

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2004 Experimental evidence for apparent competition in a tropical forest food web. \emph{Nature} \bold{428}, 310--313.

Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2005 Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449--462.

Müller, C. B., Adriaanse, I. C. T., Belshaw, R. and Godfray, H. C. J. 1999 The structure of an aphid-parasitoid community. \emph{Journal of Animal Ecology} \bold{68}, 346--370

}

\author{ Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}}

\examples{
data(Safariland)
PAC(Safariland)
}

\keyword{ package }


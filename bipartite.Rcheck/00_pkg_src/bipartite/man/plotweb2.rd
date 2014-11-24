\encoding{UTF-8}
\name{plotweb2}
\alias{plotweb2}

\title{Visualize a tripartite interaction matrix (e.g. a tritrophic foodweb)}

\description{
Two two dimensional matrix are plotted as a tripartite graph.
}

\usage{
plotweb2(web, web2, method = "cca", empty = FALSE, labsize = 1, ybig = 1,
    y_width = 0.1, spacing = 0.05, arrow="no", col.interaction="grey80",
    col.pred = "grey10", col.prey="grey10", lab.space=1, lablength = NULL, 
    sequence=NULL, low.abun=NULL,low.abun.col="green", high.abun=NULL, 
    high.abun.col="red", method2 = "cca", empty2 = TRUE, spacing2 = 0.05,  
    arrow2="no", col.interaction2="grey80", col.pred2 = "grey30", 
    col.prey2="grey20", lablength2 = NULL, sequence.pred2=NULL,low.abun2=NULL,
    low.abun.col2="green", high.abun2=NULL, high.abun.col2="red")
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{web2}{The other web to be included.}
  \item{method}{ Default method is \option{cca}, which leads to as few crossings of interactions as possible. The other option is \option{normal}, which leaves order as given by the matrix.}
  \item{empty}{logical; should empty columns or empty rows be omitted from plotting; defaults to true}
  \item{labsize}{ factor for size of labels, default is 1 }
  \item{ybig}{ vertical distance between upper and lower boxes, default is 1}
  \item{y_width}{ width of upper and lower boxes, default is 0.1 }
  \item{spacing}{ horizonatal distance between upper and lower boxes, default is 0.05}
  \item{arrow}{ display type of connection between upper and lower boxes, options are \option{up}, \option{down}, \option{both} and \option{no}, default is \option{no}, which is a polygonal connection between boxes. }
  \item{col.interaction}{color of interaction, default is grey80. }
  \item{col.pred}{color of upper boxes, default is grey10.}
  \item{col.prey}{color of lower boxes, default is grey10.}
  \item{lab.space}{sometimes it is neccessary to add additional space for labels below and above of the boxes, so all labels are shown, default is 1.}
  \item{lablength}{number of characters of labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{sequence}{list of two with two names vectors: \code{seq.pred} and \code{seq.prey}, which specify the order in which species are plotted. Cannot be set for \option{method="cca"}. Defaults to \code{NULL}, where the sequence remains as given or is determined by the CCA internally.}
  \item{low.abun}{Vector with independent abundance estimates for the lower trophic level, NULL if none exists.}
  \item{low.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to green.}
  \item{high.abun}{Vector with independent abundance estimates for the higher trophic level, NULL if none exists.}
  \item{high.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to red.}

  \item{method2}{ Default method is \option{cca}, which leads to as few crossings of interactions as possible. The other option is \option{normal}, which leaves order as given by the matrix.}
  \item{empty2}{logical; should empty columns or empty rows be omitted from plotting; defaults to true}
  \item{spacing2}{ horizonatal distance between upper and lower boxes, default is 0.05}
  \item{arrow2}{ display type of connection between upper and lower boxes, options are \option{up}, \option{down}, \option{both} and \option{no}, default is \option{no}, which is a polygonal connection between boxes. }
  \item{col.interaction2}{color of interaction, default is grey80. }
  \item{col.pred2}{color of upper boxes, default is grey10.}
  \item{col.prey2}{color of lower boxes, default is grey10.}
  \item{lablength2}{number of characters of labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{sequence.pred2}{list of two with two names vectors: \code{seq.pred} and \code{seq.prey}, which specify the order in which species are plotted. Cannot be set for \option{method="cca"}. Defaults to \code{NULL}, where the sequence remains as given or is determined by the CCA internally.}
  \item{low.abun2}{Vector with independent abundance estimates for the lower trophic level, NULL if none exists.}
  \item{low.abun.col2}{Colour for depicting the abundance estimates for the lower trophic level; defaults to green.}
  \item{high.abun2}{Vector with independent abundance estimates for the higher trophic level, NULL if none exists.}
  \item{high.abun.col2}{Colour for depicting the abundance estimates for the lower trophic level; defaults to red.}
}
\value{
  Returns a window with a tripartite graph of a food web.
}

\author{ Bernd Gruber \email{bernd.gruber@canberra.edu.au} }

\seealso{ For a different plot of food webs see  \code{\link{visweb}} and \code{\link{plotweb}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }



\encoding{UTF-8}
\name{plotweb}
\alias{plotweb}

\title{Visualize a bipartite interaction matrix (e.g. a foodweb)}
\description{
A two dimensional matrix is plotted as a bipartite graph.
}
\usage{
plotweb(web,
    method = "cca", empty = TRUE, labsize = 1, ybig = 1,  y.width.low = 0.1, 
    y.width.high = 0.1, low.spacing = NULL, high.spacing = NULL,
    arrow="no",  col.interaction="grey80", col.high = "grey10", 
    col.low="grey10",  bor.col.interaction ="black", bor.col.high="black", 
    bor.col.low="black", high.lablength = NULL, low.lablength = NULL,
    sequence=NULL, low.abun=NULL, low.abun.col="green", 
    bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
    bor.high.abun.col="black", text.rot=0, text.high.col="black", 
    text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
    low.y=0.5, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE, 
    high.plot=TRUE, high.xoff = 0, low.xoff = 0, high.lab.dis = NULL, 
    low.lab.dis = NULL, abuns.type="additional")
}

\arguments{
  \item{web}{ Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{method}{ Default method is \option{cca}, which leads to as few crossings of interactions as possible. The other option is \option{normal}, which leaves order as given by the matrix.}
  \item{empty}{logical; should empty columns or empty rows be omitted from plotting; defaults to true}
  \item{labsize}{ factor for size of labels, default is 1 }
  \item{ybig}{ vertical distance between upper and lower boxes, default is 1}
  \item{y.width.low}{ width of lower boxes, default is 0.1 }
  \item{y.width.high}{ width of upper boxes, default is 0.1 }
  \item{low.spacing}{ distance between lower boxes, default is NULL, so automatically spaced that length of upper and lower boxes is the same. Be aware if set to any value that x.lim may has to be adjusted to ensure that the network is not plotted outside the plotting region}
  \item{high.spacing}{ distance between upper boxes, default is is NULL, so automatically spaced that length of upper and lower boxes is the same. Be aware if set to any value that x.lim may has to be adjusted to ensure that the network is not plotted outside the plotting region}

  \item{arrow}{ display type of connection between upper and lower boxes, options are \option{up}, \option{down}, \option{both}, \option{up.center}, \option{down.center}, \option{both.center} and \option{no}, default is \option{no}, which is a polygonal connection between boxes. }
  \item{col.interaction}{color of interaction, default is grey80. }
  \item{col.high}{color of upper boxes, default is grey10.}
  \item{col.low}{color of lower boxes, default is grey10.}
  \item{bor.col.interaction}{border color of interaction, default is black}
  \item{bor.col.high}{border color of upper boxes, default is black}
  \item{bor.col.low}{border color of lower boxes, default is black}
  \item{high.lablength}{number of characters of upper labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{low.lablength}{number of characters of lower labels that should be plotted. If zero no labels are shown, default is NULL which plots the complete labels.}
  \item{sequence}{list of two with two names vectors: \code{seq.high} and \code{seq.low}, which specify the order in which species are plotted. Cannot be set for \option{method="cca"}. Defaults to \code{NULL}, where the sequence remains as given or is determined by the CCA internally.}
  \item{low.abun}{\bold{Named} vector with independent abundance estimates for the lower trophic level, NULL if none exists. See Notes!}
  \item{low.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to green.}
  \item{bor.low.abun.col}{border color for depicting the abundance estimates for the lower trophic level, default is black}
  \item{high.abun}{\bold{Named} vector with independent abundance estimates for the higher trophic level, NULL if none exists. See Notes!}
  \item{high.abun.col}{Colour for depicting the abundance estimates for the lower trophic level; defaults to red.}
  \item{bor.high.abun.col}{border color for depicting the abundance estimates for the higher trophic level, default is black}
  \item{text.rot}{orientation of labels in the plot (to avoid overlapping of horizontal labels if dimension of the webs are high), default is 0 for horizontal labels, use text.rot=90 for vertical labels.}
  \item{text.high.col}{Colour for text labels of higher trophic level, a vector of colours can be given}
  \item{text.low.col}{Colour for text labels of lower trophic level. A vector of colours can be given}
  \item{adj.high}{Adjust upper labels. See adj in ?text how to adjust labels }
  \item{adj.low}{Adjust upper labels. See adj in ?text how to adjust labels }
  \item{plot.axes}{axis are plotted. Sometimes useful to place boxes in higher order plots. Defaults to FALSE}
  \item{low.y}{Position of the higher boxes on the y-axis. Defaults to 1.5 }
  \item{high.y}{Position of the higher boxes on the y-axis. Defaults to 1.5 }
  \item{add}{if set to TRUE a new bipartite network is added to the previous plot. In this way multitrophic webs can be plotted, see examples below. Defaults to FALSE.}
  \item{y.lim}{To set limits for y-axis. Useful if labels are plotted outside the plotting region and for multitrophic plots, see examples below }
  \item{x.lim}{To set limits for x-axis. Useful if labels are plotted outside the plotting region and for multitrophic plots, see examples below }
  \item{low.plot}{Defines if lower boxes should be drawn. Use in multitrophic plots to avoid plotting boxes of some trophic levels - see examples below. Defaults to TRUE}
  \item{high.plot}{Defines if higher boxes should be drawn. Use in multitrophic plots to avoid plotting boxes of some trophic levels - see examples below. Defaults to TRUE}
  \item{high.xoff}{allows to set an offset to upper boxes. Useful if high.spacing is used to centre boxes manually. Use plot.axes=TRUE for easy centring}
  \item{low.xoff}{allows to set an offset to lower boxes. Useful if low.spacing is used to centre boxes manually. Use plot.axes=TRUE for easy centring}
  \item{high.lab.dis}{normally labels are staggered to avoid plotting over themselves. if set to 0, higher labels are all on one horizontal line. By using ad.low the position of the labels can be adjusted. If set to any other value labels are staggered with this distance. Defaults to NULL}
  \item{low.lab.dis}{normally labels are staggered to avoid plotting over themselves. if set to 0, lower labels are all on one horizontal line. By using ad.low the position of the labels can be adjusted. If set to any other value labels are staggered with this distance. Defaults to NULL}
  \item{abuns.type}{How to plot abundances - are they \option{independent} (e.g. flower cover) measurements or are they \option{additional} (e.g. unparasitised hosts) measurements? Defaults to \option{additional}. Option \option{none} is interpreted in the same way as \option{additional}. See Notes!}
}
\value{
  Returns a window with a bipartite graph of a food web.  For all colours vectors can be used (which are recycled if length differs. Now more trophic webs can be plotted by using \code{plotweb} and the \option{add} switch, which allows to add more webs and staggering them on top of each other. Preferred option is here to order webs by yourself and use \option{method="normal"} to keep your preferred order. See examples on three and four trophic networks.
}

\note{
%Note that in the current implementation, \option{low.abun} is actually treated as 'additional abundances' (e.g. unparasitised hosts); we are yet to find a solution to plot fully independent abundance estimates; a workaround for fully independent abundance estimates (e.g. flower densities) is to first scale \option{low.abun} such that \code{all(low.abun > rowSums(web))} and 2. then use \code{plotweb(web, low.abun= low.abun - rowSums(web))}.

Note that in previous implementations, \option{low.abun} and \option{high.abun} was actually treated as 'additional abundances' (e.g. unparasitised hosts). We added the parameter \option{abuns.type} to switch between the classic function (values \option{none} and \option{additional}, the default for backward compatibility) and an alternative function that plots fully independent abundance estimates (\option{abuns.type="independent"}). To change box (species) colors with \option{abuns.type='independent'}, use \option{low.abun.col} and \option{high.abun.col}. The function will likely be revised again in the future, so feedback welcome.

If you have total abundance measure but want to plot used and unused resources, use \option{abuns.type='additional'} and \option{low.abun= your.low.abun - rowSums(web)}.
}

\references{
Tylianakis, J. M., Tscharntke, T. and Lewis, O. T. (2007) Habitat modification alters the structure of tropical host-parasitoid food webs. \emph{Nature} \bold{445}, 202--205
 }
\author{ Bernd Gruber \email{bernd.gruber@canberra.edu.au} }

\seealso{ For a different plot of food webs see  \code{\link{visweb}} }
\examples{
data(Safariland)
plotweb(Safariland)

# shorter labels
plotweb(Safariland, high.lablength=3, low.lablength=0, arrow="down")

# centered triangles for displaying interacions
plotweb(Safariland, text.rot=90, arrow="down.center", col.interaction="wheat2",
	y.lim=c(-1,2.5))

#orginal sequence, up arrows and different box width
plotweb(Safariland, method="normal", arrow="up", y.width.low=0.3, low.lablength=4)
# interactions as lines
plotweb(Safariland, arrow="both", y.width.low=0.05, text.rot=90, col.high="blue", 
	col.low="green")

# add an abundance vector for lower trophic species 
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction="purple", 
	y.width.low=0.05, y.width.high=0.05)

plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction ="red", 
	bor.col.interaction="red", arrow="down")

# now vectors for all colours can be given, to mark certain species or 
# interactions. Colour vectors are recycled if not of appropriate length
plotweb(Safariland,col.high=c("orange","green"))
plotweb(Safariland,col.low=c("orange","green"),col.high=c("white","grey","purple"),
	text.high.col=c("blue","red"), col.interaction=c("red",rep("green",26),rep("brown",242)),
	bor.col.interaction=c(rep("green",26),rep("brown",242)),method="normal", 
	text.rot=90, low.lablength=10, high.lablength=5)


#example one (tritrophic)
plotweb(Safariland,y.width.low=0.1, y.width.high=0.05,method="normal", 
	y.lim=c(0,3), arrow="up", adj.high=c(0.5,1.5), col.high="orange",
	high.lablength=3,high.lab.dis=0)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal",
	add=TRUE,low.y=1.5,high.y=2.5, col.low="green", text.low.col="red", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1),low.plot=FALSE)

#example two (4 trophic with abundance)
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, high.abun=low.abun, col.interaction="purple", 
	y.lim=c(0,4.5), high.lablength=0, arrow="up", method="normal", 
	y.width.high=0.05)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=1.7,high.y=2.7, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4, 
	high.lablength=0)

plotweb(Safariland,y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=2.95, high.y=3.95, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4)

# now some examples with the abuns.type-option:
plotweb(Safariland, abuns.type='independent',arrow="down.center")
plotweb(Safariland, abuns.type='additional',arrow="down.center")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }



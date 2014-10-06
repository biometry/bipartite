\encoding{UTF-8}

\name{plotPAC}

\alias{plotPAC}

\title{
Function to draw a circular plot to visualise potential apparent competition (PAC)
}

\description{
Visualises the Potential for Apparent Competition as a circular graph with species represented by circles and shared parasitoids/predators/pollinators as connecting lines. Area of circles is proportional to species' abundance, filling of circles proportional to self-loops (e.g. parasitoids emerging from this species and then attacking it again) and width of connecting lines is proportional to ``export'' of parasitoids. The function is modelled after the example given in Morris et al. (2005)
}

\usage{
plotPAC(web, scaling = 1, plot.scale = 1, fill.col = rgb(0.2, 0.2, 0.2, 0.5), 
arrow.col = rgb(0.5, 0.5, 0.5, 0.5), outby = 1, text = TRUE, 
circles = FALSE, radius = 1)
}
\arguments{
  \item{web}{A community matrix with lower trophic level in rows and higher trophic level in columns.}
  \item{scaling}{A factor scaling the size of species-circles. The default of 1 may cause overlap when there are many species or some are very large. Smaller values yield smaller circles.}
  \item{plot.scale}{Scales the whole plot to the plotting region. If you want to add labels, you may want more space around the plot and hence choose smaller values for plot.scale.}
  \item{fill.col}{Colour of the species-circles. Note that the default is using transparency, which is the fourth parameter in the \code{rgb}-function.}
  \item{arrow.col}{Colour of the arrows (in fact polygons) connecting species. As for \option{fill.col}, the default uses transparency.}
  \item{outby}{A factor determining by how much the text labels should be moved out from the species-circles. Values smaller than 1 will move them inside the circle-plot.}
  \item{text}{Logical; if TRUE (default), a number will be plotted next to each species-circle. There usually is no space to put the full name next to it, so there is currently no option to do so.}
  \item{circles}{Logical; shall species labels be put into a circle (as in the original plot of Morris et al. 2005)? Defaults to FALSE.}
  \item{radius}{A factor modifying the size of the label-circles.}
}

\details{
The function is modelled after Morris et al. (2005). The whole idea and application is explained there, too. The only change to their plotting is the default choice of transparency for clarity of the visual effect.
}

\value{
None. The function is invoked for its side effects (i.e. printing).
}
\references{
Morris, R. J., Lewis, O. T. and Godfray, H. C. J. 2005 Apparent competition and insect community structure: towards a spatial perspective. \emph{Annales Zoologica Fennici} \bold{42}, 449–462.
}

\author{
Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}
}

\seealso{
See also comments and notes in \code{\link{PAC}}!
}

\examples{
\dontrun{
data(kevan1970)
plotPAC(kevan1970)

data(Safariland)
plotPAC(Safariland, plot.scale=1, fill.col="red", arrow.col="orange", 
	circles=TRUE, radius=1)
}
}

\keyword{ package }

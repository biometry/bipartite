\encoding{UTF-8}
\name{plotweb_v2}
\alias{plotweb_v2}

\title{Visualize a bipartite interaction matrix (e.g. a foodweb)}
\description{
A two dimensional matrix is plotted as a bipartite graph. 
}
\usage{
plotweb_v2 <- function(web,
                       sorting = "normal",
                       empty = FALSE,
                       upper_abundances = NULL,
                       lower_abundances = NULL,
                       add_upper_abundances = NULL,
                       add_lower_abundances = NULL,
                       scaling = "relative",
                       font = NULL,
                       family = NULL,
                       srt = 0,
                       text_size = 1,
                       box_size = 0.1,
                       upper_italic = FALSE,
                       lower_italic = FALSE,
                       x_lim = c(0, 1),
                       y_lim = c(0, 1),
                       lab_distance = 0.05,
                       upper_color = "black",
                       upper_border = "same",
                       upper_add_color = "red",
                       higher_text_color = "black",
                       lower_color = "black",
                       lower_border = "same",
                       lower_add_color = "red",
                       lower_text_color = "black",
                       horizontal = FALSE,
                       abbr_names = FALSE,
                       link_color = "upper",
                       link_border = "same",
                       link_alpha = 0.5,
                       style = "line",
                       arrow = "no",
                       spacing = "auto",
                       plot_axes = FALSE,
                       add = FALSE,
                       mar = c(1, 1, 1, 1),
                       mai = NULL)
}

\arguments{
  \item{web}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). 
             Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{higher_abundances}{Vector with independent abundance estimates for the upper trophic level, NULL if none exist. See Notes!}
  \item{lower_abundances}{Vector with independent abundance estimates for the lower trophic level, NULL if none exist. See Notes!}
  \item{add_higher_abundances}{Vector with additional abundances (e.g. unparasitised hosts) for the upper trophic level, NULL if none exist. See Notes!}
  \item{add_lower_abundances}{Vector with additional abundances (e.g. unparasitised hosts) for the lower trophic level, NULL if none exist. See Notes!}
  \item{scaling}{Choice of either "relative" or "absoulte". If set to "relative" both }
  \item{font}{}
  \item{family}{}
  \item{srt}{Orientation of the labels in degrees. If horizontal = TRUE an automatic offset of 90 is added. (Defualt: 0)}
  \item{text_size}{Factor by which the size of the label text should be scaled. Values < 1 makes them smaller, values > 1 makes them greater. (Default: 1)}
  \item{box_size}{Size of the boxes. Numerical with either c(upper_box_size, lower_box_size) or one value if both should be the same. (Default: 0.1)}
  \item{upper_italic}{If set to TRUE, the labels of the upper species are plotted in italics. (Default: FALSE)}
  \item{lower_italic}{If set to TRUE, the labels of the lower species are plotted in italics. (Default: FALSE)}
  \item{x_lim}{Vector with 2 numerical values. The range in which to plot on the x-axis. (Default: c(0, 1))}
  \item{y_lim}{Vector with 2 numerical values. The range in which to plot on the y-axis. (Default: c(0, 1))}
  \item{lab_distance}{Distance between the boxes and the labels in \strong{inches}}
  \item{higher_color}{A single color or a list of colors used to plot the higher species (columns) boxes. (Default: "black")}
  \item{higher_border}{Either "same" or a single color or a list of colors used to plot the higher species (columns) boxes borders.
                      If "same", the borders have the same color as the corresponding boxes.}
  \item{higher_add_color}{Only used when higher additional abundances are defined. 
                          A single color or a list of colors used to plot the higher (column) additional boxes. (Default: "red")}
  \item{lower_color}{A single color or a list of colors used to plot the lower species (rows) boxes. (Default: "black")}
  \item{lower_border}{Either "same" or a single color or a list of colors used to plot the lower species (rows) boxes borders.
                      If "same", the borders have the same color as the corresponding boxes.}
  \item{lower_add_color}{Only used when lower additional abundances are defined. 
                         A single color or a list of colors used to plot the lower (row) additional boxes. (Default: "red")}
  \item{horizontal}{Boolean indicating whether the connection between the boxes and the label texts should be plotted horizontally. (Default: FALSE)}
  %\item{abbr_names}{}
  \item{link_color}{"lower", "upper", or single color value. If "lower" or "upper" the link color is taken from the corresponding lower or upper box. Otherwise all links have the defined color. (Default: "upper")}
  \item{link_border}{Either "same" or a single color or a list of colors used to plot the link borders.
                      If "same", the borders have the same color as the corresponding boxes.}
  \item{link_alpha}{Number between 0.0 and 1.0 setting the transparency of the links. 
                    0 means full transparency 1 means full opacity. (Default: 0.5)}
  \item{style}{"line" or "curve" indicating whether the links should be straight or curved.}
  \item{arrow}{Display type of connection between upper and lower boxes, options are 
    \itemize{
      \item "no" (Default)
      \item "up"
      \item "down"
      \item "both"
      \item "up.center"
      \item "down.center"
      \item "both.center"
    }
  }
  \item{spacing}{Either "auto", a single value or a numerical vector of the form c(lower, upper). Defines the total space between boxes. See Notes!  (Default: "auto")}
  \item{plot_axes}{Indicates whether both axes should be drawn on the plot. (Default: FALSE)}
  \item{add}{If set to TRUE a new bipartite network plot is added to the previous plot. (Default: FALSE)}
  \item{mar}{A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. 
             The default is NULL which takes the global setting from par(). }
  \item{mai}{A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
             Please note, that the margin is increased so that the labels are always visible in the plot.}
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
\author{ Tobias Bauer \email{tobias.bauer-2@uni-hamburg.de} }

\seealso{ For a different plot of food webs see  \code{\link{visweb}} }
\examples{
data(Safariland)
plotweb_v2(Safariland)

# To let the plot look a bit more modern
plotweb_v2(Safariland, horizontal = TRUE, style = "curve")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }



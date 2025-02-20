\encoding{UTF-8}
\name{plotweb}
\alias{plotweb}

\title{Visualize a bipartite interaction matrix (e.g. a foodweb)}
\description{
This is the help file for the new plotweb implementation. 
For help with the old, deprecated plotweb function see \link{plotweb_deprecated}. 

A two dimensional matrix is plotted as a bipartite graph.
}
\usage{
plotweb <- function(web,
                    sorting = "normal",
                    empty = FALSE,
                    higher_abundances = NULL,
                    lower_abundances = NULL,
                    add_higher_abundances = NULL,
                    add_lower_abundances = NULL,
                    scaling = "relative",
                    font = NULL,
                    family = NULL,
                    srt = 0,
                    text_size = 1,
                    box_size = 0.1,
                    higher_italic = FALSE,
                    lower_italic = FALSE,
                    x_lim = c(0, 1),
                    y_lim = c(0, 1),
                    lab_distance = 0.05,
                    higher_color = "black",
                    higher_border = "same",
                    higher_add_color = "red",
                    higher_text_color = "black",
                    lower_color = "black",
                    lower_border = "same",
                    lower_add_color = "red",
                    lower_text_color = "black",
                    horizontal = FALSE,
                    abbr_names = FALSE,
                    link_color = "higher",
                    link_border = "same",
                    link_alpha = 0.5,
                    curved_links = FALSE,
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
  \item{higher_abundances}{Vector with independent abundance estimates for the higher trophic level, NULL if none exist. See Notes!}
  \item{lower_abundances}{Vector with independent abundance estimates for the lower trophic level, NULL if none exist. See Notes!}
  \item{add_higher_abundances}{Vector with additional abundances (e.g. unparasitised hosts) for the higher trophic level, NULL if none exist. See Notes!}
  \item{add_lower_abundances}{Vector with additional abundances (e.g. unparasitised hosts) for the lower trophic level, NULL if none exist. See Notes!}
  \item{scaling}{Choice of either "relative" or "absoulte". For more information see Scaling.}
  \item{font}{See \link[graphics]{par}}
  \item{family}{See \link[graphics]{par}}
  \item{srt}{Orientation of the labels in degrees. (Default: 0)}
  \item{text_size}{Factor by which the size of the label text should be scaled. Values < 1 makes them smaller, values > 1 makes them greater. (Default: 1)}
  \item{box_size}{Size of the boxes. Numerical with either c(higher_box_size, lower_box_size) or one value if both should be the same. (Default: 0.1)}
  \item{higher_italic}{If set to TRUE, the labels of the higher species are plotted in italics. (Default: FALSE)}
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
  \item{link_color}{"lower", "higher", or single color value. If "lower" or "higher" the link color is taken from the corresponding lower or higher box. Otherwise all links have the defined color. (Default: "higher")}
  \item{link_border}{Either "same" or a single color or a list of colors used to plot the link borders.
                      If "same", the borders have the same color as the corresponding boxes.}
  \item{link_alpha}{Number between 0.0 and 1.0 setting the transparency of the links. 
                    0 means full transparency 1 means full opacity. (Default: 0.5)}
  \item{curved_links}{If FALSE links are plotted as straight lines otherwise they are curved. (Default: FALSE)}
  \item{arrow}{Display type of connection between higher and lower boxes, options are 
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
  \item{spacing}{Either "auto", a single value or a numerical vector of the form c(lower, higher). Defines the total space between boxes. See Notes!  (Default: "auto")}
  \item{plot_axes}{Indicates whether both axes should be drawn on the plot. (Default: FALSE)}
  \item{add}{If set to TRUE a new bipartite network plot is added to the previous plot. (Default: FALSE)}
  \item{mar}{A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. 
             The default is NULL which takes the global setting from par(). }
  \item{mai}{A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
             Please note, that the margin is increased so that the labels are always visible in the plot.}
}
\value{
  Returns a window with a bipartite graph of the given web. 
}

\section{Abundances}{
  Both independent as well as additional abundances should be passed as a named vector. 
  
  If no independent abundances are given the row and column sums are used to calculate the size of the boxes. 
}

\section{Scaling}{
  The scaling option is only relevant, if either independent or additional or both have been defined.

  In case \code{"relative"} is selected the boxes of higher and lower abundances are scaled so that they both utilize the whole space. 

  In case \code{"absolute"} is selected only the boxes of the side with the higher total abundances uses the full space the other side is then scaled relatively to that side. 
  So that species with the same abundance are plotted with boxes of the same size.

  Be extra careful about the \code{spacing} option when using \code{"absolute"} scaling. 
}

\section{Spacing}{
  The \code{"auto"} spacing option can be very handy for quickly creating a plot. 
  However in certain cases it might not produce the best results. 
  Especially in cases where there are plenty species with rather long names or multiple species with low abundances, 
  the labels can tend to overlap. 
  In those cases defining the \code{spacing} and \code{text_size} manually or rotating the labels with the \code{srt} argument might produce better results.

  Special care must be taken if two different values are passed to the argument for the higher and lower species. 
  The total size of the boxes will not be comparable between the sides in that case. 
  So especially using the argument \code{scaling = "absolute"} might not yield the intended result.
}

\note{
Note the difference between relative and absolute scaling is mainly important in the case either additional or independent abundances have been defined.
In case relative is selected the boxes of higher and lower abundances are scaled so that they both utilize the whole space.
In case absolute is selected only the boxes of the side with the higher total abundances uses the full space the other side is then scaled relatively to that side. So that species with the same abundance are plotted with boxes of the same size.
}

\author{ Tobias Bauer \email{tobias.bauer-2@uni-hamburg.de} }

\seealso{ 
  For a different plot of food webs see  \code{\link{visweb}}

  For a comparison of 2 webs see \code{\link{plot2webs}}
}
\examples{
data(Safariland)
plotweb(Safariland)

# To let the plot look a bit more modern
plotweb(Safariland, horizontal = TRUE, style = "curve")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }



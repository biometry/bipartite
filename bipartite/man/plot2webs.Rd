\encoding{UTF-8}
\name{plot2webs}
\alias{plot2webs}

\title{Visualize two bipartite interaction matrices (e.g. foodwebs at two locations)}
\description{
Two 2-D matrices are plotted as a two-layer bipartite graph.  
}
\usage{
plot2webs <- function(web1,
                      web2,
                      middle = "higher",
                      sort_by = 1,
                      sorting = "normal",
                      higher_abundances = NULL,
                      lower_abundances = NULL,
                      add_higher_abundances = NULL,
                      add_lower_abundances = NULL,
                      scaling = "absolute",
                      font = NULL,
                      family = NULL,
                      srt = 0,
                      higher_italic = FALSE,
                      lower_italic = FALSE,
                      text_size = 1,
                      x_lim = c(0, 1),
                      y_lim = c(0, 1),
                      lower_color = "black",
                      lower_border = "same",
                      add_lower_color = "same",
                      higher_color = "black",
                      higher_border = "same",
                      add_higher_color = "same",
                      higher_text_color = "black",
                      lower_text_color = "black",
                      horizontal = FALSE,
                      abbr_names = FALSE,
                      link_color = "lower",
                      link_border = "same",
                      link_alpha = 0.6,
                      curved_links = FALSE,
                      arrow = "no",
                      spacing = "auto",
                      space_lower = 0.2,
                      space_higher = 0.2,
                      mar = c(1, 1, 1, 1),
                      mai = NULL,
                      plot_axes = FALSE)
}

\arguments{
  \item{web1}{Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). 
              Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
  \item{web2}{The second layer (e.g. temporal or spatial) web. Common species between the webs must have the exact same name.}
  \item{middle}{"higher" or "lower". Defining which trophic level should be plotted in the middle of the graph.}
  \item{sort_by}{}
  \item{sorting}{}
  \item{higher_abundances}{}
  \item{lower_abundances}{}
  \item{add_higher_abundances}{}
  \item{add_lower_abundances}{}
  \item{scaling}{}
  \item{font}{}
  \item{family}{}
  \item{srt}{}
  \item{higher_italic}{}
  \item{lower_italic}{}
  \item{text_size}{}
  \item{x_lim}{}
  \item{y_lim}{}
  \item{lower_color}{}
  \item{lower_border}{}
  \item{add_lower_color}{}
  \item{higher_color}{}
  \item{higher_border}{}
  \item{add_higher_color}{}
  \item{higher_text_color}{}
  \item{lower_text_color}{}
  \item{horizontal}{}
  \item{abbr_names}{}
  \item{link_color}{}
  \item{link_border}{}
  \item{link_alpha}{}
  \item{curved_links}{}
  \item{arrow}{}
  \item{spacing}{}
  \item{space_lower}{}
  \item{space_higher}{}
  \item{mar}{}
  \item{mai}{}
  \item{plot_axes}{}
}

\value{
  Returns a window with the two bipartite graphs aligned on top or next to each other.
}

\note{
}

\author{ Tobias Bauer \email{tobias.bauer-2@uni-hamburg.de} }

\seealso{ 
  For a different plot of food webs see  \code{\link{visweb}}

  For plotting a single web see \code{\link{plotweb}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }



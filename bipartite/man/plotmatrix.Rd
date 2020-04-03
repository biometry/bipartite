\name{plotmatrix}
\alias{plotmatrix}

\title{ Plot a matrix organised by topology }

\description{
Plotmatrix is a function to plot binary and weighted matrices
}

\usage{plotmatrix(x, background_color="white", base_color=NULL, between_color="black", 
  border_color="black", modules_colors=NULL, within_color = "black", border = FALSE, 
  row_partitions=NULL, col_partitions=NULL, binary=TRUE, plot_labels=FALSE, xlab=NA, 
  ylab=NA, offset = 0.4, ...)}

\arguments{
  \item{x}{the first argument is an interaction matrix (rows and columns are nodes; cells are links between nodes) or the list returned by \code{\link{sortweb}}. The matrix may be binary (only 0s and 1s) or weighted.}
  \item{background_color}{color of the background.}
  \item{base_color}{define the base color for the gradient in weighted matrices. If \option{NULL}, defaults to \option{background_color}.}
  \item{between_color}{color of links in the regions between modules. If the matrix is weighted, links between modules are colored following a gradient from the \option{base_color} to the \option{between_color}. Only applied if partitions are provided.}  
  \item{border_color}{color of the border plotted around modules (only if \option{border = TRUE}).}
  \item{modules_colors}{optional vector with separate colors for modules. Its lengths must match the number of partitions. If the matrix is weighted, links in each module are colored following a gradient from the \option{base_color} to the colors in \option{modules_color}. Only applied if partitions are provided.}
  \item{within_color}{color of links in the regions within modules. If the matrix is weighted, links within modules are colored following a gradient from the \option{base_color} to the \option{within_color}. If partitions are not provided, defines the color of links for the entire network.}
  \item{border}{logical; if \option{TRUE}, a border is plotted around each module. Only applied if partitions are provided.}
  \item{row_partitions}{optional vector with partitions for rows. The length of the vector must be the number of rows in the matrix, each value is the partition of the respective row, following the row sequence in the matrix. Partition can be defined by numeric or character values.}
  \item{col_partitions}{same as \option{row_partitions}, but for columns.}
  \item{binary}{logical; plot the binary or weighted information of the matrix. If the matrix is binary, must be \option{TRUE}. If the matrix is weighted and binary is \option{TRUE}, plot the binary structure of the matrix.}
  \item{plot_labels}{logical: shall row and column names of the matrix be plotted?}
  \item{xlab}{label for the column axis}
  \item{ylab}{label for the row axis.}
  \item{offset}{defines the size of each filled cell, compared to the border of each cell. Values higher than 0.4 may result in overlapping of cells.}
  \item{...}{arguments passed on to \code{axis} ornamenting the image plot.}
}

\value{Invoked for its side effect of plotting the matrix.}

\author{Rafael Pinheiro & Gabriel Felix, https://rbppinheiro.wordpress.com/, e-mail: rafael-bpp@hotmail.com}

\note{Taken, with permission, from https://github.com/pinheirorbp/Rfunctions1 }

\seealso{The output of \code{sortmatrix} is used by \code{\link{plotmatrix}}.}

\examples{
S <- sortmatrix(Safariland, topology = "nested", sort_by = "weights")
plotmatrix(S)
plotmatrix(S$matrix, binary=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
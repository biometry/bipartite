\name{sortmatrix}
\alias{sortmatrix}

\title{ Organise matrix by topology }

\description{
Sortmatrix is used to organize matrix rows and columns in order to highlight one of the three available topologies: nested, modular, and compound
}

\usage{sortmatrix(matrix, topology="compound", sort_by="degrees", row_partitions=NULL, 
      col_partitions=NULL, mod_similarity=FALSE)}

\arguments{
  \item{matrix}{the interaction matrix (rows and columns are nodes; cells are links between nodes). The matrix may be binary (only 0s and 1s) or weighted.}
  \item{topology}{defines for which topology the visualization is optimized. We implemented 3 topologies: nested, modular, and compound (modular with internally nested modules). Defaults to \option{"nested"}. If \option{topology = "nested"}, network is sorted by decreasing marginal totals. If \option{topology = "modular""}, network is sorted based on row and column partitions. If \option{topology = "compound"}, network is first sorted on row and column partitions, then sorted by decreasing marginal totals within each partition.}
  \item{sort_by}{defines whether the marginal totals used for sorting are the binary marginal totals (\option{sort_by = "degrees"}) or the weighted marginal totals (\option{sort_by = "weights"}). If the matrix is binary, sort_by must be set to \option{"degrees"}, which is also the default.}
  \item{row_partitions}{vector with row partitions to be applied in modular and compound topologies. The length of the vector must be the number of rows in the matrix, each value is the partition of the respective row, following the row sequence in the matrix. Partition can be defined by numeric or character values. Not used if topology is \option{"nested"}.}
  \item{col_partitions}{same as row_partitions, but for columns.}
  \item{mod_similarity}{logical; if \option{mod_similarity=TRUE}, the order of modules is defined based on the similarity between them. So that more similar modules are closer in the plot. Similarity is based on Euclidian distances.}
}

\value{The output of sortmatrix is a list with 5 elements:
 \item{matrix}{is the sorted matrix}
 \item{row_partitions}{a vector with the partitions for rows, in the new order}
 \item{col_partitions}{a vector with the partitions for columns, in the new order}
 \item{order_row}{a vector with the position of each row of the input matrix in the sorted matrix. order_col: a vector with the position of each column of the input matrix in the sorted matrix}
 
}

\details{More than a visualization choice, the use of binary or weighted marginal totals is related to the definition of weighted nestedness. One of the most used index, WNODF (Almeida-Neto and Ulrich 2011), requires binary nestedness to account for weighted nestedness. If this index and the definition of nestedness it reflects, is used, it may be more logical to define sort_by as degrees. Other indices, as WNODA (Pinheiro et al. 2019) do not require binary nestedness, thus, there is no reason to sort by degrees rather than weights.}

\references{
Almeida-Neto, M. and Ulrich, W. 2011. A straightforward computational approach for measuring nestedness using quantitative matrices. \emph{Environ. Model. Softw.} \bold{26}, 173--178

Lewinsohn, T.M. et al. 2006. Structure in plant-animal interaction assemblages. \emph{Oikos} \bold{113},  174--184

Pinheiro, R.B.P. et al. 2019. A new model explaining the origin of different topologies in interaction networks. \emph{Ecology} \bold{100}, 1--30
}

\author{ Rafael Pinheiro & Gabriel Felix, https://rbppinheiro.wordpress.com/, email: rafael-bpp@hotmail.com}

\note{Taken, with permission, from https://github.com/pinheirorbp/Rfunctions1 }

\seealso{The output of \code{sortmatrix} is used by \code{\link{plotmatrix}}.}

\examples{
sortmatrix (Safariland, topology = "nested", sort_by = "weights")

# see example in help for "plotmatrix"
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ networks }
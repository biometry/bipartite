\encoding{UTF-8}
\name{nest.smdm}
\alias{nest.smdm}
\alias{NODF}
\alias{WNODF}
\alias{WNODA}
\alias{module2constraints}
\title{ Computes nestedness of a matrix as WNODA (and NODF and WNODF) }
\description{
  Calculates network nestedness, also within and between modules, i.e. separate nestednesses for nodes belonging to the same module and between nodes belonging to different modules. Three nestedness metrics are implemented in the function: NODF, WNODF and WNODA.
}
\usage{
  nest.smdm(x, constraints=NULL, weighted=FALSE, decreasing="fill", sort=TRUE)
  module2constraints(mod)
}
\arguments{
  \item{x}{an interaction matrix; typically with rows for lower and columns for higher level species;}
  \item{constraints}{a vector with modules for vertices of the matrix. The vector indicates first rows modules and then columns modules, following the sequence of input matrix x. If no constraints are provided, the function calculates nestedness for the entire matrix, not taking into account possible network modules. }
\item{weighted}{Logical. Indicate whether to calculate binary or weighted version of the metrics. If set to FALSE for a weighted input matrix, the binary metric is calculated, but a warning is returned.}
\item{decreasing}{The matrix has to be sorted before computation of (W)NODA. This can be done either by \option{"fill"}, i.e. sum of cells with non-zero values, or \option{"abund"}, i.e. sum of cell values. For a binary matrix decreasing needs to be \option{"fill"} (as no abundances are available). For weighted matrices, the argument may be \option{"fill"} (yielding WNODF) or \option{"abundance"} (yielding WNODA).}
\item{sort}{Logical. Should columns and rows of the matrix be sorted, in order to maximize nestedness index?}
\item{mod}{Output of a \code{\link{computeModules}}-object. Returns a vector to be used as input for \option{constraints}.}
}

\value{Function returns a list with elements
  \item{WNODArow}{Nestedness for rows}
  \item{WNODAcol}{Nestedness for columns}
  \item{WNODAmatrix}{Nestedness for entire matrix}
  If constraints are provided, e.g. based on \code{\link{computeModules}}, output additionally includes:
  \item{WNODA_SM_row}{Nestedness for rows belonging to the same modules}
  \item{WNODA_DM_row}{Nestedness for rows belonging to different modules}
  \item{WNODA_SM_col}{Nestedness for columns belonging to the same modules}
  \item{WNODA_DM_col}{Nestedness for columns belonging to different modules}
  \item{WNODA_SM_matrix}{Nestedness for nodes (rows and columns) belonging to the same modules}
  \item{WNODA_DM_matrix}{Nestedness for nodes (rows and columns) belonging to different modules}
}

\author{ Rafael Barros Pereira Pinheiro \email{rafael-bpp@hotmail.com} }

\references{

Almeida-Neto M, Guimaraes PR, Guimaraes PR Jr, Loyola RD, Ulrich W (2008) A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}: 1227--1239

Almeida-Neto, M. & Ulrich, W. (2011). A straightforward computational approach for measuring nestedness using quantitative matrices. \emph{Environ. Model. Softw.} \bold{26}: 173--178

Felix, G.M., Pinheiro, R.B.P., Poulin, R., Krasnov, B.R. & Mello, M.A.R. (2017). The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. \emph{bioRxiv}

Flores, C.O., Valverde, S. & Weitz, J.S. (2013). Multi-scale structure and geographic drivers of cross-infection within marine bacteria and phages. \emph{ISME J.} \bold{7}: 520--532
}

\seealso{\code{vegan::nestedNODF}, \code{\link{computeModules}}}

\examples{
nest.smdm(Safariland)
nest.smdm(Safariland, weighted=TRUE)
nest.smdm(Safariland, weighted=TRUE, decreasing="abund")
nest.smdm(Safariland, weighted=TRUE, decreasing="abund", sort=FALSE)
# identify modules using computeModules:
mod <- computeModules(Safariland)
const <- module2constraints(mod)
nest.smdm(Safariland, constraint=const)
nest.smdm(Safariland, constraint=const, weighted=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ networks }
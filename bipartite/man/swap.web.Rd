\encoding{UTF-8}

\name{swap.web}

\alias{swap.web}

\title{ Creates null model for bipartite networks}

\description{
Function to generate null model webs under the following constraints: 1. marginal 
totals are identical to those observed (as in \code{r2dtable}), 2. connectance is as observed (as in \code{shuffle.web}.)
}

\usage{swap.web(N, web, verbose=FALSE, c.crit=1e4)}

\arguments{
  \item{N}{ Number of desired null model matrices. }
  \item{web}{An interaction matrix.}
  \item{verbose}{Should various verbal outputs of this function be shown? Defaults to FALSE, since it was
  mainly used during the debugging period.}
  \item{c.crit}{Sometimes the algorithm gets stuck in a very sparse matrix. Then \option{c.crit} sets the number of
  swaps it shall attempt before giving up and starting over on a new random matrix. Defaults to 10000.}
}

\details{
 This function is designed to behave similar to \code{r2dtable}, i.e. it returns a list of randomised matrices.  In addition to \code{r2dtable} is also keeps the connectance constant! 
 
 This function is thought of as a very constrained null model for the analysis of bipartite webs.  It keeps two web properties constant: The marginal totals (as in \code{r2dtable} and the  number of links (and hence connectance). A comparison of \code{swap.web}- and \code{r2dtable}-based  webs allows to elucidate the effect of evolutionary specialisation, since the unrealised connections may represent \dQuote{forbidden links}. 
 
 This null model is similar to the one employed by Vázquez et al. But while Vázquez starts by assigning 1s to the  allowed connections and then fills the web, \code{swap.web} starts with an \code{r2dtable}-web and successively  \dQuote{empties} it. The two approaches should result in very similar null models, since both constrain marginal  totals and connectance.
 
 A few words on the way \code{swap.web} works. It starts with a random web created by \code{r2dtable}. Then, it finds, randomly, 2x2 submatrices with entries all larger than 0. Next, it subtracts the minimum value from the diagonal and adds it to the off-diagonal (minor diagonal). Thereby one cell becomes 0, but the column and row sums do not change. This idea is adapted from the swap-algorithm used in various binary null models by Nick Gotelli. If the random web has  too few 0s (which is I have yet to encounter), then the opposite strategy is applied.
 
 The algorithm in our implementation has some variations on finding the submatrix and constraining the number of  unsuccessfull trials before starting on a new random matrix, but they are only for speeding up the process.
 
 
}

\value{
 A list of N randomised matrices with the same dimensions as the initial web.
}

\references{
 Artzy-Randrup, Y., and Stone, L. (2005) Generating uniformly distributed random networks. \emph{Physical Review E} \bold{72}, 1--7

 Miklós, I. and Podani, J. (2004) Randomization of presence-absence matrices: comments and new algorithms. \emph{Ecology} \bold{85}, 86--92
 
 Vázquez, D. P., and M. A. Aizen (2003) Null model analyses of specialization in plant-pollinator interactions. \emph{Ecology} \bold{84}, 2493--2501

 Vázquez, D. P., C. J. Melian, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin (2007) Species abundance and asymmetric interaction strength in ecological networks. \emph{Oikos} \bold{116}, 1120--1127

 For a very nice and thorough overview of null models in general see:

 Gotelli, N. J., and G. R. Graves (1996) Null Models in Ecology. Smithsonian Institution Press, Washington D.C.
}

\author{ Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de> }

\note{ 
 Long stories can be told about the swap algorithm. I am not the right person to do so, but for a much more detailed coverage of the topic, for many more ways to implement null models for \bold{binary} matrices, with various flavours of the swap and possible alternatives, first brew yourself a cup of tea and then check out the help pages of \code{simulate} in \pkg{vegan}. (As usual, Jari Oksanen has spend considerable care to implement even the most bizarre and abstruse way to move 0s and 1s around. His ecological advise between the lines make his package worthwhile already! I, personally, would use \option{method="quasiswap"}, as is done in the example to \code{\link{discrepancy}}.)
  
 When comparing the \code{swap.web} algorithm with that proposed by Vázquez et al. (2007, implemented in \code{\link{vaznull}}), we found that \code{swap.web} contains a certain bias. The subtraction of the swap will reduce the value of low-value cells, and increase that of high-value cells. As a consequence, it produces somewhat of a dichotomy between very high and very low values. Using e.g. H2' to quantify this pattern, \code{swap.web} will produce very specialised networks (around 0.5), while the Vazquez-algorithm yields lower H2' values and a more even distribution of interactions within cells.
 The ramifications are that \code{swap.web} will predict higher-than-necessary expectations. (Date of this entry: 15.1.2010)
 
 In fact, Artzy-Randrup & Stone (2005) have shown that the swap algorithm is \emph{fundamentally biased}, because some swaps are more likely than others. This applies to this version of the swap as well as to the one implemented in \code{\link{vaznull}}. So, despite heavy citation, the approach of Miklós & Podani (2004) is thus also not ideal, as often claimed.

  \code{swap.web} is a very constraint null model. You need to consider if it is the right one for your application!
}

\seealso{ \code{\link{r2dtable}}, \code{\link{vaznull}} and \code{\link{shuffle.web}}}

\examples{

swap.web(Safariland, N=2)

}
\keyword{ package }


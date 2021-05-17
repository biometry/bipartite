\encoding{UTF-8}
\name{frame2webs}

\alias{frame2webs}
\alias{edgelist2webs}

\title{Converts a table of observations into a network matrix}

\description{
Convenience function to convert a table of observations (i.e. an "edge list", typically compiled in a spreadsheet programm) into a network matrix for further use in bipartite.
}

\usage{
frame2webs(dframe, varnames = c("lower", "higher", "webID", "freq"), type.out = 
"list", emptylist = TRUE)
}

\arguments{
  \item{dframe}{table (i.e. data.frame) of observations, typically simply the data read by one of the \code{read.*} functions; see Details for details!}
  \item{varnames}{a vector of characters giving the column names in the table (data.frame) that correspond to
  	\itemize{
      \item species names in lower trophic level (e.g. plants),
      \item species names in higher trophic level (e.g. pollinators),
      \item the grouping factor for a web (e.g. the site on which interactions were observed)
      \item optionally, a fourth column indicating the number of times an interaction was observed (or another link weight that can be summed);
    }
    If only three varnames are provided, it is assumed that each observed interaction is in a separate row and frequencies are calculated from the number of rows per link. By default, the names \dQuote{lower}, \dQuote{higher}, \dQuote{webID} and \dQuote{freq} are used.}
  \item{type.out}{defines the type of output; it could be \option{list} (default) or an \option{array} with a separate slice for each network. Lists have the advantage that different webs do not have to include all species names, i.e. they can be of different dimensions (ragged). As such they are better suited for webs with non-comparable species sets. An array has the advantage that dimensions are the same, and therefore (i) comparisons can be easily made among webs and (ii) webs can be pooled more easily. The \option{array}-option is more suitable for sets of networks from the same community.}
  \item{emptylist}{logical, defaults to \option{TRUE}; should, for each network in the list of networks, empty columns and rows be deleted? Since this function first builds an array and from that a list, each network contains all possible links between species of the higher and lower level. When set to \option{TRUE}, all species not observed for a given network are removed.}
}

\details{
  This function supports the easy handling of typical recording that are used to compose a network. The assumed structure of the underlying table is two columns for the names of the lower and higher level species, respectively, one column for a network ID (e.g. site, observer or any other grouping code) and, optionally, a number that indicates how often this interaction was observed. If not given, the number of interactions is simply computed from the number of times the same interaction occurred for each network ID. See example code below. 
  
  Typically, data are recorded (in a field book or external data logging device), read into a spreadsheet software (such as MS Excel or Open/LibreOffice Spreadsheet), where names are checked, typos corrected and so forth. This table is then imported into R and aggregated into one or more webs (or an array or a list of webs) using \code{frame2webs}.
  
  Each link can have multiple entries or a single entry!
}

\value{
  A list or an array of networks, each to be used as input for other \pkg{bipartite} functions.
}

\note{
Great care should be taken in checking species names, slightly different spellings will be assumed to be different species! Common problems occur because a space was added to the name or lower and higher case letters were mixed. For example, R is case sensitive and observes white space, while ``Homo sapiens'', ``Homo sapiens '' (with a trailing space) and ``homo sapiens'' are not recognized as different in typical spreadsheet software (e.g. Excel).
}

\author{ Jochen Fründ }

\seealso{\code{\link{empty}}}

\examples{
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"), 
lower = c("plant1","plant2","plant1","plant2","plant3","plant4"), 
webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
frame2webs(testdata,type.out="array")
sapply(frame2webs(testdata,type.out="list"), networklevel, index=c("connectance", "C score"))
}

\keyword{ package }


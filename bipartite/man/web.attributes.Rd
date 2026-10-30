\encoding{UTF-8}
\name{web_attributes}
\alias{higher_attributes}
\alias{higher_attributes<-}
\alias{lower_attributes}
\alias{lower_attributes<-}
\alias{meta_attributes}
\alias{meta_attributes<-}

\title{Access and replacement functions for the attributes of bipartite web and bipartite web array objects}
\description{
    With \code{higher_attributes} and \code{lower_attributes} you can 
    access / replace either the whole data.frames, specific columns or values 
    of the 
}

\usage{
higher_attributes(x, web_id = NULL)
higher_attributes(x, web_id = NULL) <- value

lower_attributes(x, web_id = NULL)
lower_attributes(x, web_id = NULL) <- value

meta_attributes(x, web_id = NULL)
meta_attributes(x, web_id = NULL) <- value
}

\arguments{
    \item{x}{bipartite_web or bipartite_webarray object}
    \item{web_id}{Only used when \option{x} is a bipartite_webarray object.}
    \item{value}{\itemize{
        \item{\code{data.frame}}
        \item{\code{list} of \code{data.frame}s}
        \item{\code{list}}
        \item{\code{list} of \code{list}s}
    }}
}

\value{
    \code{higher_attributes} and \code{lower_attributes} return a data.frame;
    \code{meta_attributes} returns a list.
}

\author{ Tobias Bauer \email{tobias.bauer-2@uni-hamburg.de} }

\seealso{See \link{bipartite_web} & \link{bipartite_webarray} for more information on how to create the S3-objects.}

\examples{
    # Access the higher_attribtues data.frame
    higher_attributes(Safariland)
    # Replacing the whole data.frame
    h_attr <- data.frame(id = 1:ncol(Safariland))
    rownames(h_attr) <- colnames(Safariland)
    higher_attributes(Safariland) <- h_attr
    # Adding a new element to the meta_attributes list
    meta_attributes(Safariland)$type <- "Pollination"
}

\keyword{ package }
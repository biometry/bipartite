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


\keyword{ package }
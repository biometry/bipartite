\encoding{UTF-8}
\name{bipartite_classes}
\alias{bipartite_web}
\alias{bipartite_webarray}

\title{Construct bipartite web and bipartite web array objects}
\description{
\code{bipartite_web()} creates an S3 object of class \code{"bipartite_web"} from a bipartite adjacency matrix.
\code{bipartite_webarray()} creates an S3 object of class \code{"bipartite_webarray"} from one or more bipartite webs,
combining them into a three-dimensional array.

The \code{"bipartite_web"} class represents a single bipartite network with node-level attributes,
while \code{"bipartite_webarray"} represents a collection of such webs.
}

\usage{
bipartite_web(web,
              higher_attributes = NULL,
              lower_attributes = NULL,
              meta_attributes = list())

bipartite_webarray(webs, ...,
                   global_higher_attributes = NULL,
                   global_lower_attributes = NULL,
                   local_higher_attributes = NULL,
                   local_lower_attributes = NULL,
                   global_meta_attributes = NULL,
                   local_meta_attributes = NULL)
}

\arguments{
    \item{web}{A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.}
    \item{higher_attributes}{A data.frame containing attributes of the higher trophic level species (columns) in the web. Note that the number of rows has to be equal to the number of columns in \option{x}. If \code{NULL}, an empty data.frame is created.}
    \item{lower_attributes}{A data.frame containing attributes of the lower trophic level species (rows) in the web. Note that the number of rows has to be equal to the number of rows in \option{x}. If \code{NULL}, an empty data.frame is created.}
    \item{meta_attributes}{A list of metadata associated with the web.}
    \item{webs}{Either a single matrix or a list of matrices.}
    \item{...}{further optional matrices provided, if x is not a list (see example).}
    \item{global_higher_attributes}{A \code{data.frame} of attributes for the higher-level nodes shared across all webs. Row names should match the column names of the array. If \code{NULL}, an empty \code{data.frame} is created automatically.}
    \item{global_lower_attributes}{A \code{data.frame} of attributes for the lower-level nodes shared across all webs. Row names should match the row names of the array. If \code{NULL}, an empty \code{data.frame} is created automatically.}
    \item{local_higher_attributes}{A named list of \code{data.frame}s containing higher-level node attributes specific to each web. Each element should have row names matching the column names of the array. If \code{NULL}, empty attribute frames are created for each web.}
    \item{local_lower_attributes}{A named list of \code{data.frame}s containing lower-level node attributes specific to each web. Each element should have row names matching the row names of the array. If \code{NULL}, empty attribute frames are created for each web.}
    \item{global_meta_attributes}{A list of metadata shared globally across the web array.
    If \code{NULL}, an empty list is created.}
    \item{local_meta_attributes}{A named list of metadata lists specific to each web.
    If \code{NULL}, empty lists are created for each web.}
}

\value{
\code{bipartite_web()} returns an object of class \code{"bipartite_web"}.

\code{bipartite_webarray()} returns an object of class \code{"bipartite_webarray"}.
}

\examples{
data(Safariland, vazquenc, vazquec)
Safariland.web <- bipartite_web(Safariland)
print(Safariland.web)
}

\keyword{ package }
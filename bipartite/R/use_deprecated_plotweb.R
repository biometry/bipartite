
# Function to switch the definition of plotweb to be either
# plotweb_v2 or plotweb_deprecated.
# The function changes the definition both locally in the "bipartite"
# namespace as well as globally. 
use_deprecated_plotweb <- function(deprecated = TRUE) {
    if (deprecated) {
        assignInNamespace("plotweb", plotweb_deprecated, ns = "bipartite")
        # Function exists in Global Namespace (package is loaded)
        if (exists("plotweb", where = 1)) {
            assign("plotweb", plotweb_deprecated, envir = .GlobalEnv)
        }
    } else {
        assignInNamespace("plotweb", plotweb_v2, ns = "bipartite")
        if (exists("plotweb", where = 1)) {
            assign("plotweb", plotweb_v2, envir = .GlobalEnv)
        }
    }
}

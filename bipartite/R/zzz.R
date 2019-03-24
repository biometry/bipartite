#.First.lib <- function(lib, pkg) {

# # library.dynam("bipartite", pkg, lib)

# vers <- paste(sessionInfo()$otherPkg$bipartite$Version,".",sep="")

.onLoad <- function(lib, pkg){
	library.dynam("bipartite", pkg, lib)	
}

.onAttach <- function(lib, pkg){	
     packageStartupMessage(" This is bipartite ",
                          utils::packageDescription("bipartite", field="Version"), ".\n For latest changes see versionlog in ?\"bipartite-package\". For citation see: citation(\"bipartite\").\n Have a nice time plotting and analysing two-mode networks.", appendLF=TRUE)}
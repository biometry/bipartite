#.First.lib <- function(lib, pkg) {

# # library.dynam("bipartite", pkg, lib)

# vers <- paste(sessionInfo()$otherPkg$bipartite$Version,".",sep="")

# cat(paste("----------------------------------------------------------\nThis is bipartite",vers,"\nFor latest additions type: ?bipartite.\nFor citation please type: citation(\"bipartite\").\nHave a nice time plotting and analysing two-mode networks.\n----------------------------------------------------------\n\n"))
# }

.onLoad <- function(lib, pkg){
	library.dynam("bipartite", pkg, lib)	
}

.onAttach <- function(lib, pkg)  {	
     packageStartupMessage(" This is bipartite ",
                          utils::packageDescription("bipartite", field="Version"), "\n For latest changes see versionlog in  ?\"bipartite-package\".\n For citation see: citation(\"bipartite\").\n Have a nice time plotting and analysing two-mode networks.\n")
}

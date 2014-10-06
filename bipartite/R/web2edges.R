web2edges <- function(web, webName=NULL, weight.column=TRUE, both.directions=FALSE, is.one.mode=FALSE, 
	out.files=c("edges", "names", "groups")[1:2], return=FALSE, verbose=FALSE){
	# function to turn a web-matrix into an edge list
	#
	# webName		name under which the files shall be saved
	# weight.column  returns a three-column file (two for the edges, third for the weights)
	# both.directions represents the edge list in both directions (e.g. 1 2 and 2 1)
	# out.files     vector of types of files to be written: edges, names, groups
	# verbose		shall a report be printed on console?
    #
	# authors: Carsten F. Dormann & Rouven Strauss
		
	edges <- which(web > 0, arr.ind=TRUE) 
	edge.weights <- web[edges]#, drop=FALSE]
		
	if (is.null(webName)) webName <- "web"	
	
	if (weight.column) {
		edges <- cbind(edges, edge.weights);
	}
	
	if (!is.one.mode) edges[,2] <- edges[,2] + nrow(web)
	# this leads to continuous numbering of species, starting with rows
	
	if (both.directions){
		sequ <- if (weight.column) c(2,1,3) else c(2,1)
		edges2 <- rbind(edges, edges[,sequ])
		edges <- edges2
	}
	
	if (return) {return(edges)} else {
		if ("edges" %in% out.files){
			write.table(edges, file=paste(webName, ".pairs", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
			if(verbose) print(paste("An edges file was created in the current directory: ", webName, ".pairs", sep=""))
		}
	
		if ("names" %in% out.files){
			NAMES <- cbind(virtual=1:sum(dim(web)), real=c(rownames(web), colnames(web)))
			write.table(NAMES, file=paste(webName, "-names.lut", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
			if(verbose) print(paste("A names file was created in the current directory: ", webName, "-names.lut", sep=""))
		}
	
		if ("groups" %in% out.files){
			GROUPS <- cbind(c(rownames(web), colnames(web)), rep(1:2, times=dim(web)))
			write.table(GROUPS, file=paste(webName, ".groups", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
			if(verbose) print(paste("A group file was created in the current directory: ", webName, ".groups", sep="")) 
		}
	}
}


npartite <- function(as.one.mode.web, index=c("connectance", "links per species", "mean degree")){
	# function to compute network statistics for n-partite networks, i.e. those without links within groups;
	# an attempt is made to ensure that these indices converge to the same value as returned for a bipartite network!
	# 
	# as.one.mode.web		a one-mode representation of a network; if this is based on several bipartite networks, i.e. species in levels NOT interacting within that level OR without interactions between separated levels (e.g. 1st and 3rd) then such forbidden links must be represented as NA! For bipartite networks this can easily be achieved using as.one.mode(., fill=NA), but obviously for bipartite networks this function is obsolete
	#
	aomw <- as.one.mode.web
	
	# check that it is a masked matrix; otherwise npartite makes no sense!
	if (attr(aomw, "one.mode") != "masked") stop("Input needs to be a 'masked' one-mode matrix. \nUse fill=NA in function as.one.mode to generate it.")
	
	results <- list()
	
	#--------
	if ("connectance" %in% index){
		results[["connectance"]] <- sum(aomw>0, na.rm=TRUE)/sum(!is.na(aomw))
	}
	
	#--------
	if ("links per species" %in% index | "mean degree" %in% index){
		results[["links per species"]] <- mean(colSums(aomw>0, na.rm=TRUE))/2
	}

	return(results)
	
	
}

#image(aomw <- as.one.mode(Safariland, fill=NA))
#npartite(aomw)
#networklevel(Safariland, index=c("connectance", "links per species"))

## needed: a function to turn a multi-level system into a masked matrix, where also interactions across more than one level are/not possible
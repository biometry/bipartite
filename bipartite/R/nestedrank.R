nestedrank <- function(web, method="NODF", weighted=TRUE, normalise=TRUE, return.matrix=FALSE){
	# following Alarcon et al. 2008
	# returns a matrix sorted by per-species contribution to nestedness, with most "nested" species first
	
	
	if (!any(c("NODF", "nodf", "binmatnest", "wine", "sort") %in% method)){
		 warning("Your choice for method has not been recognised. Will use 'NODF' instead.")
		 method <- "NODF"
	}
	
	if (any(dim(web) < 2)){
		warning("You web is too small for a meaningful computation of nestedrank (and probably other indices)!")
		out <- list("lower level"=rep(NA, NROW(web)), "higher level"=rep(NA, NCOL(web)))
	} else {
	
		# if the web has no names (e.g. null models), give them names:
		if (is.null(rownames(web))) rownames(web) <- paste0("L", seq.int(nrow(web)))
		if (is.null(colnames(web))) colnames(web) <- paste0("H", seq.int(ncol(web)))
			
		if (weighted == FALSE) web.for.wine <- (web>0)*1 else web.for.wine <- web
		
		nestedcomm <- switch(method, 
			NODF = nestednodf(web, weighted=weighted)$comm,
			nodf = nestednodf(web, weighted=weighted)$comm,
			binmatnest = nestedtemp(web)$comm,
			#{nn <- nestedness(web, null.models=FALSE); (web[nn$pack.order.row, nn$pack.order.col]>0)*1},
			wine = sortweb(wine(web.for.wine)$dij.w[nrow(web):1, ncol(web):1]),
			# a word of explanation: wine only returns a sorted-by-binary matrix; thus, we use the actual entries in this matrix to re-sort it, with the species with the highest sum of w_ij now being the most generalist
			sort = sortweb(web.for.wine)
		)	
		
		row.seq <- match(rownames(web), rownames(nestedcomm))
		names(row.seq) <- rownames(web)
		col.seq <- match(colnames(web), colnames(nestedcomm))
		names(col.seq) <- colnames(web)
		if (normalise){
				row.seq <- (row.seq-1)/(length(row.seq)-1)
				col.seq <- (col.seq-1)/(length(col.seq)-1)
		}
		out <- list("lower level"=row.seq, "higher level"=col.seq)
	}
	if (return.matrix) out$"nested.matrix" <- nestedcomm
	return(out)
}

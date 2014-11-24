linklevel <- function(web, index=c("dependence", "endpoint")){
	out <- list()
	
	if ("dependence" %in% index){
    	depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
	    depH <- web/matrix(colSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE)

		out$"HL dependence" <- depH
		out$"LL dependence" <- depL
	}
	
	if ("endpoint" %in% index)	out$endpoint <- endpoint(web)
	return(out)
}

#linklevel(Safariland)
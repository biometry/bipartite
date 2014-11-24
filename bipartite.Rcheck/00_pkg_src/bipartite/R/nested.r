nested <- function(web, method="binmatnest2", rescale=FALSE, normalised=TRUE){
  # a wrapper function to call any of the currently implemented measures of nestedness
 
  if (! any(method %in% c("binmatnest", "discrepancy", "binmatnest2", "discrepancy2", "NODF", "NODF2", "weighted NODF", "wine", "C.score", "checker", "ALL"))) stop("Typo? Unknown method!")
  if ("ALL" %in% method) index <- c("binmatnest", "discrepancy", "binmatnest2", "discrepancy2", "NODF", "NODF2", "weighted NODF", "wine", "C.score", "checker") else index <- method

  out <- NULL
	if ("binmatnest2" %in% index){ 
		nessy <- try(nestedtemp(web)$statistic, silent=TRUE)
		nessy.value <- if (inherits(nessy, "try-error")) NA else nessy
		out <- c(out, binmatnest2 = nessy.value)
	} 

	if ("binmatnest" %in% index){ # NA occur if web is full (i.e. no 0s)
		nessy <- try(nestedness(web, null.models = FALSE)$temperature, silent=TRUE)
		nessy.value <- if (inherits(nessy, "try-error")) NA else nessy
		out <- c(out, binmatnest = nessy.value)
	}
  
  if ("discrepancy2" %in% index) {
  	out <- c(out, "discrepancy2"=nesteddisc(web)$statistic)
  	}

  if ("discrepancy" %in% index) out <- c(out, "discrepancy"=unname(discrepancy(web)))
     
  if ("C.score" %in% index) out <- c(out, "C.score"=C.score(web, normalised=normalised))
  
  if ("checker" %in% index) out <- c(out, "checker"=nestedchecker(web)$C.score) 
  # identical to C.score(., FALSE)
 
  if ("NODF2" %in% index) out <- c(out, "NODF2"=unname(nestednodf(web, order=TRUE)$statistic[3])) # the nestedness of the row/column-sorted matrix (probably makes more sense)

  if ("NODF" %in% index) out <- c(out, "NODF"=unname(nestednodf(web, order=FALSE)$statistic[3])) # the "original", as I had implemented it, too (in NODF)

  if ("weighted NODF" %in% index) out <- c(out, "weighted NODF"=unname(nestednodf(web, order=TRUE, weighted=TRUE)$statistic[3]))
	
  if ("wine" %in% index) out <- c(out, "wine"=wine(web)$wine)
  
  if (rescale & ! "ALL" %in% method) warning("You requested rescaling, but you won't get it (unless you use method='ALL')!") 
  
  if (rescale & "ALL" %in% method) out <- abs(c(100,100,0,0,0,0,0,0,0,0) - out)
  
  out
  
}
# example:
#data(Safariland)
#nested(Safariland, "ALL")
#nested(Safariland, "ALL", rescale=TRUE)
#
#nested(Safariland, c("C.score", "checker"), rescale=FALSE)
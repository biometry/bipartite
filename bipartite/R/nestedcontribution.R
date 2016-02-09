nestedcontribution <- function(web, nsimul = 99){
	# following Saavedra et al. 2011
	# returns a list of species' contributions to nestedness
	
    # create a binary matrix
	web <- ifelse(web>0,1,0)

	# if the web has no names (e.g. null models), give them names:
	if (is.null(rownames(web))) rownames(web) <- paste0("L", seq.int(nrow(web)))
	if (is.null(colnames(web))) colnames(web) <- paste0("H", seq.int(ncol(web)))

	# build containers for the different contribution values
	lower.out <- data.frame(row.names=rownames(web))
	lower.out$nestedcontribution <- NA

	higher.out <- data.frame(row.names=colnames(web))
	higher.out$nestedcontribution <- NA

	if (any(dim(web) < 2)){
		warning("Your web is too small for a meaningful computation of nestedcontrrank (and probably other indices)!")
	}else{
		# compute the original nestedness
		nested.orig <- vegan::nestednodf(web)$statistic["NODF"]

		# compute the nestedness contribution of each row
		for(i in rownames(web)){
			message(i)
			# determine the probability of each interaction according to Bascompte et al. 2003
			probs <- (rowSums(web)[i] / ncol(web) + colSums(web) / nrow(web)) / 2.

			# calculate the nestedness for an ensemble of networks with row i's interactions drawn at random
			nested.null <- sapply(1:nsimul, function(x) {web.null <- web; web.null[i,] <- rbinom(ncol(web),1,probs); vegan::nestednodf(web.null)$statistic["NODF"]})

			# the below is the swap model
			#nested.null <- sapply(1:nsimul, function(x) {web.null <- web; web.null[i,] <- sample(web[i,]); vegan::nestednodf(web.null)$statistic["NODF"]})

			lower.out[i,"nestedcontribution"] <- (nested.orig - mean(nested.null))/sd(nested.null)
		}

		# compute the nestedness contribution of each column
		for(i in colnames(web)){
			message(i)
			# determine the probability of each interaction according to Bascompte et al. 2003
			probs <- (rowSums(web) / ncol(web) + colSums(web)[i] / nrow(web)) / 2.

			# calculate the nestedness for an ensemble of networks with row i's interactions drawn at random
			nested.null <- sapply(1:nsimul, function(x) {
				web.null <- web; 
				web.null[,i] <- rbinom(nrow(web),1,probs); 
				vegan::nestednodf(web.null)$statistic["NODF"]
				})
			
			# the below is the swap model
			#nested.null <- sapply(1:nsimul, function(x) {web.null <- web; web.null[,i] <- sample(web[,i]); vegan::nestednodf(web.null)$statistic["NODF"]})
			
			higher.out[i,"nestedcontribution"] <- (nested.orig - mean(nested.null))/sd(nested.null)
		}
	}

	out <- list("higher level"=higher.out, "lower level"=lower.out)	
	return(out)
}

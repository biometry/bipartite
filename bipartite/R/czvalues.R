czvalues <- function(moduleWebObject, weighted=FALSE, level="higher"){
	# function to compute c and z values of module members according to Guimera & Amaral (2005, Nature); formulae taken from Olesen et al. (2006, PNAS)
	# mod		an object produced by computeModules
	# weighted  logical; if TRUE computes c and z from quantitative (=weighted) data, based on strength, rather than degrees; 
	# level		"higher" or "lower" trophic level to compute c and z values for; defaults to "higher"
	#
	#	z = (k.is - ks.bar) / SD.ks  # within-module degree
	#  c = 1 - sum( (k.it/k.i)^2)    # among-module connectivity = participation coefficient P in Guimera & Amaral
	#
	# author: C.F. Dormann 19 Mar 2012
	#
	# note: these indices were developed for one-mode networks; we'll have to see whether they make sense for bipartite networks, too!
	#
	# k.is = number of links of i to other species in its own module s
	# ks.bar = average k.is of all species in module s
	# SD.ks = standard deviation of k.is of all species in module s
	# k.it = number of links of species i to module t
	# k.i = degree of species i
	
	# Note: for modules with only one species from the lower trophic level, the z-values will be NaN, since SD.ks is 0!
	# I decided to SET these values to 0, since they only occur when all species in that module will have the same number of links (which is obviously the case when there is only one lower-level species). Then the numerator is also 0. Thus, the value of 0 indicates that this species has no deviation from the rest of the module members (which is what I think z is supposed to represent).
	
	if(!isCorrectModuleWebObject(moduleWebObject)) stop("This function cannot be applied to this type of object!")
	
	## zvalues
	zvalues <- function(web, weighted=FALSE){
		# for the columns (typically higher trophic level)!!!
		# this function simply standardises the number of interactions within a module
		if (weighted){
			# computes "strength":
			depL <- web/matrix(rowSums(web), nrow = NROW(web), ncol = NCOL(web),  byrow = FALSE)
			k.is <- colSums(depL)
		} else {
			# computes "degrees":
			k.is <- colSums(web>0)
		}
		out <- (k.is - mean(k.is))/sd(k.is)
		# if there is only one species in a module:	
		#if (length(k.is) == 1) out <- nrow(web) 
		out
	}
	
	modInfo <- listModuleInformation(moduleWebObject)
	modules <- modInfo[[2]]
	nModules <- length(modInfo[[2]])
	#moduleList <- list()
	if (nModules < 2) stop("This web has no modules.") 
	#for (i in 1:nModules){
	#	moduleList[[i]] <- mod@originalWeb[modules[[i]][[1]], modules[[i]][[2]], drop=FALSE]
	#}

	web <- moduleWebObject@originalWeb
	
	if (level == "lower") web <- t(web) # simply transposes the web for the analysis
	
	if (weighted){
		# based on "strength":
		depL <- web/matrix(rowSums(web), nrow = NROW(web), ncol = NCOL(web),  byrow = FALSE)
		k.i <- colSums(depL)
	} else {
		#based on degrees:
		k.i <- colSums(web>0) # degrees of all species
	}
	# for number of links of species i in module t simply access only the subset of the full web containing these lower-trophic level species:
	z.values.all <- web[1,] # gets a vector automatically with species names in it
	k.it.for.t <- matrix(NA, ncol=NCOL(web), nrow=nModules) # contains the (k.it/k.i)^2-values for each species for each module
	for (t in 1:nModules){
		h.level.species.names <- if (level =="higher") modules[[t]][[2]] else modules[[t]][[1]]
		l.level.species.names <- if (level =="higher") modules[[t]][[1]] else modules[[t]][[2]]
		#k.it <- colSums(web[l.level.species.names, , drop=FALSE] > 0)
		if (weighted){
			# based on strength; compute strength for this module only (but for all species of the higher trophic level):
			depL.mod1 <- web[l.level.species.names,  , drop=FALSE]/matrix(rowSums(web[l.level.species.names,  , drop=FALSE]), nrow = NROW(web[l.level.species.names,  , drop=FALSE]), ncol = NCOL(web),  byrow = FALSE)# the part after the division (/) is "only" a matrix with the row totals in each column! It standardises the number of interactions per LTL species, so that they eventually (across all modules) add to 1.
			k.it <- colSums(depL.mod1)
		} else {
			 k.it <-  colSums(web[l.level.species.names, , drop=FALSE] > 0)
		}
		k.it.for.t[t,] <- (k.it/ k.i)^2
		# now also compute the z-values:
		z.values.all[h.level.species.names] <- zvalues(web[l.level.species.names, h.level.species.names , drop=FALSE], weighted=weighted)
	}
	c.values.all <- 1 - colSums(k.it.for.t)
	names(c.values.all) <- names(z.values.all)
	z.values.all[is.nan(z.values.all)] <- 0 ### THIS IS MY (CFD) decision. Otherwise these species had no z-value. I think 0 makes sense (see note).
	
	return(list("c"=c.values.all, "z"=z.values.all))
	
}

#moduleWebObject <- mod
#web <- moduleWebObject@originalWeb[modules[[1]][[1]], modules[[1]][[2]]]


# data(memmott1999)
# mod <- computeModules(memmott1999, steps=1E5)
# cz <- czvalues(mod)
# plot(cz[[1]], cz[[2]], pch=16, xlab="c", ylab="z", cex=0.8, xlim=c(0,1), las=1)
# abline(v=0.62) # threshold of Olesen et al. 2007
# abline(h=2.5)   # dito
# text(cz[[1]], cz[[2]], names(cz[[1]]), pos=4, cex=0.7)
#
# same for lower trophic level using also weights:
# cz <- czvalues(mod, level="lower", weighted=TRUE)
# plot(cz[[1]], cz[[2]], pch=16, xlab="c", ylab="z", cex=0.8, xlim=c(0,1), las=1)
# abline(v=0.62) # threshold of Olesen et al. 2007
# abline(h=2.5)   # dito
# text(cz[[1]], cz[[2]], names(cz[[1]]), pos=4, cex=0.7)


# multiCM <- function(web, N=100, steps=1000){
	# replicate(n=N, try(computeModules(web, ...)))
# }

#res <- multiCM(N=100, web=memmott1999, steps=100)
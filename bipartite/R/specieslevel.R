# `specieslevel` 
specieslevel <-
function(web, index="ALLBUTD", level="both", logbase=exp(1), low.abun=NULL, high.abun=NULL, PDI.normalise=TRUE, PSI.beta=c(1,0), nested.method="NODF", nested.normalised=TRUE, nested.weighted=TRUE, empty.web=TRUE) {
    # function to calculate bipartite web indices at the species level
    #
    # Carsten Dormann, Jochen Fruend & Denis Lippok, April/May 2007 - 2013


    # m <- matrix(c(4,7,0,0,9,1,2,0,5), 3, byrow=TRUE)
    if (!is.numeric(logbase)) stop("logbase must be numeric, e.g. 2 or exp(1) or 10!")

#! JFedit: we need an option empty.web just as in networklevel
  # must still be checked for case empty.web=FALSE (web.e not yet used anywhere)
    if(empty.web) {web <- empty(web)}   # delete unobserved species
    if(empty.web==FALSE) warning("empty.web=FALSE not yet supported in specieslevel")
    web.e <- empty(web) # emptied web for some indices 

    allindex <- c("degree", "normalised degree", "species strength", "nestedrank", "interaction push pull", "PDI", "resource range", "species specificity", "PSI", "NSI", "betweenness", "closeness", "Fisher alpha", "partner diversity", "effective partners", "d", "dependence", "proportional generality", "proportional similarity")
  #JFedit:
      # added "proportional similarity" (or PS), "proportional generality" (or "effective resource range" or "effective proportional resource use")
      # note that proportional generality can be higher than 1, (when a species selects for diversity)
#! JFedit: synonyms e.g. "paired differences index" only work if there is also an index name from the list above

    if ("ALL" %in% index) index <- allindex
    if ("ALLBUTD" %in% index) index <- allindex[-which(allindex=="dependence")]
    #out <- list()
    
    # only if indices are not given explicitly:
    if (length(index) == 1 & !all(index %in% allindex)){                        
    	index <- switch(index,
        	"ALL" = allindex,
        	"ALLBUTD" = allindex[-which(allindex=="dependence")],
         	stop("Your index is not recognised! Typo? Check help for options!", call.=FALSE) #default for all non-matches
          )
    }
    
    # catch cases where an invalid index is demanded:
    if (length(which(!(index %in% allindex))) > 0) warning(paste0("Index '", index[which(!(index %in% allindex))], "' not recognised and hence not computed!"))
    
    higher.out <- list()
    lower.out <- list()

    # sort out which level to compute indices for:
    if (level == "both") {for.higher <- TRUE; for.lower <- TRUE}
    if (level == "higher") {for.higher <- TRUE; for.lower <- FALSE}
    if (level == "lower") {for.higher <- FALSE; for.lower <- TRUE}
    if (!(level %in% c("both", "higher", "lower"))) stop("Please choose a valid level: 'both', 'higher' or 'lower'.")    
    
    
    # !!!!!! functions are split up over higher and lower level computations !!!!!!
    
    ############## helper functions (alphabetically) ########################
    
    BCC_weighted <- function(web, method="sum", level="higher", index="betweenness"){
      # computes weighted betweenness as proportion of shorted paths through this species
      # notice that these results are identical to the function "BC" when using binary data (web>)!
      # In weighted networks, some paths are actually LOST!
      #
      # Will return a warning when web is comparted!
      
      if (level == "higher") web <- t(web) 
      CO <- compart(web)
      
      if (CO$n.compart >= sum(web>0)){
        # projection fails when all species have their own compartment!
        out <- rep(NA, NCOL(web))
        return(out)
      }
          
      el <- web2edges(web, return=TRUE)
      # el <- symmetrise_w(el) # needed for undirected webs
      suppressWarnings(proj <- projecting_tm(el, method=method) )
      if (any(dim(proj) < 2)){
          warning("Web contains too few nodes to compute closeness or betweenness!")
          return(rep(NA, NROW(web))) # closeness/betweenness cannot be computed with only one link!
      } 
      
      if (index == "betweenness") {
        b <- betweenness_w(proj)[,2]
        if (length(b) != NROW(web)){
          b <- c(b, rep(NA, (NROW(web) - length(b))))
        }
        if (!is.null(rownames(web))) names(b) <- rownames(web)
        if (sum(b, na.rm=TRUE) == 0) out <- b else out <- b/sum(b, na.rm=TRUE)
      }
      
      if (index == "closeness") {
        ## closeness returns only the larger compartment, thus I first find the largest compartment and only compute the stuff for that;
        # the problem with closeness_w is that the names are lost!
        # find which is the largest compartment:
        which.is.where <- apply(CO[[1]], 1, function(x) sort(unique(x))[1] )			
        group <- names(sort(table(which.is.where), decreasing=TRUE))[1]
        keep <- which(CO[[1]] == group, arr.ind=TRUE)
        subweb <- web[unique(keep[,1]), unique(keep[,2]), drop=FALSE]
        # now compute closeness_w for that:
        el2 <- web2edges(subweb, return=TRUE)
        #if (any(dim(el2)) < 2) return(rep(NA, NROW(web))) # closeness cannot be computed with only one link!
        proj2 <- projecting_tm(el2, method=method)
        cc <- closeness_w(proj2, directed=TRUE, gconly=TRUE)[,2]
        # now pad these results with NAs:
        cc.full <- rep(NA, NROW(web))
        if (!is.null(rownames(web))){
          names(cc.full) <- rownames(web)
          cc.full[match(rownames(subweb), names(cc.full)) ] <- cc
          out <- cc.full
        } else {
          cc.full[1:length(cc)] <- cc
          out <- cc.full
        }
      }
      
      out
    }
    
    CoV <- function(web){
      	# variability of interactions, "species specificity index", following Poisot et al. (2012) using a normalisation to values between 0 and 1 (originally by Julliard et al. (2006)).
		numerator <- sqrt(colSums((web - matrix(colMeans(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE))^2))
		R <- NROW(web)
		denominator <- colSums(web) * sqrt((R-1)/R) 
		numerator/denominator
   	}

    PropSimilarity <- function(species_use, abuns){
    	# by Jochen Fruend
        p_i <- species_use / sum(species_use)
        q_i <- abuns / sum(abuns)
        sum(pmin(p_i, q_i))  
    }

   PSI <- function(web, beta=c(1,0)){
      # calculates the average contribution per visit for each pollinator species
      # (which in itself depends on the specialisation and abundance of the bees,
      # as well as the abundance of the plant species)
      #
      # web   a pollination web, with plants in rows
      # beta  a parameter accounting for the fact that two repeated landings are
      #       needed to transfer pollen: one for the source, one for the sink; a
      #       value of 2 would assume no memory of plant species in the pollinator;
      #  	a value of 1 implies infinitely storing pollen from source to sink plant
      #
      # developed by Dormann, Bluethgen & Gruber, 3 May 2007
      #
      # example:
      # m <- matrix(c(4,4,0,4,1,7), nrow=3, byrow=TRUE)
      # PSI(m, beta=1)
      
      Wi. <- matrix(rep(colSums(web), NROW(web)), nrow=NROW(web), byrow=TRUE)
      W.j <- matrix(rep(rowSums(web), NCOL(web)), ncol=NCOL(web), byrow=FALSE)
      
      PSImat <- (web/W.j)^beta * web/Wi.
      (PSI <- colSums(PSImat))
    }
   
   shannon <- function(x, base=exp(1)) {# shannon's diversity index
   		Pvec <- x/sum(x)
   		-sum(Pvec*log(Pvec, base=base), na.rm=TRUE)
   	} 
    
    ############### overall computations ####################

    if ("normalised degree" %in% index){
      nds <- ND(web)
    }
    
    if ("nestedrank" %in% index){
    	nested.ranks <- nestedrank(web, method=nested.method, normalise=nested.normalised, weighted=nested.weighted)
    }
    
    if (any(c("species strength", "dependence", "interaction push pull") %in% index)){
      depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
      depH <- web/matrix(colSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE)
      Dij <- depH-depL  # positive values indicate a stronger effect of i (=plants) on j (bees) than vice versa
    }
    
    if ("NSI" %in% index){
      NS <- nodespec(web)
    }
    
    if ("betweenness" %in% index){
      bcs <- BC(web)
    }
    
    if ("closeness" %in% index){
      ccs <- CC(web)
    }
    
    if (any(c("partner diversity", "effective partners", "proportional generality") %in% index)){
      preytot.mat <- matrix(rep(colSums(web), NROW(web)), NROW(web), byrow=TRUE)
      preyprop.mat <- web/preytot.mat  # = b_ik/b_.k in the first formula
      #H_Nk is the diversity index of inflow (diversity of flower visits for each pollinator)
      predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), NROW(web), byrow=FALSE)
      predprop.mat <- web/predtot.mat  # = b_kj/b_.k in the second formula
      H_Nk <- apply(preyprop.mat, 2, shannon, base=logbase)
      #H_Pk is the diversity index of pollinators for each plant species
      H_Pk <- apply(predprop.mat, 1, shannon, base=logbase)
      # next, we need the reciprocals of this
      # note that the ifelse is only needed if the web contains prey that is
      # not eaten or predators that don't eat ...
      n_Nk <- ifelse(colSums(web) != 0, logbase^H_Nk, 0)
      n_Pk <- ifelse(rowSums(web) != 0, logbase^H_Pk, 0)
    }
    
    
    
    ################################ higher level computations #########################
    if (for.higher){ 
      
    # species degrees:
    if ("degree" %in% index){
      sdH <- colSums(web>0)
      higher.out$"degree" <- unlist(sdH)
    }
    
    # normalised degrees:
    if ("normalised degree" %in% index){
      higher.out$"normalised degree" <- nds[[2]]
    }
       
    # dependence values, following the lead by Bascompte et al. 2006 (Science) and modifications suggested by Bluethgen et al. 2007 (Current Biology)
    if (any(c("species strength", "interaction push pull") %in% index)){ #"dependence", 
      #if ("dependence" %in% index){ # moved to linklevel!
      #  higher.out$"dependence" <- depH
      #}
      # species strength:
      if ("species strength" %in% index){ # strength = sum of dependences for a species (referenced in Bascompte et al. 2006)               
        SH <- colSums(depL) # accordingly ...
        higher.out$"species strength" <- SH
      } 
      # Interaction asymmetry (Vazquez et al. 2007, Oikos); rather similar to dependence above, really
      if ("interaction push pull" %in% index) {
        Aihigh <- colSums(-Dij)/colSums(web>0)
        higher.out$"interaction push pull" <- Aihigh
      } 
    } # end strength/dependence/interaction

    # nested ranks:
    if ("nestedrank" %in% index){
    	higher.out$"nestedrank" <- nested.ranks[["higher level"]]
    }

#! JFedit: this is PDI not PSI; before it returned PairedDifferenceIndex when asking for pollination support index...; name changed below as well
    # Poisot's paired differences index
      # comments should be adjusted here and below, though
    if (any(c("PDI", "paired differences index") %in% index)){
      higher.out$"PDI" <- PDI(web, normalise=PDI.normalise, log=FALSE)
    }
    
    # resource range
    if ("resource range" %in% index){
    	higher.out$"resource range" <- PDI(web>0)	
    }
    
    # species specificity (Juillard) = Coefficient of Variation (Poisot)
    if ("species specificity" %in% index){
    	higher.out$"species specificity index" <- CoV(web)
    }

    
    # Pollination webs only: Pollination service index
    if (any(c("PSI", "pollination service index") %in% index)){    
      if (for.lower & for.higher & length(PSI.beta) < 2) stop("You need to provide 2 values (in a vector) for PSI.beta.")
      higher.out$"PSI" <- PSI(web, beta=PSI.beta[1])
    }
    
    # node specialisation:
    if ("NSI" %in% index){
      higher.out$"node specialisation index NSI" <- NS$higher
    }
    
    # betweenness (incl. weighted):
    if ("betweenness" %in% index){
      higher.out$"betweenness" <- bcs[[2]]
      higher.out$"weighted betweenness" <- BCC_weighted(web, level="higher")
    }
    
    #closeness
    if ("closeness" %in% index){
      higher.out$"closeness" <- ccs[[2]]
      higher.out$"weighted closeness" <- BCC_weighted(web, level="higher", index="closeness")
    }
    
    # Fisher alpha:
    if ("Fisher alpha" %in% index){
      ff.high <- try(suppressWarnings(fisher.alpha(web, MARGIN=2)), silent=TRUE)
      higher.out$"Fisher alpha" <- if (!inherits(ff.high, "try-error")) ff.high else rep(NA, times=NCOL(web))
    }
    
    # diversity:
    if ("partner diversity" %in% index){
      higher.out$"partner diversity" <- H_Nk
    }
    
    # effective number of partners:
    if ("effective partners" %in% index){
      higher.out$"effective partners" <- n_Nk
    }
    
    # proportional generality
    # by Jochen Fruend March 2013
    if ("proportional generality" %in% index){
      mylow.abun <- if (is.null(low.abun)) mylow.abun <- rowSums(web) else mylow.abun <- low.abun
      pgenH <- n_Nk / logbase^(shannon(mylow.abun, base=logbase))      # uses the same abuns as dprime!
      higher.out$"proportional generality" <- pgenH
    }

    # proportional similarity (Jochen Fruend 2013)
    if ("proportional similarity" %in% index){
   	  # Feinsinger P, Spears EE, Poole RW: A simple measure of niche breadth. Ecology 1981, 61:27-32.
      mylow.abun <- if (is.null(low.abun)) mylow.abun <- rowSums(web) else mylow.abun <- low.abun
      psH <- apply(web, 2, PropSimilarity, abuns=mylow.abun) # uses the same abuns as dprime!
      higher.out$"proportional similarity" <- psH
    }
    
    # species-level standardised diversity index d
    if ("d" %in% index){
      dsH <- dfun(t(web), abuns=low.abun)[[1]]
      higher.out$d <- dsH
    }

    
    } # end condition "for.higher"
    
    
    
    ################################ lower level computations #########################
    if (for.lower){
    
    # species degrees:
    if ("degree" %in% index){
      sdL <- rowSums(web>0)
      lower.out$"degree" <- sdL
    }

    # normalised degrees:
    if ("normalised degree" %in% index){
      lower.out$"normalised degree" <- nds[[1]]
    }
    
    # dependence values, following the lead by Bascompte et al. 2006 (Science) and modifications suggested by Bluethgen et al. 2007 (Current Biology)
    if (any(c("species strength", "dependence", "interaction push pull") %in% index)){
      if ("dependence" %in% index){
        lower.out$"dependence" <- depL
      }
      # species strength:
      if ("species strength" %in% index){
        # strength = sum of dependences for a species (referenced in Bascompte et al. 2006)
        SL <- rowSums(depH) # a plant's strength is the sum of the dependencies of all its pollinators
        lower.out$"species strength" <- SL
      } 
      # Interaction asymmetry (Vazquez et al. 2007, Oikos); rather similar to dependence above, really
      if ("interaction push pull" %in% index) {
        Ailow <- rowSums(Dij)/rowSums(web>0)
        lower.out$"interaction push pull" <- Ailow
      }
    } # end strength/dependence/interaction

    # nested ranks:
    if ("nestedrank" %in% index){
    	lower.out$"nestedrank" <- nested.ranks[["lower level"]]
    }
    
    # Poisot's paired differences index
    if (any(c("PDI", "paired differences index") %in% index)){
      lower.out$"PDI" <- PDI(t(web), normalise=PDI.normalise, log=FALSE)
    }
    
    # species specificity (Juillard) = Coefficient of Variation (Poisot)
    if ("species specificity" %in% index){
    	lower.out$"species specificity index" <- CoV(t(web))
    }

    # resource range
    if ("resource range" %in% index){
    	lower.out$"resource range" <- PDI(t(web)>0)	
    }

    # Pollinator support index:
    if (any(c("PSI", "pollinator support index") %in% index)){    
      if (for.lower & for.higher & length(PSI.beta) < 2) stop("You need to provide 2 values (in a vector) for PSI.beta.")
      # need to accomodate a length-1-vector when only single level output is requested (not needed for PSI.higher):
      if (length(PSI.beta) == 1) PSI.beta.lower <- PSI.beta[1] else PSI.beta.lower <- PSI.beta[2]
      lower.out$"PSI"   <- PSI(t(web), beta=PSI.beta.lower)
    }
    
    # node specialisation:
    if ("NSI" %in% index){
      lower.out$"node specialisation index NSI" <- NS$lower
    }
    
    # betweenness (incl. weighted)
    if ("betweenness" %in% index){  
      lower.out$"betweenness" <- bcs[[1]]
      lower.out$"weighted betweenness" <- BCC_weighted(web, level="lower") 
    }
    
    # closeness:
    if ("closeness" %in% index){
      lower.out$"closeness" <- ccs[[1]]
      lower.out$"weighted closeness" <- BCC_weighted(web, level="lower", index="closeness") 
    }
    
    # Fisher alpha:
    if ("Fisher alpha" %in% index){
      ff.low <- try(suppressWarnings(fisher.alpha(web, MARGIN=1)), silent=TRUE)
      lower.out$"Fisher alpha" <- if (!inherits(ff.low, "try-error")) ff.low else rep(NA, times=NROW(web))
    }
    
    # diversity:
    if ("partner diversity" %in% index){
      lower.out$"partner diversity" <- H_Pk
    }
    
    # effective number of partners:
    if ("effective partners" %in% index){
      lower.out$"effective partners" <- n_Pk
    }

    #JFedit:
    # proportional similarity
    if ("proportional similarity" %in% index){
      myhigh.abun <- if (is.null(high.abun)) myhigh.abun <- colSums(web) else myhigh.abun <- high.abun
      psL <- apply(web, 1, PropSimilarity,abuns=myhigh.abun) # uses the same abuns as dprime!
      lower.out$"proportional similarity" <- psL
    }
    # Feinsinger P, Spears EE, Poole RW: A simple measure of niche breadth. Ecology 1981, 61:27-32.
    
    # proportional generality
    # by Jochen Fruend March 2013
    if ("proportional generality" %in% index){
      myhigh.abun <- if (is.null(high.abun)) myhigh.abun <- colSums(web) else myhigh.abun <- high.abun
      pgenL <- n_Pk / logbase^(shannon(myhigh.abun, base=logbase))  # uses the same abuns as dprime!
      lower.out$"proportional generality" <- pgenL
    }

    # species-level standardised diversity index d
    if ("d" %in% index){
      dsL <- dfun(web, abuns=high.abun)[[1]]
      lower.out$d <- dsL
    }
    
    
    } # end condition "for.lower"
    
    
    #---------------------------------------------------------------------------
    if (!("dependence" %in% index)) {
        higher.out <- as.data.frame(higher.out)
        lower.out <- as.data.frame(lower.out)
    }
    
    if (level == "lower") out <- lower.out
    if (level == "higher") out <- higher.out
    if (level == "both") out <- list("higher level"=higher.out, "lower level"=lower.out)  
    
    out
}


#specieslevel(bezerra2009, index="ALLBUTD", level="both", PSI.beta=c(1,1))
#specieslevel(bezerra2009, index=c("degree", "PSI"), level="lower")
#specieslevel(Safariland)
#specieslevel(Safariland, index="dependence", level="lower")
#specieslevel(Safariland, index=c("dependence", "d", "effective partners"), level="lower")


#JFedit:
 # below here some tests of the new functions, can be deleted most likely
#specieslevel(bezerra2009, level="higher", index=c("proportional similarity", "proportional generality"))

#plot(specieslevel(web=Safariland, level="higher", PSI.beta=c(1,1), index=c("proportional similarity", "proportional generality")))
#cor(specieslevel(Safariland, level="higher", PSI.beta=c(1,1),index=c("proportional similarity", "proportional generality", "d")))
#cor(specieslevel(Safariland, level="higher", PSI.beta=c(1,1),index=c("proportional similarity", "proportional generality", "d"))[colSums(Safariland)>10,])

#plot(specieslevel(Safariland, level="higher", PSI.beta=c(1,1),index=c("proportional similarity"))[,1])
#plot(specieslevel(Safariland, level="higher", PSI.beta=c(1,1),index=c("proportional similarity"))[,1], dfun(t(Safariland))$d)
#specieslevel(Safariland, level="both", PSI.beta=c(1,1),index=c("proportional similarity"))
#specieslevel(Safariland, level="both", PSI.beta=c(1,1),index=c("proportional generality"))

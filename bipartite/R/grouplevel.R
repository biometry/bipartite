grouplevel <- function(web, index="ALLBUTDD", level="both", weighted=TRUE, empty.web=FALSE, dist="horn", CCfun=mean, logbase="e", normalise=TRUE,  extinctmethod="r", nrep=100, fcdist="euclidean", fcweighted=TRUE){
    if (level == "both")   for.higher <- for.lower <- TRUE
    if (level == "higher") {for.higher <- TRUE; for.lower=FALSE}
    if (level == "lower")  {for.lower  <- TRUE; for.higher=FALSE}
    
     if ("vulnerability" %in% index & level == "higher") warning("You requested 'vulnerability' for the higher level, although it is not a higher level index! You will get 'generality' instead (same thing, really).", call.=FALSE)
      if ("generality" %in% index & level == "lower") warning("You requested 'generality' for the lower level, although it is not a lower level index! You will get 'vulnerability' instead (same thing, really).", call.=FALSE)
    
    
    if (for.higher) outh <- one.grouplevel(web, level="higher", index=index, weighted=weighted, empty.web=empty.web, dist=dist, CCfun=CCfun, logbase=logbase, normalise=normalise, extinctmethod=extinctmethod, nrep=nrep, fcdist=fcdist, fcweighted=fcweighted)

    if (for.lower) outl <- one.grouplevel(web, level="lower", index=index, weighted=weighted, empty.web=empty.web, dist=dist, CCfun=CCfun, logbase=logbase, normalise=normalise, extinctmethod=extinctmethod, nrep=nrep, fcdist=fcdist, fcweighted=fcweighted)
    
    if (level == "higher") return(outh)
    if (level == "lower")  return(outl)
    if (!("degree distribution" %in% index)) {
        out <- c(outh, outl)
        SEQ <- seq(1, length(out), by=2)
        out <- out[order(c(SEQ, SEQ+1))] # that needed some hard thinking, getting those bloody numbers align properly ...
    }
    
    if (level == "both"){
        if (!("degree distribution" %in% index)) {
            out <- c(outh, outl)
            SEQ <- seq(1, length(out), by=2)
            out <- out[order(c(SEQ, SEQ+1))] # that needed some hard thinking, getting those bloody numbers align properly ...
        } else {out <- list("HL"=outh, "LL"=outl)}
        return(out) # leave all the naming issues to function one.grouplevel
    }
}



one.grouplevel <- function(web, index="ALLBUTDD", level="higher", weighted=TRUE, empty.web=FALSE, dist="horn", CCfun=mean, logbase="e", normalise=TRUE, extinctmethod="r", nrep=100, fcdist="euclidean", fcweighted=TRUE){
  ### computes indices for one group level (i.e. higher/lower trophic level)
  
  web <- as.matrix(web)   

  ######
  if (level == "lower") web <- t(web)
  ######
  
  if(empty.web) {web <- empty(web)}
  web.e <- empty(web) # emptied web for some indices 
  if (NROW(web) < 2 | NCOL(web) <2) warning("Web is really too small to calculate any reasonable index. You will get the values nonetheless, but I wouldn't put any faith in them!")
                
   allindex <- c( "number of species", "mean number of links", "mean number of shared partners",
           "cluster coefficient", "degree distribution", "niche overlap",
           "togetherness", "C score", "V ratio", "discrepancy", "extinction slope", "robustness",    
           #quantitative series:
#! JFedit: effective partners removed here; do we really have fd and fc?; probably needs some work throughout
          "weighted cluster coefficient", "generality", "vulnerability", "partner diversity", "functional diversity", "functional complementarity")         
           # GONE: "mean interaction diversity", "effective partners" [syn. with generality/vulnerability, name still allowed],
  
       # only if indices are not given explicitly:
       if (length(index) == 1 & !all(index %in% allindex)){                        
       index <- switch(index,
               "ALL" = allindex,
               "ALLBUTDD" = allindex[-which(allindex=="degree distribution")],
               stop("Your index is not recognised! Typo? Check help for options!", call.=FALSE) #default for all non-matches
               )
           }
  
           
  # compute weights: sum of observed interactions per species (i.e., a list of 2)
  if (weighted){
    Wl <- apply(web, 1, sum)/sum(web)
    Wh <- apply(web, 2, sum)/sum(web)
  } else {
    Wl <- rep(1, NROW(web))
    Wh <- rep(1, NCOL(web))
  }
   
  out <- list()
    
  if (!(level %in% c("higher", "lower"))) stop("Please choose a valid level: 'higher' or 'lower'.")    
  
  ######################### computations for both levels ##########################
  
  # mean partner diversity
  if (any(c("partner diversity", "effective partners", "generality", "vulnerability") %in% index)){
      preytot.mat <- matrix(rep(colSums(web), NROW(web)), NROW(web), byrow=TRUE)
      preyprop.mat <- web/preytot.mat  # = b_ik/b_.k in the first formula
      #H_Nk is the diversity index of inflow (diversity of flower visits for each pollinator)
      predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), NROW(web), byrow=FALSE)
      predprop.mat <- web/predtot.mat  # = b_kj/b_.k in the second formula
      if (logbase==2 | logbase=="2"){
        H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x*log2(x), na.rm=TRUE))
        #H_Pk is the diversity index of pollinators for each plant species
        H_Pk <- apply(predprop.mat, 1, function(x) -sum(x*log2(x), na.rm=TRUE))
        # next, we need the reciprocals of this
        # note that the ifelse is only needed if the web contains prey that is
        # not eaten or predators that don't eat ...
        n_Nk <- ifelse(colSums(web)!=0, 2^H_Nk, 0)
        n_Pk <- ifelse(rowSums(web)!=0, 2^H_Pk, 0)
      }
      if (logbase=="e"){ # same code as above, just with "e"
        H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x*log(x), na.rm=TRUE))
        H_Pk <- apply(predprop.mat, 1, function(x) -sum(x*log(x), na.rm=TRUE))
        n_Nk <- ifelse(colSums(web)!= 0, exp(H_Nk), 0)
        n_Pk <- ifelse(rowSums(web)!= 0, exp(H_Pk), 0)
      }
  }  
   
   # degree distribution fits:
   if ("degree distribution" %in% index){
        dd <- suppressWarnings(try(degreedistr(web, plot.it=FALSE, pure.call=FALSE), silent=TRUE))
   }
             
  ###### for the higher level ######
  
    if ("number of species" %in% index) {
      spH <- c("number of species"=NCOL(web))
      out$"number of species" <- spH
    }
    
    ### 1. indices aggregated from individual species using mean or weighted.mean  #############
    
    # mean number of links
    if ("mean number of links" %in% index) {
      out$"mean number of links" <- weighted.mean(colSums(web>0), w=Wh)
    }

    if ("mean number of shared partners" %in% index) {
      # NO WEIGHTING HERE, since it works on distances between species.
      out$"mean number of shared partners"  <- mean(designdist(t(web)>0, method="J", terms="minimum"))     # HL
    }
    
    # cluster coefficient:
    if ("cluster coefficient" %in% index){
        Ci.high <- colSums((web>0))/nrow(web)
        out$"cluster coefficient" <- weighted.mean(Ci.high, w=Wh)
    }
    
    if ("weighted cluster coefficient" %in% index){  
        # compute the weighted cluster coefficient using tnet:
        edgelist <- web2edges(web, return=TRUE)
        wcc <- try(unname(clustering_tm(edgelist)["am"]), silent=TRUE) #uses arithmetic mean!
        out$"weighted cluster coefficient" <- if (inherits(wcc, "try-error")) "NA" else wcc
    }
    
    # mean niche overlap
    # mean similarity of niches (niche overlap, sensu Krebs, Ecological Methodology)
    # vegdist demands "sites" to be in rows, therefore the web has to be transposed
    # to calculate dissimilarity between higher level species; similarity is simply
    # 1-dissimilarity:
    if ("niche overlap" %in% index) { # no weighting!
      NOhigher <- mean(1-vegdist(t(web.e), method=dist))
      out$"niche overlap" <- NOhigher
    }  

  # togetherness (not symmetric)
  if ("togetherness" %in% index){
    out$"togetherness" <- togetherness(web, normalise=normalise, na.rm=TRUE)
  }
  
  # C.score (not symmetric)
  if ("C score" %in% index){
  	CS <- try(C.score(web, normalise=normalise, na.rm=TRUE), silent=TRUE)
    out$"C score" <- if (inherits(CS, "try-error")) NA else CS
  }
  
  # V.ratio (not symmetric)
  if ("V ratio" %in% index){
    out$"V ratio" <- V.ratio(web)
  }
  
  # discrepancy (not symmetric)
  if ("discrepancy" %in% index){
    out$discrepancy <- as.integer(unname(discrepancy(web)))
  }

  # degree distribution fits:
  if ("degree distribution" %in% index){
    if (class(dd)=="try-error"){
      print("Panic!")
      dd <- list()
      dd[[2]] <- NA
    }
    out$"degree distribution" <- dd[[2]] # HTL
  }
  # secondary extinction slope (not symmetric)
  if (any(c("extinction slope", "robustness") %in% index)){
    extL <- try(second.extinct(web=web, method=extinctmethod, nrep=nrep, participant="lower"), silent=TRUE)       
    if ("extinction slope" %in% index){
      slopeH <- try(slope.bipartite(extL, col="green", pch=16, type="b", plot.it=FALSE), silent=TRUE)
      out$"extinction slope"=suppressWarnings(as.numeric(slopeH))            
    }
    
    if ("robustness" %in% index) {
      if (inherits(extL, "try-error")){
        robustH <- NA
      } else {
        robustH <- try(robustness(extL), silent=TRUE)                
      }
      rH <- if (inherits(robustH, "try-error")) NA else robustH          
      out$"robustness" = as.numeric(rH)
    }
  }
  # Devoto's fd:
#! JFedit: fd changed to fc here (also fddist and fdweighted replaced throughout script); old: "functional diversity", 
  if (any(c("fc", "functional complementarity") %in% index)){   
    out$"functional complementarity" <- fc(t(web), dist=fcdist, method="average", weighted=fcweighted)
  }
      
    # diversity:
    if ("partner diversity" %in% index){
      out$"partner diversity" <- weighted.mean(H_Nk, w=Wh)
    }

#! JFedit: effective partners and generality/vulnerability now reduced to one index (before, they differed because only effective partners responded to "weighted")  
  # effective number of partners:
  if (any(c("generality", "vulnerability", "effective partners") %in% index)){
    G <- weighted.mean(n_Nk, w=Wh)
    if (level == "higher") out$generality <- G
    if (level == "lower") out$vulnerability <- G	
  }
      
      #  if (any(c("generality", "vulnerability") %in% index)){     
      #    # for weighted=T, this should be the same as G:
      #    # marginal totals-weighted mean exp(Shannon diversity)
      #    G <- sum(colSums(web)/sum(web)*n_Nk)
      #    
      #    if (level == "higher") out$generality <- G 
      #    if (level == "lower") out$vulnerability <- G
      #  }    
  
  
  ################ NOT included #################
  # PDI (no interpretation at group level)
  # PSI (no interpretation at group level)
  # NSI (no interpretation at group level)
  # betweenness, closeness (no interpretation at group level)
  # Fisher alpha (no interpretation at group level)

           
    #---------------------------------------------------------------------------
           
    if (level == "higher") names(out) <- paste(names(out), ".HL", sep="")
    if (level == "lower") names(out) <- paste(names(out), ".LL", sep="")
    
    if (!("degree distribution" %in% index)) {
      out <- unlist(as.data.frame(out))
      rownames(out) <- NULL
    }
  
    return(out)
           
}

# memmott takes ages!
#grouplevel(Safariland, level="both")
#grouplevel(Safariland, level="lower")
#grouplevel(Safariland, index="ALLBUTDD", level="lower")
#grouplevel(Safariland, index="degree distribution", level="lower")
#grouplevel(Safariland, index="degree distribution", level="both")
#grouplevel(Safariland, index="ALLBUTDD", level="higher", weighted=FALSE)
# outh <- one.grouplevel(Safariland, index="degree distribution", level="higher", weighted=FALSE)
#grouplevel(Safariland, index="extinction slope", level="lower")
#grouplevel(Safariland, index=c("degree distribution"))
#grouplevel(Safariland, index="ALL")
#grouplevel(Safariland, index="generality", level="higher")
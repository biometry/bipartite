`networklevel` <- function(web, index="ALLBUTDD", level="both", weighted=TRUE, ISAmethod="Bluethgen", SAmethod="Bluethgen", extinctmethod="r", nrep=100, CCfun=median, dist="horn", normalise=TRUE, empty.web=TRUE, logbase="e", intereven="prod", H2_integer=TRUE, fcweighted=TRUE, fcdist="euclidean", legacy=FALSE){
    ##
    ## web         interaction matrix, with lower trophic level in rows, higher in columns
    ## legacy      se to TRUE allows to run networklevel in its old form
#! JFedit: fddist and fdweighted replace throughout by fcdist and fcweighted
  
    if(empty.web) {web <- empty(web)}
    web.e <- empty(web) # emptied web for some indices 
    if (NROW(web) < 2 | NCOL(web) <2) warning("Web is really too small to calculate any reasonable index. You will get the values nonetheless, but I wouldn't put any faith in them!")
    
    allindex <- c( #descriptive:
        "number of species", "connectance", "web asymmetry", 
        #binary based and related:
        "links per species", "number of compartments", "compartment diversity",
        "cluster coefficient", "degree distribution", "mean number of shared partners",
        "togetherness", "C score", "V ratio", "discrepancy", "nestedness", 
        "weighted nestedness",
        #miscelleneous:
        "ISA", "SA", "extinction slope", "robustness", "niche overlap",   
        #quantitative series:
        "weighted cluster coefficient", "weighted NODF", "partner diversity", "generality", "vulnerability", "linkage density", "weighted connectance", "Fisher alpha",  "interaction evenness", "Alatalo interaction evenness", "Shannon diversity", "functional complementarity", "H2" )
    # GONE: "mean interaction diversity", "effective partners", 
    
    index <- unique(index) # enforces that each index name is used only once
    
    wrong.name <- which(is.na(pmatch(index, c(allindex, "ALL", "ALLBUTDD", "info","quantitative", "binary", "topology", "networklevel"))))
    if (length(wrong.name) > 0) stop("You selected an index that is not available: ", paste(index[wrong.name], collapse=", "))
    
    # only if indices are not given explicitly:
    if (length(index)==1 & !all(index %in% allindex)){
        index <- switch(index,
                        "ALL" = allindex,
                        "ALLBUTDD" = allindex[-which(allindex=="degree distribution")],
                        "info" = c("number of species", "connectance", "web asymmetry", "links per species", "number of compartments"),
                        # logic: only rough information on the network's general structure
                        "quantitative" = c("weighted cluster coefficient", "weighted nestedness", "weighted NODF", "functional complementarity", "partner diversity", "effective partners", "H2", "diversity","linkage density", "weighted connectance", "niche overlap"), #"mean interaction diversity", 
                        # logic: the "quantitative series"
                        "binary" = c("connectance", "links per species", "nestedness", "mean number of partners","cluster coefficient",  "C-score", "Fisher alpha"),
                        # logic: metrics for binary networks
                        "topology" = c("connectance", "cluster coefficient", "degree distribution", "togetherness", "nestedness"),
                        # logic: more abstract, topological metrics for binary networks
                        "networklevel" = c("connectance", "web asymmetry", "links per species", "number of compartments", "compartment diversity", "cluster coefficient", "nestedness", "weighted NODF", "ISA", "SA", "linkage density", "Fisher alpha", "diversity", "interaction evenness", "Alatalo interaction evenness", "H2"),
                        # only the truly networky indices
                        stop("Your index is not recognised! Typo? Check help for options!", call.=FALSE) #default for all non-matches
        )
    }
    
    if (legacy == FALSE){
        
        out <- list()
        # first: all the true networklevel indices
        #--------------------
        if ("connectance" %in% index){
            # connectance: "the fraction of all possible links that are realized in a network", (p. 12917 in Dunne et al. 2002)
            suppressWarnings(out$connectance <- sum(web>0)/prod(dim(web)))
        }
        #--------------------
        if ("web asymmetry" %in% index) out$"web asymmetry" <- (NCOL(web)-NROW(web))/sum(dim(web))     # web asymmetry (Bluethgen et al. 2007, Fig. S2)
        #--------------------
        if ("links per species" %in% index){
            L <- sum(web>0)/sum(dim(web))
            out$"links per species"=L
        }
        #--------------------
        if (any(c("number of compartments", "compartment diversity") %in% index)){
            CD <- function(co){
                if (co$n.compart>1){
                    no <- NA
                    for (i in 1:co$n.compart){
                        comp <- which(abs(co$cweb)==i, arr.ind=TRUE) # who is in a compartment?
                        no[i] <- length(unique(comp[,1])) + length(unique(comp[,2])) # how many species
                    }
                    no <- no/sum(dim(web)) # standardise for number of species in the web
                    CD <- exp(-sum(no*log(no)))
                }  else {CD <- NA} #; warning("only one compartment")}
                CD
            }
            
            comps <- try(compart(web.e), silent=TRUE)
            if (class(comps)=="try-error") {
                ncompart <- compdiv <- NA
            } else  {
                ncompart <- comps$n.compart
                compdiv <- CD(comps)
            }
            if ("number of compartments" %in% index) out$"number of compartments" <- as.integer(ncompart)
            if ("compartment diversity" %in% index) out$"compartment diversity" <- compdiv
        }
        #------------------------
        if ("cluster coefficient" %in% index){
            cluster.coef <- function(web, full=FALSE, FUN=mean){
                # calculate cluster coefficient
                # web   a bipartite web
                # full  logical; return cluster coefficients for each species in the network?
                # FUN   give a function to summarise individual cluster coefficients, defaults to 'mean'.
                # The concept was introduced by Watts & Strogatz (1998) and is described under in Wikipedia under http://en.wikipedia.org/w/index.php?title=Clustering_coefficient
                # Its main idea was to help identifying small-world networks, which should
                # have high cluster coefficients but low path lengths.
                # Ci of species i is simply the number of realised links devided by the number of possible links. At the network level, the CC for a network is the average CC of its members. This is a little bit fishy, since the CCs are log-normal
                # distributed (in pollinations networks at least). Therefore with FUN also
                # other summary measures can be employed.
                # Because within 'bipartite' we look at 2-mode networks, mean C is:
                # C_lowertrophiclevel = C_highertrophiclevel = C_entirenetwork.
                #
                # Literature: Watts DJ, Strogatz S (1998) Collective dynmics of 'small-world' networks. Nature 393:440-442
                
                # notice that this CC differs from the binary version of Tore Opsahl's. This here is the same as the one-mode cluster coefficient, Tore's is the "corrected" two-mode!
                web <- as.matrix(web)
                Ci.high <- colSums((web>0))/nrow(web)
                Ci.low <- rowSums((web>0))/ncol(web)
                CC <- FUN(Ci.high)
                if (full) out <- list("cluster coefficient"=CC, "CC values higher"=Ci.high,
                                      "CC values lower"=Ci.low) else out <- c("cluster coefficient"=CC)
                out
            }
            out$"cluster coefficient"=as.numeric(cluster.coef(web, FUN=CCfun, full=FALSE))
        }
        #-------------------
        if ("nestedness" %in% index){
            nest <- try(nestedtemp(web)$statistic, silent=TRUE)
            out$nestedness <- ifelse(class(nest)=="try-error", NA, nest)
            # a fast implementation of nestedness by Jari Oksanen
            #old: nestedness(web, null.models=FALSE)$temperature
        }
        #-------------------
        if ("weighted nestedness" %in% index){
            wine.res <- try(wine(web.e, nreps=nrep)$wine, silent=TRUE)
            out$"weighted nestedness" <- if (!inherits(wine.res, "try-error")) {wine.res} else {NA}
        }
        #-------------------
        if ("weighted NODF" %in% index){
			NODF <- try(unname(nestednodf(web, order=TRUE, weighted=TRUE)$statistic[3]), silent=TRUE)
            out$"weighted NODF" <- if (inherits(NODF, "try-error")) NA else NODF
        }
        #------------------
        if (any(c("ISA", "interaction strength asymmetry", "dependence asymmetry") %in% index)){
            # Dependence asymmetry (Bascompte et al. 2006; Bluethgen et al. 2007, Fig. S2)
            depL <- web.e/matrix(rowSums(web.e), nrow=NROW(web.e), ncol=NCOL(web.e), byrow=FALSE)
            depH <- web.e/matrix(colSums(web.e), nrow=NROW(web.e), ncol=NCOL(web.e), byrow=TRUE)
            
            if (ISAmethod=="Bascompte" & "ISA" %in% index) {
                #depMax <- depL
                #greaterindex <- depL < depH
                #depMax[greaterindex] <- depH[greaterindex]
                
                out$"dependence asymmetry"=mean(abs(depL-depH)/pmax(depL, depH), na.rm=TRUE)
            }
            if (ISAmethod=="Bluethgen" & "ISA" %in% index) {
                web2 <- web
                # delete cells for species encountered only once (Bluethgen, pers. comm.):
                web2[, which(colSums(web)==1)] <- 0
                web2[which(rowSums(web)==1), ] <- 0
                rowsummat <- matrix(rowSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=FALSE)
                colsummat <- matrix(colSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=TRUE)
                depL <- web2/rowsummat
                depH <- web2/colsummat
                
                depL[depL<=0] <- NA
                depH[depH<=0] <- NA
                # now we need a correction to account for the fact that links with few (e.g. 2) observations will have a minimum depL of 1/2: all on one species: depL=1, one on each of two: depL=0.5
                depLprime <- (depL - 1/rowsummat)/(1 - 1/rowsummat) 
                # assumes depLmin = 1/web2 and depLmax=1
                depHprime <- (depH - 1/colsummat)/(1 - 1/colsummat)
                
                out$"interaction strength asymmetry"=mean(as.matrix(depHprime-depLprime), na.rm=TRUE) #ranges from -1 to 1 /sum(depLprime, depHprime, na.rm=TRUE)
            }
        }
        #------------------
        if ("SA" %in% index){
            # Specialisation asymmetry (Bluethgen et al. 2007, Fig. S2)
            # 2 options for calculating the "mean" SA:
            # either as Bluethgen et al: average weighted by number of interactions in the 
            # cell or as mean of logarithms (since the dependencies follow a lognormal
            # distribution)
            di <- dfun(web)$dprime  # plants
            dj <- dfun(t(web))$dprime # pollinators
            if (SAmethod=="log"){
                lgmeani <- mean(log(di[di>0])); lgmeanj <- mean(log(dj[dj>0]))
                SA <- (lgmeanj-lgmeani)/sum(lgmeani, lgmeanj)  
                # ij-sequence changed because log changes sequence, too
            }
            if (SAmethod=="Bluethgen"){
                wmeani <- sum(di*rowSums(web.e))/sum(web.e)
                wmeanj <- sum(dj*colSums(web.e))/sum(web.e)
                SA <- (wmeanj-wmeani)/sum(wmeani, wmeanj) 
                # positive values indicate more specialisation in the higher trophic level
            }
            out$"specialisation asymmetry" <- SA 
        }
        #------------------
        if (any(c("linkage density","weighted connectance") %in% index)){
            # for formula see Tylianakis et al. (2006), supplement.
            # N refers to prey, P to predators
            
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
            # mean number of predators per prey
            V <- sum(rowSums(web)/sum(web)*n_Pk)
            # marginal totals-weighted mean exp(Shannon diversity)
            G <- sum(colSums(web)/sum(web)*n_Nk)
            # linkage density
            LD_q <- 0.5*(V+G)
            #------------------
            if ("linkage density" %in% index) out$"linkage density" <- LD_q
            if ("weighted connectance" %in% index) out$"weighted connectance" <- LD_q/sum(dim(web))
            #LD_qs <- LD_q/(NROW(web)+NCOL(web)) # "weighted food web connectance", according to Jason's appendix
            #------------------
            # #We found no reference to this metric and saw little use for it. It is very similar to vulnerability/generality and can easily be computed from the output of \code{\link{specieslevel}} as \code{mean(specieslevel(web, index="diversity"))}        
            # if ("mean interaction diversity" %in% index){
            # out$"interaction diversity LTL" <- mean(H_Nk)
            # out$"interaction diversity HTL" <- mean(H_Pk)
            #}
        }
        #------------------
        if ("Fisher alpha" %in% index) {
            fish <- try(fisherfit(web)$estimate, silent=TRUE) #in vegan        
            if (inherits(fish, "try-error")) {
                out$"Fisher alpha" <- NA
            } else {
                out$"Fisher alpha" <- fish
            }
        }
        #----------------------------- evenness & diversity ----------------------
        if (any(c("interaction evenness", "Alatalo interaction evenness", "Shannon diversity") %in% index)){
            # interaction evenness
            p_i.mat <- web/sum(web)
            #---------------        
            SH <- -sum(p_i.mat*log(p_i.mat), na.rm=TRUE)
            if ("Shannon diversity" %in% index) out$"Shannon diversity" <- SH
            IE <- ifelse(intereven=="prod", SH/log(prod(dim(web))), SH/log(sum(web>0)))
            #---------------
            if ("interaction evenness" %in% index) out$"interaction evenness" <- IE
            #---------------
            if ("Alatalo interaction evenness" %in% index){
            evenness <- function(web){
                # calculates evenness of the numbers of individual of different species in
                # a community, NOT according to formula in Mueller et al. (1999, 
                # J. Anim. Ecol), but according to the original formula in Alatalo 
                # (1981, Oikos) 
                # can be extended at some point to more indices ...
                pk <- web/sum(web)
                (Alatalo <- (1/sum(pk^2) -1) / (exp(-sum(pk * log(pk), na.rm=TRUE)) -1))
            }
            E <- evenness(web)
            out$"Alatalo interaction evenness" <- E
            }
        }      
        #---------------
        # Bluethgen's H2'
        if ("H2" %in% index){
            is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # from the help of is.integer!
            if (any(is.wholenumber(web)==FALSE)) H2_integer <- FALSE # turns H2_integer off if values are not integers
            H2 <- as.numeric(H2fun(web, H2_integer=H2_integer)[1]) #1.element is the standardised H2 prime
            out$"H2"= ifelse(H2<0, 0, H2)
        }
        #----------------------- now: grouplevel -------------------
        # a list of network indices (which should not be called through grouplevel):
#! JFedit: weighted connectance added here
        netw.index <- match(c("connectance", "web asymmetry", "links per species", "number of compartments", "compartment diversity", "nestedness", "weighted nestedness", "weighted NODF", "ISA", "SA", "interaction evenness", "Alatalo interaction evenness", "Fisher alpha", "H2", "Shannon diversity", "linkage density", "weighted connectance"), index)
        exclude.index <- netw.index[!is.na(netw.index)]
        gindex <- if (length(exclude.index)==0) index else index[-exclude.index] # exclude NAs from this vector
        if (length(gindex) > 0) outg <- grouplevel(web, index=gindex, level=level, weighted=weighted, extinctmethod=extinctmethod, nrep=nrep, CCfun=CCfun, dist=dist, normalise=normalise, empty.web=empty.web, logbase=logbase, fcweighted=fcweighted, fcdist=fcdist)
        if (exists("outg")){
            # coerce potential list of HL/LL to one list:
            if (is.list(outg)){
                SEQ <- seq(1, 2*length(outg[[1]]), by=2)
                sorted.outg <- c(outg[[1]], outg[[2]])
                outg <- sorted.outg[order(c(SEQ, SEQ+1))]
            } 
            
            out <- c(unlist(out), outg) 
        } else {out <- unlist(out)} 
    } # end legacy = FALSE      
    
    if (legacy == TRUE){ # for the old way of running networklevel
        out <- .networklevel(web, index=index, ISAmethod=ISAmethod, SAmethod=SAmethod, extinctmethod=extinctmethod, nrep=nrep, plot.it.extinction=FALSE, plot.it.dd=FALSE, CCfun=CCfun, dist=dist, normalise=normalise, empty.web=empty.web, logbase=logbase, fcweighted=fcweighted, fcdist=fcdist)
        
    } # end legacy condition
    
    return(out)
}


#! JFedit: generality / vs. effective partners etc. NOT changed here! should it be? the changes are not really "legacy"?
.networklevel <- function(web, index="ALLBUTDD", ISAmethod="Bluethgen", SAmethod="Bluethgen", extinctmethod="r", nrep=100, plot.it.extinction=FALSE, plot.it.dd=FALSE, CCfun=median, dist="horn", normalise=TRUE, empty.web=TRUE, logbase="e", intereven="sum", H2_integer=TRUE, fcweighted=TRUE, fcdist="euclidean"){
    
    ##
    ## web         interaction matrix, with lower trophic level in rows, higher in columns
    ## legacy      se to TRUE allows to run networklevel in its old form
    
    if(empty.web) {web <- empty(web)}
    web.e <- empty(web) # emptied web for some indices 
    if (NROW(web) < 2 | NCOL(web) <2) warning("Web is really too small to calculate any reasonable index. You will get the values nonetheless, but I wouldn't put any faith in them!")
    
    allindex <- c( #descriptive:
        "number of species", "connectance", "web asymmetry", 
        #binary based and related:
        "links per species", "number of compartments", "compartment diversity",
        "cluster coefficient", "degree distribution", "mean number of shared partners",
        "togetherness", "C score", "V ratio", "discrepancy", "nestedness", 
        "weighted nestedness",
        #miscelleneous:
        "ISA", "SA", "extinction slope", "robustness", "niche overlap",   
        #quantitative series:
        "weighted cluster coefficient", "weighted NODF", "generality", "vulnerability", "linkage density", "Fisher alpha",  "interaction evenness", "Alatalo interaction evenness", "diversity", "effective partners", "functional complementarity", "H2" )
    
    # GONE: "mean interaction diversity", 
    
    # only if indices are not given explicitly:
    if (length(index)==1 & !all(index %in% allindex)){                        
        index <- switch(index,
                        "ALL" = allindex,
                        "ALLBUTDD" = allindex[-which(allindex=="degree distribution")],
                        "info" = c("number of species", "connectance", "web asymmetry", "links per species", "number of compartments"),
                        # logic: only rough information on the network's general structure
                        "quantitative" = c("weighted cluster coefficient", "weighted nestedness", "weighted NODF", "functional complementarity", "H2", "diversity", "effective partners", "mean interaction diversity", "linkage density"),
                        # logic: the "quantitative series"
                        "binary" = c("connectance", "links per species", "nestedness", "cluster coefficient",  "C-score"),
                        # logic: metrics for binary networks
                        "topology" = c("connectance", "cluster coefficient", "degree distribution", "togetherness", "nestedness"),
                        # logic: more abstract, topological metrics for binary networks
                        stop("Your index is not recognised! Typo? Check help for options!", call.=FALSE) #default for all non-matches
        )
    }
    out <- list()

# set up enough panels for plotting:
if (plot.it.extinction) {m=2; n=1; par(mfrow=c(m,n), mar=c(5,5,4,1))} else m=1
if (plot.it.dd) {n=2; par(mfrow=c(m,n), mar=c(5,5,4,1))} else n=1

#-------------------
if ("number of species" %in% index) {
    out$"number of lower trophic species" <- as.integer(NROW(web))
    out$"number of higher trophic species" <- as.integer(NCOL(web))
}
#--------------------
if ("connectance" %in% index){
    # connectance: "the fraction of all possible links that are realized in a network", (p. 12917 in Dunne et al. 2002)
    out$connectance <- sum(web>0)/prod(dim(web))
}
#--------------------
if ("web asymmetry" %in% index) out$"web asymmetry" <- (NCOL(web)-NROW(web))/sum(dim(web))     # web asymmetry (Bluethgen et al. 2007, Fig. S2)
###---###---###---###---###---###---###---###---###---###---###---###---###
#--------------------
if ("links per species" %in% index){
    L <- sum(web>0)/sum(dim(web))
    out$"links per species"=L
}
#--------------------
if (any(c("number of compartments", "compartment diversity") %in% index)){
    CD <- function(co){
        if (co$n.compart>1){
            no <- NA
            for (i in 1:co$n.compart){
                comp <- which(abs(co$cweb)==i, arr.ind=TRUE) # who is in a compartment?
                no[i] <- length(unique(comp[,1])) + length(unique(comp[,2])) # how many species
            }
            no <- no/sum(dim(web)) # standardise for number of species in the web
            CD <- exp(-sum(no*log(no)))
        }  else {CD <- NA} #; warning("only one compartment")}
        CD
    }
    
    comps <- try(compart(web.e), silent=TRUE)
    if (class(comps)=="try-error") {
        ncompart <- compdiv <- NA
    } else  {
        ncompart <- comps$n.compart
        compdiv <- CD(comps)
    }
    if ("number of compartments" %in% index) out$"number of compartments" <- as.integer(ncompart)
    if ("compartment diversity" %in% index) out$"compartment diversity" <- compdiv
}
#------------------------
if ("cluster coefficient" %in% index){
    cluster.coef <- function(web, full=FALSE, FUN=mean){
        # calculate cluster coefficient
        # web   a bipartite web
        # full  logical; return cluster coefficients for each species in the network?
        # FUN   give a function to summarise individual cluster coefficients, defaults to 'mean'.
        # The concept was introduced by Watts & Strogatz (1998) and is described under in Wikipedia under http://en.wikipedia.org/w/index.php?title=Clustering_coefficient
        # Its main idea was to help identifying small-world networks, which should
        # have high cluster coefficients but low path lengths.
        # Ci of species i is simply the number of realised links devided by the number of possible links. At the network level, the CC for a network is the average CC of its members. This is a little bit fishy, since the CCs are log-normal
        # distributed (in pollinations networks at least). Therefore with FUN also
        # other summary measures can be employed.
        # Because within 'bipartite' we look at 2-mode networks, mean C is:
        # C_lowertrophiclevel = C_highertrophiclevel = C_entirenetwork.
        #
        # Literature: Watts DJ, Strogatz S (1998) Collective dynmics of 'small-world' networks. Nature 393:440-442
        
        # notice that this CC differs from the binary version of Tore Opsahl's. This here is the same as the one-mode cluster coefficient, Tore's is the "corrected" two-mode!
        web <- as.matrix(web)
        Ci.high <- colSums((web>0))/nrow(web)
        Ci.low <- rowSums((web>0))/ncol(web)
        CC <- FUN(Ci.high)
        if (full) out <- list("cluster coefficient"=CC, "CC values higher"=Ci.high,
                              "CC values lower"=Ci.low) else out <- c("cluster coefficient"=CC)
        out
    }
    out$"cluster coefficient"=as.numeric(cluster.coef(web, FUN=CCfun, full=FALSE))
}

if ("weighted cluster coefficient" %in% index){  
    # compute the weighted cluster coefficient using tnet:
    edgelist <- web2edges(web, return=TRUE)
    wcc <- try(unname(clustering_tm(edgelist)["am"]), silent=TRUE) #uses arithmetic mean!
    out$"weighted cluster coefficient" <- if (inherits(wcc, "try-error")) "NA" else wcc
}

#--------------------------
# degree distribution fits:
if ("degree distribution" %in% index){
    dd <- suppressWarnings(try(degreedistr(web, plot.it=plot.it.dd, pure.call=FALSE), silent=TRUE))
    if (class(dd)=="try-error"){
        dd <- list()
        dd$"dd fits LTL" <- NA
        dd$"dd fits HTL" <- NA
    }
    out$"degree distribution LTL" <- dd$"dd fits LTL"
    out$"degree distribution HTL" <- dd$"dd fits HTL"
}
#-------------------
# The Stone & Roberts's indices: S, T and C score:
#-------------------
if ("mean number of shared partners" %in% index) {
    out$"mean number of shared HTL species" <- 
        mean(designdist(t(web)>0, method="J", terms="minimum"))
    out$"mean number of shared LTL species" <- 
        mean(designdist(web>0, method="J", terms="minimum"))
} 
#-------------------
if ("togetherness" %in% index){
    out$togetherness <- togetherness(web, normalise=normalise, na.rm=TRUE)
}
#-------------------
if ("C score" %in% index){
    out$"C score" <- C.score(web, normalise=normalise, na.rm=TRUE)
}

#-------------------
if ("V ratio" %in% index){
    out$"V ratio" <- V.ratio(web)
}
#-------------------
if ("discrepancy" %in% index){
    out$discrepancy <- as.integer(unname(discrepancy(web)))
}
#-------------------
if ("nestedness" %in% index){
    nest <- try(nestedtemp(web)$statistic, silent=TRUE)
    out$nestedness <- ifelse(class(nest)=="try-error", NA, nest)
    # a fast implementation of nestedness by Jari Oksanen
    #old: nestedness(web, null.models=FALSE)$temperature
}
#-------------------
if ("weighted nestedness" %in% index){
    wine.res <- try(wine(web.e, nreps=nrep)$wine, silent=TRUE)
    out$"weighted nestedness" <- if (!inherits(wine.res, "try-error")) {wine.res} else {NA}
}
#-------------------
if ("weighted NODF" %in% index){
    out$"weighted NODF" <- unname(nestednodf(web, order=TRUE, weighted=TRUE)$statistic[3])
}

###---###---###---###---###---###---###---###---###---###---###---###---###
if ("ISA" %in% index){#(any(c("ISA", "interaction strength asymmetry", "dependence asymmetry")) %in% index){
    # Dependence asymmetry (Bascompte et al. 2006; Bluethgen et al. 2007, Fig. S2)
    depL <- web.e/matrix(rowSums(web.e), nrow=NROW(web.e), ncol=NCOL(web.e), byrow=FALSE)
    depH <- web.e/matrix(colSums(web.e), nrow=NROW(web.e), ncol=NCOL(web.e), byrow=TRUE)
    
    if (ISAmethod=="Bascompte" & "ISA" %in% index) {
        #depMax <- depL
        #greaterindex <- depL < depH
        #depMax[greaterindex] <- depH[greaterindex]
        
        out$"dependence asymmetry"=mean(abs(depL-depH)/pmax(depL, depH), na.rm=TRUE)
    }
    if (ISAmethod=="Bluethgen" & "ISA" %in% index) {
        web2 <- web
        # delete cells for species encountered only once (Bluethgen, pers. comm.):
        web2[, which(colSums(web)==1)] <- 0
        web2[which(rowSums(web)==1), ] <- 0
        rowsummat <- matrix(rowSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=FALSE)
        colsummat <- matrix(colSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=TRUE)
        depL <- web2/rowsummat
        depH <- web2/colsummat
        
        depL[depL<=0] <- NA
        depH[depH<=0] <- NA
        # now we need a correction to account for the fact that links with few (e.g. 2) observations will have a minimum depL of 1/2: all on one species: depL=1, one on each of two: depL=0.5
        depLprime <- (depL - 1/rowsummat)/(1 - 1/rowsummat) 
        # assumes depLmin = 1/web2 and depLmax=1
        depHprime <- (depH - 1/colsummat)/(1 - 1/colsummat)
        
        out$"interaction strength asymmetry"=mean(as.matrix(depHprime-depLprime), na.rm=TRUE) #ranges from -1 to 1 /sum(depLprime, depHprime, na.rm=TRUE)
    }
}

#---------------------------------------------------------
# Specialisation asymmetry (Bluethgen et al. 2007, Fig. S2)
# 2 options for calculating the "mean" SA:
# either as Bluethgen et al: average weighted by number of interactions in the 
# cell or as mean of logarithms (since the dependencies follow a lognormal
# distribution)
if ("SA" %in% index){#(any(c("SA", "specialisation asymmetry")) %in% index){
    di <- dfun(web)$dprime  # plants
    dj <- dfun(t(web))$dprime # pollinators
    if (SAmethod=="log"){
        lgmeani <- mean(log(di[di>0])); lgmeanj <- mean(log(dj[dj>0]))
        SA <- (lgmeanj-lgmeani)/sum(lgmeani, lgmeanj)  
        # ij-sequence changed because log changes sequence, too
    }
    if (SAmethod=="Bluethgen"){
        wmeani <- sum(di*rowSums(web.e))/sum(web.e)
        wmeanj <- sum(dj*colSums(web.e))/sum(web.e)
        SA <- (wmeanj-wmeani)/sum(wmeani, wmeanj) 
        # positive values indicate more specialisation in the higher trophic level
    }
    out$"specialisation asymmetry" <- SA 
}

#--------------------------
# species extinction curve:
if (any(c("extinction slope", "robustness") %in% index)){
    extL <- try(second.extinct(web=web, method=extinctmethod, nrep=nrep, participant="lower"), silent=TRUE)       
    extH <- try(second.extinct(web=web, method=extinctmethod, nrep=nrep, participant="higher"), silent=TRUE)
    if ("extinction slope" %in% index){
        slopeL <- try(slope.bipartite(extL, col="green", pch=16, type="b", plot.it=plot.it.extinction), silent=TRUE)
        out$"extinction slope LTL"=suppressWarnings(as.numeric(slopeL))            
        
        slopeH <- try(slope.bipartite(extH, col="green", pch=16, type="b", plot.it=plot.it.extinction), silent=TRUE)
        out$"extinction slope HTL"=suppressWarnings(as.numeric(slopeH))
    }
    
    if ("robustness" %in% index) {
        if (inherits(extL, "try-error")){
            robustL <- robustH <- NA
        } else {
            robustL <- try(robustness(extL), silent=TRUE)
            robustH <- try(robustness(extH), silent=TRUE)
            
            rL <- if (inherits(robustL, "try-error")) NA else robustL		  
            rH <- if (inherits(robustH, "try-error")) NA else robustH		  
        }
        out$"robustness higher exterminated" = as.numeric(rH)
        out$"robustness lower exterminated" = as.numeric(rL)
    }
}
#--------------------------
# mean similarity of niches (niche overlap, sensu Krebs, Ecological Methodology)
# vegdist demands "sites" to be in rows, therefore the web has to be transposed
# to calculate dissimilarity between higher level species; similarity is simply
# 1-dissimilarity:
if ("niche overlap" %in% index) {
    NOhigher <- mean(1-vegdist(t(web.e), method=dist))
    NOlower <- mean(1-vegdist(web.e, method=dist))
    out$"niche overlap LTL" <- NOlower
    out$"niche overlap HTL" <- NOhigher
}
###---###---###---###---###---###---###---###---###---###---###---###---###
if (any(c("links per species", "linkage density", "vulnerability", "generality", "Fisher alpha", "mean interaction diversity") %in% index)){
    # for formula see Tylianakis et al. (2006), supplement.
    # N refers to prey, P to predators
    
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
    # mean number of predators per prey
    V <- sum(rowSums(web)/sum(web)*n_Pk)
    # marginal totals-weighted mean exp(Shannon diversity)
    G <- sum(colSums(web)/sum(web)*n_Nk)
    # linkage density
    LD_q <- 0.5*(V+G)
    if ("vulnerability" %in% index) out$"vulnerability" <- V
    #------------------
    if ("generality" %in% index) out$"generality"=G
    #------------------
    if ("linkage density" %in% index) out$"linkage density"=LD_q
    #LD_qs <- LD_q/(NROW(web)+NCOL(web)) # "weighted food web connectance", according to Jason's appendix
    #------------------
    if ("Fisher alpha" %in% index) {
        fish <- try(fisherfit(web)$estimate, silent=TRUE) #in vegan        
        if (inherits(fish, "try-error")) {
            out$"Fisher alpha" <- NA
        } else {
            out$"Fisher alpha" <- fish
        }
    }
    #------------------
    # #We found no reference to this metric and saw little use for it. It is very similar to vulnerability/generality and can easily be computed from the output of \code{\link{specieslevel}} as \code{mean(specieslevel(web, index="diversity"))}        
    # if ("mean interaction diversity" %in% index){
    # out$"interaction diversity LTL" <- mean(H_Nk)
    # out$"interaction diversity HTL" <- mean(H_Pk)
    #}
}
#----------------------------- evenness & diversity ----------------------
if (any(c("interaction evenness", "Alatalo interaction evenness", "Shannon diversity") %in% index)){
    # interaction evenness
    p_i.mat <- web/sum(web)
    SH <- -sum(p_i.mat*log(p_i.mat), na.rm=TRUE)
    IE <- ifelse(intereven=="prod", SH/log(prod(dim(web))), SH/log(sum(web>0)))
    #---------------
    if ("interaction evenness" %in% index) out$"interaction evenness" <- IE
    #---------------
    if ("Alatalo interaction evenness" %in% index){
        evenness <- function(web){
            # calculates evenness of the numbers of individual of different species in
            # a community, NOT according to formula in Mueller et al. (1999, 
            # J. Anim. Ecol), but according to the original formula in Alatalo 
            # (1981, Oikos) 
            # can be extended at some point to more indices ...
            pk <- web/sum(web)
            (Alatalo <- (1/sum(pk^2) -1) / (exp(-sum(pk * log(pk), na.rm=TRUE)) -1))
        }
        
        E <- evenness(web)
        out$"Alatalo interaction evenness" <- E
    }    
    #---------------        
    if ("Shannon diversity" %in% index) out$"Shannon diversity" <- SH    
    
}
#---------------
# Devoto's fd:
#! JFedit: changed fd to fc here (but only the function and arguments because in legacy function; this may not be the best decision)
if (any(c("fc", "functional complementarity") %in% index)){
    out$"Functional complementarity LTL" <- fc(t(web), dist=fcdist, method="average", weighted=fcweighted)		
    out$"Functional complementarity HTL" <- fc(web, dist=fcdist, method="average", weighted=fcweighted)
}

#---------------
# Bluethgen's H2'
if ("H2" %in% index){
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # from the help of is.integer!
    if (any(is.wholenumber(web)==FALSE)) H2_integer <- FALSE # turns H2_integer off if values are not integers
    H2 <- as.numeric(H2fun(web, H2_integer=H2_integer)[1]) #1.element is the standardised H2 prime
    out$"H2"= ifelse(H2<0, 0, H2)
}

#---------------------------------------------------------------------------
if (!("degree distribution" %in% index)) out <- unlist(out)
    
return(out)
    
}


#networklevel(Safariland, index="H2", legacy=F)
#networklevel(Safariland, index="ALLBUTDD", legacy=F)
#networklevel(Safariland, index=c("Alatalo interaction evenness", "H2"), legacy=F)
#networklevel(Safariland, index=c("extinction slope", "Alatalo interaction evenness", "H2"), legacy=F)
#networklevel(Safariland, index=c("H2", "degree distribution" ,"ISA", "SA", "clustering coefficient", "partner diversity"))
#networklevel(Safariland, index=c("H2", "extinction slope" ,"ISA", "SA", "cluster coefficient", "partner diversity", "diversity"))
#networklevel(Safariland, index="ALL")
#networklevel(Safariland, index=c("ISA", "SA"), level="both", weighted=FALSE)
#networklevel(Safariland,level="higher")
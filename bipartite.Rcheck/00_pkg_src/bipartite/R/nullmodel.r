nullmodel <- function(web, N=1000, method="r2d", ...){
    # function to generate null models based on a web presented
    # only a small helper function to use null models
    #
    # web   a binary or quantitative network
    # N     number of null model replicates wanted
    # method  number or name of the null model type: 1/"r2dtable", 2/swap.web, 3/vaznull, 4/shuffle.web, 5/mgen; partial match of names; methods 1 to 4 work for quantitative webs, 4 and 5 for binary.
    # ...   arguments to be passed on to null model function (see e.g. swap.web and mgen)
    #
    # shuffle.web can be used for binary and quantitative webs: for binary, it uses the function commsimulator with method "quasiswap", for quantitative, it uses the function shuffle.web.
    #
    # by Carsten F. Dormann 31-Jul-2008
    
    methods <- c("r2dtable", "swap.web", "vaznull", "shuffle.web", "mgen")
    if (is.numeric(method)){
        m <- method 
    } else {
        m <- pmatch(method, methods)
    }
    if (is.na(m)) stop("Abbreviated name does not uniquely identify method.")
    
    if (m == 1){ #r2dtable nullmodel
        if (all(web < 2))#{
            warning("This seems to be a binary web. Only methods shuffle.web and mgen should be used!\n  I proceeded nonetheless. Read the note in the help file!")
#            m <- 5
#        } else {
           rs <- rowSums(web)
           cs <- colSums(web)
           out <- r2dtable(N, r=rs, c=cs)
#	}
    }
    
    if (m == 2){# swap.web
        if (all(web < 2))#{
            warning("This seems to be a binary web. Only methods shuffle.web and mgen should be used!\n  I proceeded nonetheless. Read the note in the help file!")
#            m <- 5
#        } else { 
		out <- swap.web(N, web, ...)     
#		}
    }
    
    if (m == 3){ #vaznull
        if (all(web < 2))#{
            warning("This seems to be a binary web. Only methods shuffle.web and mgen should be used!\n  I proceeded nonetheless. Read the note in the help file!")
#            m <- 5
#        } else { 
			out <- vaznull(N, web)
#		 }
    }
    
    if (m == 4){ #shuffle.web
        if (any(web > 1)) out <- shuffle.web(web, N, ...)
        if (all(web < 2)) out <- replicate(n=N, expr=unname(commsimulator(web, method="quasiswap", ...)), simplify=FALSE) 
    }

    if (m == 5){ #mgen
       #if (any(web > 1)) warning("Discarding quantitative information! Using only the binary version of the web!")
       #binweb <- web > 0
       #pweb <- outer(rowSums(binweb)/sum(binweb), colSums(binweb)/sum(binweb), FUN="*")
       # mgen.web <- function(binweb, pweb){
          # rbinweb <- matrix(0, nrow=nrow(binweb), ncol=ncol(binweb))
          # rbinweb[sample(1:prod(dim(binweb)), size=sum(binweb), prob=pweb)] <- 1
          # out <- empty(rbinweb)
          # #if (dim(out) < dim(web)) warning("Some null models are rank-deficient! See ?mgen for details.")
          # out
       # }
       out <- replicate(n=N, unname(mgen(web, autotransform="equiprobable")), simplify=FALSE)    
    }
    
    if (!(m %in% 1:5)) stop("Please choose a valid method.")

    return(out)

} 

## example:
#nullmodel(Safariland, N=2, method=1)
#nullmodel(Safariland>0, N=2, method=4)
## analysis example:
#obs <- unlist(networklevel(Safariland, index="weighted nestedness"))
#nulls <- nullmodel(Safariland, N=1000, method=1)
#null <- unlist(sapply(nulls, networklevel, index="weighted nestedness")) #takes a while ...
#
#plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))), main="comparison of observed with null model Patefield")
#abline(v=obs, col="red", lwd=2)    
#
#praw <- sum(null>obs) / length(null)
#ifelse(praw > 0.5, 1-praw, praw)    # P-value
#
## comparison of null model 3 and 4 for binary:
#nulls3 <- nullmodel(Safariland>0, N=200, method=3)
#nulls4 <- nullmodel(Safariland>0, N=200, method=4)
#null3 <- unlist(sapply(nulls3, networklevel, index="weighted nestedness"))
#null4 <- unlist(sapply(nulls4, networklevel, index="weighted nestedness"))
#
#
#plot(density(null3), xlim=range(c(null3, null4)), lwd=2, main="comparison of null models")
#lines(density(null4), col="red", lwd=2)
#legend("topright", c("shuffle", "mgen"), col=c("black", "red"), lwd=c(2,2), bty="n", cex=1.5)
#abline(v=networklevel(Safariland>0, index="weighted nestedness"))
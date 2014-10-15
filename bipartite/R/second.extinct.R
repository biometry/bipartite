`second.extinct` <-
function(web, participant="higher", method="abun", nrep=10, details=FALSE, ext.row=NULL, ext.col=NULL){
    # returns a matrix quantifying secondary extinctions following the deletion of
    # a species from the web
    # web           an interaction web
    # participant   "high" or "low" or "both"
    # method        deletion following "random" or "abundance"
    # nrep          only for "random": number of repititions of random extinction sequence

    # if (details==TRUE & pmatch(participant, c("both", "lower", "higher"))==1)
    # {
          # warning("\nFor random extinctions of both participants extinction sequences
            # will differ in length. Simply averaging sequences can hence not be used. Thus,
            # option 'details' will be set to FALSE internally.\n")
          # details <- TRUE
    # }
    if (participant=="both" & method=="external") stop("Sorry, that won't work. When you specify the sequence, you have to choose one of the two levels. 'both' won't work.")
    if (!is.null(ext.row) & length(ext.row) != NROW(web)) stop("The length of the external row vector is different from the numbers of rows in the network!")
    if (!is.null(ext.col) & length(ext.col) != NCOL(web)) stop("The length of the external col vector is different from the numbers of cols in the network!")
    if (participant == "higher" & method=="external" & is.null(ext.col)) stop("You need to provide an external sequence of extinction for the higher trophic level!")
    if (participant == "lower" & method=="external" & is.null(ext.row)) stop("You need to provide an external sequence of extinction for the lower trophic level!")

  
    one.second.extinct <- function(web=web, participant=participant, method=method, ext.row=ext.row, ext.col=ext.col){
		    dead <- matrix(nrow=0, ncol=3)
        colnames(dead) <- c("no", "ext.lower", "ext.higher")
        m2 <- web
        i <- 1
        repeat {
            # extinct a species and count secondary extinctions:
            n <- extinction(m2, participant=participant, method=method, ext.row=ext.row, ext.col=ext.col)
        #    rowSums(n)
# test only:  (n <- extinction(m2, participant=participant, method=method, ext.row, ext.col))
            dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count=TRUE))$empty))
            if (participant == "lower" & NROW(m2) < 2) break;
            if (participant == "higher" & NCOL(m2) < 2) break;
            if (participant == "both" & min(dim(m2)) < 2) break;
            if (any(dim(n) == 1)) break;
            if (method=="external") {
            	# unelegant: subtract 1 from each value larger than the exterminated species
            	# (more elegant, but requiring more re-writing: change web according to sequence, then delete from top to bottom/left to right; this way also "degree" and "abundance" can be handled much more conveniently; TO DO!)
            	ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > ext.col[1]] - 1
            	ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > ext.row[1]] - 1
				# now get rid of the already deleted row/column:
            	ext.row <- ext.row[-1]
            	ext.col <- ext.col[-1]
            	} # remove the already removed species from the sequence
            i <- i + 1
        }
        #m2
        dead2 <- rbind(dead, c(NROW(dead)+1, NROW(m2), NCOL(m2))) # counts extinction knock-on for the last species
       	#dead2
       	# there is a completely mystifying bug somewhere, occassionally (roughly 1/50) producing a 2 for NROW(m2) but only under the settings participants="lower" and method="degree".
       	# This is a fix for these rare situations:
       	if (participant == "lower" & method== "degree"){
       		if (length(table(dead[,2])) > 1) dead2[,2] <- 1
       	}
		
      	# If I use a random sequence, sometimes it takes only 23 steps to let all pollinators go extinct, sometimes 24. So, the matrix "dead" will have different dimensions. Correct this:
        if (nrow(dead)+1 != nrow(dead2)) stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
        if (participant == "lower") supposed.length <- NROW(web) 
        if (participant == "higher") supposed.length <- NCOL(web) 
        if (participant == "both") supposed.length <- NROW(dead2)#supposed.length <- sum(dim(web)) ### was max; but obviously can be the sum of both levels

        if (NROW(dead2) != supposed.length) { # Is dead of the right length?
        	missing <- supposed.length - NROW(dead2) 
        	addit1 <- (NROW(dead2)+1):(NROW(dead2)+missing)
        	addit2n3 <- rep(0, times=missing)
        	dead2 <- rbind(dead2, as.matrix(data.frame(addit1, addit2n3, addit2n3)))
          # # ci <- which(sapply(1:3, function(x) all(dead[-nrow(dead),x]==1))==TRUE)  # a clumsy way to find out the participant subject to primary extinction: for this column, there is a 1 in every row;
          # # BETTER: 
            # ci <- pmatch(participant, c("both", "lower", "higher")) 
            # dead2 <- matrix(0, nrow=ifelse(ci==2, nrow(web), ncol(web)), ncol=3)  # new dead-matrix
            # colnames(dead2) <- colnames(dead) #use names from "dead"
            # dead2[,1] <- 1:nrow(dead2) # step numbers
            # dead2[1:nrow(dead),2:3] <- dead[1:nrow(dead),2:3] # plug matrix "dead" into "dead2", leaving the potentially longer matrix dead2 filled with zeros below nrow(dead).

        }
        return(dead2)
    }

    if (is.vector(method)) sequence = method ### In case someone provides a sequence to the method, rather than ext.row or ext.col!
    if (pmatch(method, c("abundance", "random", "degree", "external")) %in% c(1,3,4)){# i.e. if "abundance", "degree" or "external"
        out <- one.second.extinct(web=web, participant=participant, method=method, ext.row=ext.row, ext.col=ext.col)
    
		#test only:	one.second.extinct(web=web, participant=participant, method=method, ext.row=sample(9), ext.col=1:27)

    } else {
        o <- replicate(nrep, one.second.extinct(web=web, participant=participant, method=method, ext.row=ext.row, ext.col=ext.col), simplify=FALSE)
        if (details){
        	out <- o
        } else {
			lengths <- sapply(o, nrow)
            z <- o[[which.max(lengths)]]
            z[,2:3] <- 0
            for (k in 1:length(o)) { 
            	nr <- nrow(o[[k]])
            	z[1:nr, ] <- z[1:nr, ] + o[[k]]
            	rm(nr) 
            }
            out <- z/length(o)
            out[,1] <- 1:max(lengths)
        }
    }

    class(out) <- "bipartite"
    attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
    out


}
#data(Safariland)
# second.extinct(Safariland, "higher", "d")
# second.extinct(web=Safariland, participant="lower", method="external", ext.row=sample(9), ext.col=sample(27))
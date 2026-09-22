`H2fun` <-
function(web, H2_integer=TRUE){
    # function to calculate a measure of specialisation in a bipartite web: H2'
    #
    # web   a bipartite interaction web, with e.g. pollinators as columns and plants as rows
    #
    # returns the normalised H2' and its subcomponents H2max, H2min and H2uncorrected
    #
    # for web implementation see also: http://itb1.biologie.hu-berlin.de/~nils/stat/
    #
    # Example:
    # set.seed(8793)
    # web <- matrix(rpois(81, 0.5), 9)
    # H2prime(web)
    
    # code proposals by Jochen to deal with non-integer data:
    if (H2_integer & any((web %% 1) != 0)) stop("web does not contain integers! maybe you should set H2_integer to FALSE")
    

    tot <- sum(web)       #saemtliche Interaktionen im Netz
    rs <- rowSums(web)    #Interaktion der Pflanze mit saemtlichen Bestaeubern
    cs <- colSums(web)    #saemtliche Interaktionen des jeweiligen Bestaebers

   #--------------- H2 uncorrected------------
    H2uncorr = -sum(web/tot*log(web/tot), na.rm=TRUE)

   #--------------- H2 max ------------
    # Key idea here is to allocate interactions one-by-one into the places where
    # they fit best according to the non-integer optimal web. Easy! (Not.)
    exexpec <- outer(rs, cs/tot) # non-integer optimal web
    if (!H2_integer) {
    	newweb <- exexpec
    	H2_max <- -sum(newweb/tot*log(newweb/tot), na.rm=TRUE)    # first attempt of finding H2_max
    } else {
    	newweb <- floor(exexpec)  # start new web
    	# NOTE: as in all earlier versions, the first pass compares against the EMPTY web
    	# (difexp = exexpec - 0), not against newweb; kept so results stay comparable.
    	difexp <- exexpec
    	webfull <- matrix(FALSE, nrow(web), ncol(web)) # cells whose row or column is already full
    	# running marginals: only one cell changes per pass, so re-scanning the whole matrix
    	# with rowSums()/colSums() every iteration was the dominant cost here
    	rs.now <- rowSums(newweb); cs.now <- colSums(newweb); total.now <- sum(newweb)
    	first.pass <- TRUE
    	nr <- nrow(newweb)
    	while (total.now < tot) {
       	webfull[rs.now == rs, ] <- TRUE # rows/cols that have reached their marginal total
       	webfull[, cs.now == cs] <- TRUE
       	OK <- !webfull # matrix of potential cells
	       smallestpossible <- newweb==min(newweb[OK])  # find cell with lowest number of interactions (e.g. 0)
       	greatestdif <- max(difexp[smallestpossible & OK]) # find cell value with largest different between "is" and "should"
       	bestone <- which(OK & smallestpossible & difexp==greatestdif ) # find cell for all three conditions
       	if (length(bestone)>1) bestone <- sample(bestone,1) # select randomly a cell, if different are possible
       	newweb[bestone] <- newweb[bestone] + 1 # put an interaction into that cell
       	r.idx <- (bestone - 1L) %% nr + 1L; c.idx <- (bestone - 1L) %/% nr + 1L
       	rs.now[r.idx] <- rs.now[r.idx] + 1; cs.now[c.idx] <- cs.now[c.idx] + 1
       	total.now <- total.now + 1
       	if (first.pass) { difexp <- exexpec - newweb; first.pass <- FALSE
       	} else difexp[bestone] <- exexpec[bestone] - newweb[bestone]
    	}
      H2_max <- -sum(newweb/tot*log(newweb/tot), na.rm=TRUE)    # first attempt of finding H2_max
      
#    ## what may happen is that there are very few too-high-entries, but many too-low-entries!!
#    difweb <- newweb - exexpec #positive values indicate too-high-entries
#    H2_max <- -sum(newweb/tot*log(newweb/tot), na.rm=TRUE)
#
#    while (range(difweb)[1] < -.5){    # CFD, 2 Aug 2009
#        difweb <- newweb - exexpec #positive values indicate too-high-entries
#        range(difweb)
#        ## now we need to allocate entries from towards the too-low-slots:
#        towards <- which.min(difweb)
#        from <- which.max(difweb)[1]
#        newweb[towards] <- newweb[towards] + 1
#        newweb[from] <- newweb[from] - 1
#        H2_max_proposal <- -sum(newweb/tot*log(newweb/tot), na.rm=TRUE)
#        H2_max <- ifelse (H2_max_proposal > H2_max, H2_max_proposal, H2_max)
#    }
# VIOLATES MARGINAL TOTALS!! NOT USEFUL AS IT IS!!!

    
#  # Hmaxfind     IMPROVED CODE FROM JOCHEN ------------------------------------
  	if (max(exexpec)>0.3679*tot) {    # 0.3679 is the proportion yielding maximal contribution
    # warning("one cell dominates too extremely, H2max can probably not be estimated correctly")
    # further modification to match expected values even better... ; great advance, but can get caught!!  introduce random step if this happens
    	H2.of <- function(m) -sum(m/tot*log(m/tot), na.rm=TRUE)
    	
    	# This search wanders rather than converges (it changed the matrix on all 500 passes
    	# on every web tested) and the old code kept whatever the LAST pass produced. H2max is
    	# a maximum, so keep the best matrix seen and stop after `patience` fruitless passes.
    	best.web <- newweb; best.H2 <- H2.of(newweb); stall <- 0; patience <- 50
    	for (tries in 1:500) {# reduces overfits, but NOT underfits!!
      		newmx <- newweb                         # "hin- und herschieben"
      		difexp <- exexpec - newmx #newmx=newweb!
      		greatestdif <- difexp==min(difexp)
      		if (length(which(greatestdif))>1) {
          		largestvalue <- newmx==max(newmx[greatestdif])  # evtl den groessten auswaehlen
          		first <- greatestdif & largestvalue
      		} else {first <- greatestdif}   # "first" is a boolean matrix
      		newmx[first][1] <- newmx[first][1] - 1                                    # remove one interaction from one cell (with largest difference to expected values; "too large")
      		throw <- which(rowSums(first)>0)[1] ;  thcol <- which(colSums(first)>0)[1]  # find row- and column-number of removed cell (not elegant!!)
      		mr <- max(difexp[throw,])   ; mc <- max(difexp[,thcol])                         # find largest difference value in row and column to expected; "too small"
      		if (mr >= mc) {
      		  scnd <- which(difexp[throw,]==mr) [1]                         # reallocation; start in row
            newmx[throw, scnd] <- newmx[throw,scnd] + 1                      # put in cell
            thrd <- which(difexp[,scnd] == min(difexp[,scnd]))[1]
            newmx[thrd,scnd] <- newmx[thrd,scnd] - 1                      # remove interaction from "too large cell" in that column
            newmx[thrd,thcol] = newmx[thrd,thcol] + 1                    # put interaction in the cell that recovers original r/c-sums
          } else {                                                      # as above, but first reallocation in column
             scnd <- which(difexp[,thcol]==mc)[1]
             newmx[scnd,thcol] <- newmx[scnd,thcol] + 1
             thrd <- which(difexp[scnd,] == min(difexp[scnd,]))[1]
             newmx[scnd,thrd] <- newmx[scnd,thrd] - 1
             newmx[throw,thrd] <- newmx[throw,thrd] + 1
          }
      		newweb <- newmx
      		H2.now <- H2.of(newweb)
      		if (H2.now > best.H2) { 
      		  best.H2 <- H2.now
      		  best.web <- newweb
      		  stall <- 0 
      		} else {
      		  stall <- stall + 1
      		}
      		if (stall >= patience) break
      }
      newweb <- best.web
  }     # end Hmaxfind
		#  Hmax=H2(newmx)    
	} # end else-part of initial exexpec
    H2_max.improved <- -sum(newweb/tot * log(newweb/tot), na.rm=TRUE)
    H2_max <- max(H2_max, H2_max.improved) # was an ifelse(); with best-tracking above this is just the better of the two
   
    

   #--------------- H2 min ------------
    # The key idea here is that allocating the col- & row-sums into the
    # web-matrix will yield the minimum H2 for the matrix (since it is the
    # maximum of aggregation possible.
    newweb <- matrix(0,length(rs),length(cs))
    rsrest=rs; csrest=cs
    while (!isTRUE(all.equal(sum(rsrest), 0))) { ## was round(sum(rsrest), 10) != 0 (Bob O'Hara's suggestion, finally applied)
        newweb[which(rsrest==max(rsrest))[1],which(csrest==max(csrest))[1]] = min(c(max(rsrest),max(csrest)))
        rsrest=rs-rowSums(newweb)
        csrest=cs-colSums(newweb)
    }
    # Now we have our new web, with maximum aggregation of observations,
    # and hence minimum H2
    Pnew <- newweb/sum(newweb)
    H2_min <- -sum(Pnew*log(Pnew), na.rm=TRUE)

   #--------------- H2'  ------------
   if (H2uncorr < H2_min) H2_min <- H2uncorr
   if (H2_max < H2uncorr) H2_max <- H2uncorr  # ideally, this should not occur
    #ranging (between 0 and 1):
    H_2prime <- (H2_max - H2uncorr) / (H2_max - H2_min)

    #output:
    c("H2"=H_2prime, "H2min"=round(H2_min,3), "H2max"=round(H2_max,3), "H2uncorr"=round(H2uncorr,3))
    
}

#H2fun(r2dtable(1, r=rowSums(Safariland), c=colSums(Safariland))[[1]])


# This is still a bit unsatisfactory: First, we build a web which contains too many too-high-values, then we correct it. Ideally, the check for BOTH, too-high and too-low values would be done in the building of the web. That would certainly save computational time!!
# Currently, the build-up is relatively complicated and slow because marginal totals need to be observed.


# New test cases (bug in H2fun, H2maxfind, up to version 2.13):
# myweb <- matrix(c(0,1,5,2), nrow=2, byrow=T) # this and similar webs don't work with H2fun!! can lead to NaN or "wrong" 0s
# H2fun(matrix(c(0,1,5,2), nrow=2, byrow=T)) # here, the maximum is not found (only the same as the min is found) so H2=NaN
# the web that it should find is
# bestweb <- matrix(c(1,0,4,3),nrow=2,byrow=T)
# H2fun(bestweb)  # yep, H2'=0, H2max is 0.974, but just because the function uses the observed web for H2max if that is better then the one found

# conclusion: yes, there is a long-standing bug in H2fun! the "improvement" introduced unknown time after bipartite0.3 actually makes things worse at least sometimes (in above example web, H2max-web is found in basic H2maxfind, and then in the second step lost again!!)

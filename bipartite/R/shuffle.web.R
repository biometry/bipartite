`shuffle.web` <-
function(web, N, legacy=TRUE){
  # shuffles entries of an interaction web, thereby maintaining connectence, but
  # changing marginal sums
  # to maintain dimensionality, interactions are first allocated to the diagonal
  # (thereby having all species in the web), then to the rest
  web <- as.matrix(web)
  web <- empty(web)
  
  if (legacy == FALSE){
  	if (dim(web)[1] > dim(web)[2]){ #test if the matrix has more rows than columns   (i.e. is 'long')
      long=TRUE  #flag for a 'long' matrix
      web=t(web)  #make it 'wide'
    } else {
    	long=FALSE #flag for a 'wide' matrix
    }  
  }

  shuffle <- function(web, legacy){
      dimdiff <- dim(web)[1]-dim(web)[2]
      # more columns than rows, i.e. dimdiff is negative
      if ((length(diag(web))+dimdiff) > sum(as.vector(web)>0))  stop("Too few entries in the web: less interactions than length of web diagonal.")
      out <- web
      out[,] <- 0
      shuf <- sample(as.vector(web)) #shuffle web entries
      nozero.index <- which(shuf!=0) #pick non-zeros
      diag(out) <- shuf[nozero.index[1:length(diag(out))]] #allocate to diagonal of new web
    
      if (legacy){
# find remaining cols/rows beyond the diagonal:
      if (dimdiff < 0) {colindex <- ((length(diag(out))+1) : dim(web)[2])} else {colindex <- 1:dim(web)[2]}
      if (dimdiff > 0) {rowindex <- ((length(diag(out))+1) : dim(web)[1])} else {rowindex <- 1:dim(web)[1]}
    
      if (dimdiff < 0) rowposition <- sample(rowindex, length(colindex), replace=TRUE)
      if (dimdiff > 0) colposition <- sample(colindex, length(rowindex), replace=TRUE)
    
      # slot in non-zeros:
      if (dimdiff < 0) {
          for (i in 1:abs(dimdiff)){
            out[rowposition[i], colindex[i]] <- shuf[nozero.index[(length(diag(web))+i)]]
          }
      }
      if (dimdiff > 0){
          for (j in 1:abs(dimdiff)){
           out[rowindex[j], colposition[j]] <- shuf[nozero.index[(length(diag(web))+j)]]
          }
      }

	  # number of interactions already allocated:
	  gone <- sum(out>0)	
    
      # allocate remaining entries:
	  remains <- shuf[nozero.index[-c(1:gone)]]
#      remains <- shuf[!(shuf %in% out)]  
      option <- which(out == 0, arr.ind=TRUE)
      out[option[sample(dim(option)[1], length(remains)),]] <- remains
      colnames(out) <- rownames(out) <- NULL
      out
	} else { # this is the second part of the "legacy"-condition:
	# Here starts the code by Paul Rabie, 8 April 2011, which is faster than the old.
    	xdiag.index=matrix(ncol=2,nrow=max(dim(web))-min(dim(web)))
       xdiag.index[,1]=sample(c(1:(min(dim(web)))),size=max(dim(web))-min(dim(web)),replace=TRUE)
       xdiag.index[,2]=(min(dim(web))+1):max(dim(web)) #xdiag is a matrix containing one cell per 'extra' columns (extra, relative to the diagonal)  This removes the need for a loop
       out[xdiag.index]=shuf[nozero.index[(length(diag(web)) + 1):(length(diag(web))+nrow(xdiag.index))]]
       gone <- sum(out > 0)
       remains <- shuf[nozero.index[-c(1:gone)]]
       option <- which(out == 0, arr.ind = TRUE)
       out[option[sample(dim(option)[1], length(remains)), ]] <- remains
       colnames(out) <- rownames(out) <- NULL
       if(long){
      	   out=t(out) #if the input matrix was long, transform the output matix  so it is also long
       }  
       out    
	} # end of legacy
  } # end of shuffle function
  replicate(N, shuffle(web, legacy=legacy), simplify=FALSE)
}
# data(Safariland)
# system.time(shuffle.web(Safariland, 100))
# system.time(shuffle.web(Safariland, 100, legacy=FALSE))

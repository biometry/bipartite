wine <- function(web, nreps=1){
       
      # function to compute WIN
      wincomp <- function(M, Nr, Nc,Nl, xi, yj){  
                M.packed <- M[order(rowSums(M > 0)), order(colSums(M > 0)), drop = FALSE]
                Pc <-   M.packed / matrix(rep(rowSums(M.packed),Nc),Nr,Nc) 
                Pr <-    M.packed / matrix(rep(colSums(M.packed),Nr),Nr,Nc, byrow=TRUE) 
                dij.w <- Pr * matrix(rep(xi, Nc), Nr, Nc) + Pc * matrix(rep(yj, Nr), Nr, Nc, byrow=TRUE) 
                win <- sum(dij.w)/Nl
                return(list(win=win, dij.w=dij.w))
      }


      # function to randomize equipropably constrained to left al margins > 0
      shuffle.equiprob <- function(M){
               
                randiag <- function(vdims){
                    mats <- diag(vdims)
                    posdiag <- 1:vdims +((1:vdims -1)*vdims) 
                    desp <- matrix(rep(1:vdims,vdims),vdims,vdims, byrow=TRUE)-(1:vdims)
                    fdesp <- sample(1:vdims)
                    move <- desp[cbind(1:vdims,fdesp)]
                    moved <- posdiag + move
                    mdesp <- matrix(0, vdims,vdims)
                    mdesp[moved] <- 1
                    return(mdesp)
                }
                M<- as.matrix(M)
                values <- M[which(M>0)] 
                lvalues <- length(values) 
                if(identical(dim(M)[1],dim(M)[2])){
                    vdims <-dim(M)[1]
                    vdimb <-dim(M)[1]
                  
                }                    
                    
                if(!identical(dim(M)[1],dim(M)[2])){
                    dims <- which(dim(M)==min(dim(M))) 
                    vdims <- dim(M)[dims] 
                    dimb <- which(dim(M)==max(dim(M))) 
                    vdimb <- dim(M)[dimb]
                }                    

                MR <- matrix(0, vdims, vdimb)  
                lMR <- vdims * vdimb  
 
                sample1 <- sample(vdimb, vdims)
                diag1 <- randiag(vdims)
                MR[,sample1] <- diag1
                sample2<- (1:vdimb)[-sample1]
                pos <- sample(vdims, length(sample2), replace=TRUE)
                MR[cbind(pos,sample2)] <-2
                MRoccupied <- which(MR>0)
                vleft <- lvalues- vdimb
                if(vleft>0) MRoccupied <- c(MRoccupied, sample((1:lMR)[-MRoccupied], vleft))
 
                MR[MRoccupied] <- sample(values)

                if(dim(MR)[1] != dim(M)[1]) MR <- t(MR)

                return(MR)
      }


      # function to compute maximum nestedness matrix
       mpack<- function(M){
                dim1 <- dim(M)[1]
                dim2 <- dim(M)[2]
                M_pack <- matrix(0,dim1,dim2)
                M_pos <- M[M>0]
                i_val <- length(M_pos)
                val<-sort(M_pos)
                for(i in 1:dim1){
                        M_pack[dim1-i+1,dim2] <- val[i_val]

                    i_val <- i_val-1
                }
                ultif <- dim1
                ultic <- dim2-1
                while(i_val>0){
                      for(i in 1:min(ultic,i_val)){
                            M_pack[ultif,ultic-i+1] <- val[i_val]
                          i_val <- i_val-1
                      }
                      ultif <- ultif-1
                      if(i_val>0){
                            for(i in 1:min(ultif,i_val)){
                                     M_pack[ultif-i+1,ultic] <- val[i_val]
                                    i_val <- i_val-1
                            }
                            ultic <- ultic-1     
                      }
                }
                return(M_pack)
      }



      # Compute WIN et al 
      M <- web
          c0 <- colSums(M)==0
          r0 <- rowSums(M)==0
          if(sum(c0,r0) > 0) {
              M<- as.matrix(M[!r0,!c0])
              warning("one or more rows and/or columns are completely made out of 0s and have been removed")
          }   
      
      Nr  <- dim(M)[1]
      Nc <- dim(M)[2]
      Nl <- length(which(M!=0))
      i <- 1:Nr
      j <- 1:Nc
      xi <- (i-1)/Nr + 1/(2*Nr)
      yj <- (j-1)/Nc + 1/(2*Nc)
       
      dw <- wincomp (M=M, Nr=Nr, Nc=Nc, Nl=Nl, xi=xi, yj=yj)
      dij.w <- dw$dij.w 
      dw <- dw$win 
      dmax <- wincomp (M=mpack(M), Nr=Nr, Nc=Nc, Nl=Nl, xi=xi, yj=yj)
      dij.max <- dmax$dij.w
      dmax <- dmax$win
      d.rnd <-NULL
      for (i in 1: nreps){
            MR <- shuffle.equiprob(M) 
            d.rnd <- c(d.rnd, wincomp (M=MR, Nr=Nr, Nc=Nc, Nl=Nl, xi=xi, yj=yj)$win)
      }
         
      z_score <- (dw-mean(d.rnd))/ sd(d.rnd)
      p_value <- 1 - pnorm(dw,mean(d.rnd, na.rm=TRUE), sd(d.rnd,na.rm=TRUE))
      eta.w <- (dw- mean(d.rnd)) / (dmax- mean(d.rnd))

    result <- list(win=dw, wine =eta.w, zscore=z_score, pvalue=p_value,
                    dmax=dmax, drnd=mean(d.rnd), dij.w=dij.w,
                    dij.max=dij.max)
    class(result) <- c("wine", class(result))
    return(result)
}

#-------------------------------------------------------------------------------
plot.wine<- function(x,...){
          w <- t(x$dij.w)
          dim1 <- dim(w)[1]
          dim2 <- dim(w)[2]
          image.plot(x=1:dim1, y=1:dim2, z=w[,dim2:1], 
            main="Weighted distances matrix",xlab="Column",ylab="Row",yaxt="n",xaxt="n",...)
            ticks_x<-1:dim1  
            ticks_y<-1:dim2
            abline(h=0)
            axis(1, at= ticks_x, labels=ticks_x) 
            abline(v=0)
            axis(2, at= ticks_y, labels=rev(ticks_y), las=2) 
}

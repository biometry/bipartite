plotmatrix <- function(x, background_color="white", base_color=NULL, between_color="black", border_color="black", modules_colors=NULL, within_color = "black", border = FALSE, row_partitions=NULL, col_partitions=NULL, binary=TRUE, plot_labels=FALSE, xlab=NA, ylab=NA, offset = 0.4, ...){
  
  if (is.list(x)){ # cfd
    # for the probably many cases where users provide the sortweb result, rather than only its matrix
    row_partitions <- x$row_partitions
    col_partitions <- x$col_partitions
    mat <- x$matrix 
  } else {
    mat <- x
  }
  mat <- as.matrix(mat)
  
  if (is.null(base_color)){
    base_color=background_color
  }
  
  if (!is.null(row_partitions)|!is.null(col_partitions)){
    if(!identical(sort(unique(row_partitions)),sort(unique(col_partitions)))){
      stop("different numbers of col and row partitions")
    }
    if(length(row_partitions)!=nrow(mat)|length(col_partitions)!=ncol(mat)){
      stop ("partitions with inapropriate length")
    }
  }

  
  if (!is.null(row_partitions)){
    index=unique(row_partitions)
    row_partitions=match(row_partitions, index)
    col_partitions=match(col_partitions, index)
  }
  
  if (!is.null(row_partitions)){
    row_partitions=as.numeric(row_partitions)
    col_partitions=as.numeric(col_partitions)
  }
  
  ## Checking inputs
  if (border==TRUE & is.null(row_partitions)){
    border=FALSE
    warning("Borders cannot be plotted if partitions are not defined.")
  }
  if (is.null(row_partitions)){modules_colors=NULL}
  if (!is.null(modules_colors) & length(unique(row_partitions)) != length(modules_colors)){
    stop("wrong number of modules colors")
  }
  
  ####
  #Defining colors
  if(!is.null(modules_colors)){within_color=NULL}
  colors=c(between_color, within_color, modules_colors)
  #mat2 = interaction matrices in which numbers defines the color of the rectangle
  mat2=mat
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      if(!is.null(row_partitions)){
        if(is.null(modules_colors)){
          if(mat[i,j]==0){mat2[i,j]=0} else{
            mat2[i,j]=ifelse(row_partitions[i]==col_partitions[j],2,1)
          }}
        if(!is.null(modules_colors)){
          if(mat[i,j]==0){mat2[i,j]=0} else{
            mat2[i,j]=ifelse(row_partitions[i]!=col_partitions[j],1,row_partitions[i]+1)}
        }
      }
    }
  }
  image(y=1:nrow(mat2), x=1:ncol(mat2), z=0*t(mat2), col = background_color, xaxt = "n", yaxt = "n", xlab=xlab, ylab=ylab)
  if (plot_labels){axis(1,at = 1:ncol(mat2), labels = colnames(mat), ...)}
  if (plot_labels){axis(2,at = 1:nrow(mat2), labels = rownames(mat)[length(rownames(mat)):1], ...)}
  box()
  
  TRmat2 <- t(mat2[ nrow(mat2):1, ] )
  
  if (!binary){
    TRmat <- t(mat[ nrow(mat):1, ] )
    TRmat=TRmat/max(TRmat)
  }
  
  offset=offset
  
  #plotting filled positions as rectangles
  for (a in 1:nrow(TRmat2)) {
    for (b in 1:ncol(TRmat2)) { 
      if (TRmat2[a,b] != 0){
        if (is.null(row_partitions)){
          if(binary){rect(a-offset, b-offset, a+offset, b+offset, col=within_color, border='NA')}
          if(!binary){
            MAX=max(mat)
            funcol=colorRamp(c(base_color,within_color))
            rect_color=funcol(TRmat[a,b])
            rect(a-offset,b-offset,a+offset,b+offset,col=rgb(rect_color[1],rect_color[2],rect_color[3],255,maxColorValue = 255),border='NA')}
          
        } else { 
          if (binary){rect(a-offset,b-offset,a+offset,b+offset,
                             col=colors[TRmat2[a,b]],border='NA')}
          if (!binary){
            MAX=max(mat)
            funcol=colorRamp(c(base_color,colors[TRmat2[a,b]]))
            rect_color=funcol(TRmat[a,b])
            rect(a-offset, b-offset, a+offset, b+offset, col=rgb(rect_color[1], rect_color[2], rect_color[3], 255, maxColorValue = 255), border='NA')
          }
        }
      }
    }
  }
  # plotting borders
  if (border) {
    c=0
    for (i in 1:length(mat)){
      if (i%%nrow(mat)==1){c=c+1}
      r=i%%nrow(mat)
      if (r==0){r=nrow(mat)}
      mod=row_partitions[r]==col_partitions[c]
      if (r==1){
        modtop=FALSE
      }else{
        modtop=row_partitions[r-1]==col_partitions[c]
      }
      if (c==1){
        modleft=FALSE
      }else{
        modleft=row_partitions[r]==col_partitions[c-1]
      }
      
      if(mod & !modtop & !modleft){
        xleft = c-0.5
        ytop = nrow(mat) - r + 1.5
        endmod = FALSE
        j=c
        while(!endmod){
          j=j+1
          if ((j-1)==ncol(mat)){
            endmod=T
          } else {
            endmod=col_partitions[j]!=col_partitions[j-1]
          }
        }
        xright = j-0.5
        endmod = FALSE
        j = r
        while(!endmod){
          j=j+1
          if((j-1)==nrow(mat)){
            endmod=TRUE
          } else {
            endmod=row_partitions[j]!=row_partitions[j-1]
          }
        }
        ybottom=nrow(mat)-j+1.5
        rect(xleft = xleft,xright = xright,ybottom = ybottom,ytop = ytop,col="NA",border=border_color)
      }
    }
    
  }
}
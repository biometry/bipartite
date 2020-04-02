sortmatrix <- function(matrix, topology="compound", sort_by="degrees", row_partitions=NULL, col_partitions=NULL, mod_similarity=FALSE){
  
  matrix <- as.matrix(matrix)
  
  if(!is.null(row_partitions)){
    index <- unique(row_partitions)
    row_partitions2 = match(row_partitions,index)
    col_partitions2 = match(col_partitions,index)
  }
  
  ## Checking inputs
  if(topology!="nested"){
    if(!is.null(row_partitions2)|!is.null(col_partitions2)){
      if(!identical(sort(unique(row_partitions2)),sort(unique(col_partitions2)))){
        stop("different numbers of col and row partitions")
      }
      if(length(row_partitions2)!=nrow(matrix)|length(col_partitions2)!=ncol(matrix)){
        stop ("partitions with inapropriate length")
      }
    }
  }
  
  if(topology!="nested" & topology!="compound" & topology!="modular"){
    stop("topology must be: nested, modular or compound")
  }
  if(sort_by!="degrees" & sort_by!="weights"){
    stop("sort_by must be: degrees or weights")
  }
  ####
  
  if (sort_by=="degrees"){
    MATRIX=1*(matrix!=0)
  }
  if (sort_by=="weights"){
    MATRIX=matrix
  }
  
  if (topology=="nested"){
    order_row=order(rowSums(MATRIX),decreasing = TRUE)
    order_col=order(colSums(MATRIX),decreasing = TRUE)
  }
  if (topology=="modular"){
    if (mod_similarity==FALSE){
      order_row=order(row_partitions2)
      order_col=order(col_partitions2)
    }
    if (mod_similarity){
      x = matrix(NA,length(unique(row_partitions2)),length(unique(row_partitions2)))
      for (rr in 1:nrow(x)){
        for (cc in 1:rr){
          x[rr,cc] = mean(c(MATRIX[row_partitions2==rr,col_partitions2==cc],MATRIX[row_partitions2==cc,col_partitions2==rr]))
        }
      }
      x1=(max(x,na.rm = TRUE) - as.dist(x))/max(x, na.rm=TRUE)
      y=hclust(x1)
      part_order=y$order
      row_partitions3=match(row_partitions2,part_order)
      col_partitions3=match(col_partitions2,part_order)
      order_row=order(row_partitions3)
      order_col=order(col_partitions3)
    }
  }
  if (topology=="compound"){
    if (!mod_similarity){
      part_row=row_partitions2[order(rowSums(MATRIX),decreasing=TRUE)]
      order_row=order(rowSums(MATRIX),decreasing=TRUE)
      order_row=order_row[order(part_row)]
      part_col=col_partitions2[order(colSums(MATRIX), decreasing=TRUE)]
      order_col=order(colSums(MATRIX), decreasing=TRUE)
      order_col=order_col[order(part_col)]
    }
    if (mod_similarity){
      x = matrix(NA, length(unique(row_partitions2)), length(unique(row_partitions2)))
      for (rr in 1:nrow(x)){
        for (cc in 1:rr){
          x[rr,cc] = mean(c(MATRIX[row_partitions2==rr, col_partitions2==cc], MATRIX[row_partitions2==cc, col_partitions2==rr]))
        }
      }
      x1 = (max(x, na.rm=TRUE)-as.dist(x))/max(x, na.rm=TRUE)
      y=hclust(x1)
      part_order=y$order
      row_partitions3=match(row_partitions2,part_order)
      col_partitions3=match(col_partitions2,part_order)
      
      part_row=row_partitions3[order(rowSums(MATRIX),decreasing = TRUE)]
      order_row=order(rowSums(MATRIX),decreasing = TRUE)
      order_row=order_row[order(part_row)]
      part_col=col_partitions3[order(colSums(MATRIX),decreasing = TRUE)]
      order_col=order(colSums(MATRIX),decreasing = TRUE)
      order_col=order_col[order(part_col)]
    }
  }
  finalmatrix=matrix[order_row,order_col]
  if(!is.null(rownames(matrix))){
    row_names= rownames(matrix)[order_row]
    rownames(finalmatrix)=row_names
  }
  if(!is.null(colnames(matrix))){
    col_names= colnames(matrix)[order_col]
    colnames(finalmatrix)=col_names
  }
  result=list(matrix=finalmatrix,
              row_partitions=row_partitions[order_row],
              col_partitions=col_partitions[order_col],
              order_row=order_row,
              order_col=order_col)
  return(result)
}

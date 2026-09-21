nest.smdm <- function(x, constraints=NULL, weighted=FALSE, decreasing="fill", sort=TRUE){
  ### Checking inputs ####
  if (!is.null(constraints)&length(unique(constraints))==1){
    warning("Only one module. Nestedness calculated only for the entire matrix")
    constraints = NULL
  }
  if(is.element(NA, constraints) | is.element(NaN, constraints)){
    warning("NA or NaN in constraints. Nestedness calculated only for the entire matrix")
    constraints = NULL
  }
  if (!is.null(constraints)&length(constraints)!=nrow(x)+ncol(x)){
    stop("constraints vector is not of the same length that network vertices")
  }
  if (weighted == FALSE & any(x != 0 & x != 1)){
    x[x>0] = 1
    warning ("binary metric applied")
  }
  if (decreasing!="fill" & decreasing!="abund"){
    stop("decreasing should be fill or abund")
  }
  if (!is.null(constraints)){constraints = as.character(constraints)}
  
  # check and potentially create dimnames:
  if (is.null(rownames(x))){
    xrnames <- paste("R", 1:nrow(x), "")
    rownames(x) <- xrnames
  }
  if (is.null(colnames(x))){
    xcnames <- paste("C", 1:ncol(x), "")
    colnames(x) <- xcnames
  }
  

  ### Vectorised pair counters #################################################
  # The six scalar loops these replace were of the form
  #     for (j in 2:n) for (i in 1:(j-1)) { S <- 0; for (k in 1:m) if (<cond>) S <- S+1 }
  # i.e. O(n^2 * m) interpreted operations. The same counts come out of one matrix
  # product (binary case) or one pass over columns (weighted case).

  # binary: number of columns where BOTH rows carry a 1
  cooccurrence.count <- function(m) tcrossprod(m == 1) * 1

  # weighted: number of columns where row j is non-zero AND strictly smaller than row i
  dominance.count <- function(m){
    S <- matrix(0, nrow(m), nrow(m))
    for (k in seq_len(ncol(m))){
      v <- m[, k]
      S <- S + outer(v, v, function(aj, ai) (aj != 0) & (aj < ai))
    }
    S
  }

  # fill the lower triangle: 0 where cmp[j] >= cmp[i], else S[j,i] * 100 / denom[j]
  fill.N <- function(S, denom, cmp){
    n <- length(denom)
    dj <- matrix(denom, n, n, byrow=FALSE)
    cj <- matrix(cmp,   n, n, byrow=FALSE)
    ci <- matrix(cmp,   n, n, byrow=TRUE)
    vals <- ifelse(cj >= ci, 0, S * 100 / dj)
    N <- matrix(NA_real_, n, n)
    lo <- lower.tri(N)
    N[lo] <- vals[lo]
    N
  }

  ### Unweighted NODF Function ####
  unweightednodf = function (x, constraints){
    # Sorting matrix order by row and collumn sums
    if (sort==TRUE){tab0=x[sort(rowSums(x), index=TRUE, decreasing=TRUE)$ix,
                        sort(colSums(x), index=TRUE, decreasing=TRUE)$ix]}
    else {tab0=x}
    
    # N for rows
    MTrow = rowSums(tab0)
    Nrow = matrix(rep(NA, times=nrow(tab0)^2), nrow(tab0), nrow(tab0))
    dimnames(Nrow)=list(rownames(tab0), rownames(tab0))
    
    Nrow[] <- fill.N(cooccurrence.count(tab0),  denom=MTrow, cmp=MTrow)
    Nrow = Nrow[rownames(x), rownames(x)]
    
    # NODF for rows
    NODFrow = mean(Nrow, na.rm = TRUE)
    
    # N for collumns
    
    MTcol = colSums(tab0)
    Ncol = matrix(rep(NA, times=ncol(tab0)^2), ncol(tab0), ncol(tab0))
    dimnames(Ncol) = list(colnames(tab0), colnames(tab0))
    
    Ncol[] <- fill.N(cooccurrence.count(t(tab0)), denom=MTcol, cmp=MTcol)
    Ncol=Ncol[colnames(x),colnames(x)]
    
    # NODF for rows
    NODFcol= mean(Ncol,na.rm = TRUE)
    
    # NODF for the entire matrix
    NODFmatrix= mean(c(Ncol,Nrow),na.rm=TRUE)
    
    #### NODF SM/DM ###
    if (!is.null(constraints)){
      # constraints for rows
      
      rowcons=cbind (rownames(x),constraints[1:nrow(x)])
      tabrcons=table(rowcons[,1],rowcons[,2])
      distrcons= dist(tabrcons,method = "binary")
      distrcons= as.matrix (distrcons)
      distrcons=distrcons[rownames(x),rownames(x)]
      rm(rowcons,tabrcons)
      
      # NODF SM/DM for rows
      SM_Nrow=0
      SM_nrow=0
      DM_Nrow=0
      DM_nrow=0
      for (i in 1:nrow(x)){
        for (j in 1:nrow(x)){
          if (!is.na(Nrow[i,j])){
            if(distrcons[i,j]==0){
              SM_Nrow=SM_Nrow+Nrow[i,j]
              SM_nrow=SM_nrow+1
            }
            else{
              DM_Nrow=DM_Nrow+Nrow[i,j]
              DM_nrow=DM_nrow+1
            }
          }
        }
      }
      NODF_SM_row= SM_Nrow/SM_nrow
      NODF_DM_row= DM_Nrow/DM_nrow
      
      # constraints for columns
      
      colcons=cbind (colnames(x),constraints[(nrow(x)+1):length(constraints)])
      tabccons=table(colcons[,1],colcons[,2])
      distccons= dist(tabccons,method = "binary")
      distccons= as.matrix (distccons)
      distccons=distccons[colnames(x),colnames(x)]
      rm(colcons,tabccons)
      
      # NODF SM/DM for collumns
      SM_Ncol=0
      SM_ncol=0
      DM_Ncol=0
      DM_ncol=0
      for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
          if (!is.na(Ncol[i,j])){
            if(distccons[i,j]==0){
              SM_Ncol=SM_Ncol+Ncol[i,j]
              SM_ncol=SM_ncol+1
            }
            else{
              DM_Ncol=DM_Ncol+Ncol[i,j]
              DM_ncol=DM_ncol+1
            }
          }
        }
      }
      NODF_SM_col= SM_Ncol/SM_ncol
      NODF_DM_col= DM_Ncol/DM_ncol
      
      # NODF SM/DM for matrix
      
      NODF_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
      NODF_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
      # return
      return(list(NODFrow=NODFrow,NODFcol=NODFcol, NODFmatrix=NODFmatrix,
                  NODF_SM_row= NODF_SM_row, NODF_DM_row=NODF_DM_row, 
                  NODF_SM_col= NODF_SM_col, NODF_DM_col=NODF_DM_col,
                  NODF_SM_matrix= NODF_SM_matrix, NODF_DM_matrix=NODF_DM_matrix))
      
    }
    else {
      return(list(NODFrow=NODFrow,NODFcol=NODFcol, NODFmatrix=NODFmatrix))}
  }
  
  ### Weighted NODF function ####
  weightednodf=function (x, constraints){
    # Sorting matrix order by row and collumn sums
    if(sort==TRUE){tab0=x[sort(rowSums(x!=0), index=TRUE, decreasing=TRUE)$ix,
                       sort(colSums(x!=0), index=TRUE, decreasing=TRUE)$ix]}
    else{tab0=x}
    
    # N for rows
    MTrow= rowSums(tab0)
    Frow= rowSums(tab0!=0)
    Nrow= matrix(rep(NA,times=nrow(tab0)^2),nrow(tab0),nrow(tab0))
    dimnames(Nrow)=list(rownames(tab0),rownames(tab0))
    
    Nrow[] <- fill.N(dominance.count(tab0),    denom=Frow, cmp=Frow)
    Nrow=Nrow[rownames(x), rownames(x)]
    
    # WNODF for rows
    NODFrow= mean(Nrow,na.rm = TRUE)
    
    # N for collumns
    
    MTcol= colSums(tab0)
    Fcol= colSums(tab0!=0)
    Ncol= matrix(rep(NA,times=ncol(tab0)^2),ncol(tab0),ncol(tab0))
    dimnames(Ncol)=list(colnames(tab0),colnames(tab0))
    
    Ncol[] <- fill.N(dominance.count(t(tab0)), denom=Fcol, cmp=Fcol)
    Ncol=Ncol[colnames(x),colnames(x)]
    
    # WNODF for rows
    NODFcol= mean(Ncol,na.rm = TRUE)
    
    # WNODF for the entire matrix
    NODFmatrix= mean(c(Ncol,Nrow),na.rm=TRUE)
    
    #### WNODF SM/DM ###
    if (!is.null(constraints)){
      # constraints for rows
      
      rowcons=cbind (rownames(x),constraints[1:nrow(x)])
      tabrcons=table(rowcons[,1],rowcons[,2])
      distrcons= dist(tabrcons,method = "binary")
      distrcons= as.matrix (distrcons)
      distrcons=distrcons[rownames(x),rownames(x)]
      rm(rowcons,tabrcons)
      
      # WNODF SM/DM for rows
      SM_Nrow=0
      SM_nrow=0
      DM_Nrow=0
      DM_nrow=0
      for (i in 1:nrow(x)){
        for (j in 1:nrow(x)){
          if (!is.na(Nrow[i,j])){
            if(distrcons[i,j]==0){
              SM_Nrow=SM_Nrow+Nrow[i,j]
              SM_nrow=SM_nrow+1
            }
            else{
              DM_Nrow=DM_Nrow+Nrow[i,j]
              DM_nrow=DM_nrow+1
            }
          }
        }
      }
      NODF_SM_row= SM_Nrow/SM_nrow
      NODF_DM_row= DM_Nrow/DM_nrow
      
      # constraints for collumns
      
      colcons=cbind (colnames(x),constraints[(nrow(x)+1):length(constraints)])
      tabccons=table(colcons[,1],colcons[,2])
      distccons= dist(tabccons,method = "binary")
      distccons= as.matrix (distccons)
      distccons=distccons[colnames(x),colnames(x)]
      rm(colcons,tabccons)
      
      # WNODF SM/DM for collumns
      SM_Ncol=0
      SM_ncol=0
      DM_Ncol=0
      DM_ncol=0
      for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
          if (!is.na(Ncol[i,j])){
            if(distccons[i,j]==0){
              SM_Ncol=SM_Ncol+Ncol[i,j]
              SM_ncol=SM_ncol+1
            }
            else{
              DM_Ncol=DM_Ncol+Ncol[i,j]
              DM_ncol=DM_ncol+1
            }
          }
        }
      }
      NODF_SM_col= SM_Ncol/SM_ncol
      NODF_DM_col= DM_Ncol/DM_ncol
      
      # WNODF SM/DM for matrix
      NODF_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
      NODF_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
      # return
      return(list(WNODFrow=NODFrow,WNODFcol=NODFcol, WNODFmatrix=NODFmatrix,WNODF_SM_row= NODF_SM_row, WNODF_DM_row=NODF_DM_row,WNODF_SM_col= NODF_SM_col, WNODF_DM_col=NODF_DM_col,WNODF_SM_matrix= NODF_SM_matrix, WNODF_DM_matrix=NODF_DM_matrix))
      
    }
    else {
      return(list(WNODFrow=NODFrow,WNODFcol=NODFcol, WNODFmatrix=NODFmatrix))}
  }
  
  ### Weighted NODA funcion ####
  weightednoda=function (x,constraints){
    # Sorting matrix order by row and collumn sums
    if (sort == TRUE){
      tab0=x[sort(rowSums(x), index=TRUE, decreasing=TRUE)$ix,
                       sort(colSums(x), index=TRUE, decreasing=TRUE)$ix]
    } else { tab0 <- x }
    
    # N for rows
    MTrow= rowSums(tab0)
    Frow= rowSums(tab0!=0)
    Nrow= matrix(rep(NA,times=nrow(tab0)^2),nrow(tab0),nrow(tab0))
    dimnames(Nrow)=list(rownames(tab0),rownames(tab0))
    
    Nrow[] <- fill.N(dominance.count(tab0),    denom=Frow, cmp=MTrow)
    Nrow=Nrow[rownames(x), rownames(x)]
    
    # WNODA for rows
    NODArow= mean(Nrow,na.rm = TRUE)
    
    # N for collumns
    
    MTcol= colSums(tab0)
    Fcol= colSums(tab0!=0)
    Ncol= matrix(rep(NA,times=ncol(tab0)^2),ncol(tab0),ncol(tab0))
    dimnames(Ncol)=list(colnames(tab0),colnames(tab0))
    
    Ncol[] <- fill.N(dominance.count(t(tab0)), denom=Fcol, cmp=MTcol)
    Ncol=Ncol[colnames(x),colnames(x)]
    
    # NODA for rows
    NODAcol= mean(Ncol,na.rm = TRUE)
    
    # NODA for the entire matrix
    NODAmatrix= mean(c(Ncol,Nrow),na.rm=TRUE)
    
    #### WNODA SM/DM ###
    if (!is.null(constraints)){
      
      # constraints for rows
      rowcons=cbind (rownames(x),constraints[1:nrow(x)])
      tabrcons=table(rowcons[,1],rowcons[,2])
      distrcons= dist(tabrcons,method = "binary")
      distrcons= as.matrix (distrcons)
      distrcons=distrcons[rownames(x),rownames(x)]
      rm(rowcons,tabrcons)
      
      # WNODA SM/DM for rows
      SM_Nrow=0
      SM_nrow=0
      DM_Nrow=0
      DM_nrow=0
      for (i in 1:nrow(x)){
        for (j in 1:nrow(x)){
          if (!is.na(Nrow[i,j])){
            if(distrcons[i,j]==0){
              SM_Nrow=SM_Nrow+Nrow[i,j]
              SM_nrow=SM_nrow+1
            }
            else{
              DM_Nrow=DM_Nrow+Nrow[i,j]
              DM_nrow=DM_nrow+1
            }
          }
        }
      }
      NODA_SM_row= SM_Nrow/SM_nrow
      NODA_DM_row= DM_Nrow/DM_nrow
      
      # constraints for collumns
      
      colcons=cbind (colnames(x),constraints[(nrow(x)+1):length(constraints)])
      tabccons=table(colcons[,1],colcons[,2])
      distccons= dist(tabccons,method = "binary")
      distccons= as.matrix (distccons)
      distccons=distccons[colnames(x),colnames(x)]
      rm(colcons,tabccons)
      
      # WNODA SM/DM for collumns
      SM_Ncol=0
      SM_ncol=0
      DM_Ncol=0
      DM_ncol=0
      for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
          if (!is.na(Ncol[i,j])){
            if(distccons[i,j]==0){
              SM_Ncol=SM_Ncol+Ncol[i,j]
              SM_ncol=SM_ncol+1
            }
            else{
              DM_Ncol=DM_Ncol+Ncol[i,j]
              DM_ncol=DM_ncol+1
            }
          }
        }
      }
      NODA_SM_col= SM_Ncol/SM_ncol
      NODA_DM_col= DM_Ncol/DM_ncol
      
      # WNODA SM/DM for matrix
      
      NODA_SM_matrix= (SM_Nrow+SM_Ncol)/(SM_nrow+SM_ncol)
      NODA_DM_matrix= (DM_Nrow+DM_Ncol)/(DM_nrow+DM_ncol)
      # return
      return(list(WNODArow=NODArow,WNODAcol=NODAcol, WNODAmatrix=NODAmatrix,
                  WNODA_SM_row= NODA_SM_row, WNODA_DM_row=NODA_DM_row, 
                  WNODA_SM_col= NODA_SM_col, WNODA_DM_col=NODA_DM_col,
                  WNODA_SM_matrix= NODA_SM_matrix, WNODA_DM_matrix=NODA_DM_matrix))
      
    }
    else {
      return(list(WNODArow=NODArow,WNODAcol=NODAcol, WNODAmatrix=NODAmatrix))}
  }
  
  ### Using functions ####
  if(decreasing=="abund"){
    return(weightednoda(x,constraints))
  }
  if (decreasing=="fill"){
    if (weighted==FALSE){
      return(unweightednodf(x,constraints))
    }
    if (weighted==TRUE){
      return(weightednodf(x,constraints))
    }
  }
}

module2constraints <- function(mod){
  # helper function to extract the module to which a species belongs from a computeModule-object
  # Note that row 1 and columns 1 and 2 of mod@modules are for book-keeping
  # returns a single vector, with first rows and then columns
  # This vector contains a number for each module, in the position of the species that belong to this module; i.e. it starts with 4, 3, 1, ... indicating that the first species belongs to module 4, the second to module 3, the third to module 1, and so forth.
    apply(mod@modules[-1, -c(1,2)], 2, function(x) which(x > 0))
}

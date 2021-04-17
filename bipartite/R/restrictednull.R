restrictednull <- function(web, Prior.Pij = "degreeprob", 
                           conditional.level = "modules", 
                           N=10, 
                           print.null = FALSE, 
                           allow.degeneration = FALSE, 
                           return.nonrm.species = FALSE, 
                           connectance = TRUE, 
                           byarea = FALSE, 
                           R.partitions = NULL, 
                           C.partitions = NULL){
  
  #### Vaznull algorithm of bipartite modified to run the restricted null model
  #### Authors: Gabriel M. Felix, Rafael B. P. Pinheiro, and Marco A. R. Mello
  
  M <- web
  # compute per-cell probabilities:
  Pij.Prob <- PosteriorProb(web = M, 
                            R.partitions = R.partitions, C.partitions = C.partitions, 
                            Prior.Pij = Prior.Pij,
                            conditional.level = conditional.level)
  
  ### Test of assumptions
  
  if (!is.matrix(M)){stop("M is not a matrix.")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated.")}
  
  if (!is.matrix(Pij.Prob)){stop("Pij is not a matrix.")}
  if (T %in% c(Pij.Prob < 0)){stop("Pij must contain only values >= 0.")}
  
  if (nrow(M) != nrow(Pij.Prob) | ncol(M) != ncol(Pij.Prob)){stop("Dimensions of M and Pij.Prob must be identical.")}
  
  if (byarea){
    if (is.null(C.partitions) | is.null(R.partitions)){stop("Partitions missing")}
    if (length(unique(c(length(R.partitions), nrow(M), nrow(Pij.Prob)))) != 1){stop("The number of elements of R.partition should be the same as the number of rows of M and Pij.prob.")}
    if (length(unique(c(length(C.partitions), ncol(M), ncol(Pij.Prob)))) != 1){stop("The number of elements of C.partition should be the same as the number of column of M and Pij.prob.")}
    if (!identical(sort(unique(R.partitions)), sort(unique(C.partitions)))){stop("The number and labels of modules in R.partition and C.partition must be the same")}
  }  
  
  if (N <= 0 | !is.numeric(N)) {stop("Number of nulls, N, must be > 0")}
  if (!is.logical(connectance)){stop("Connectance must be logical (T or F)")}
  if (!is.logical(allow.degeneration)){stop("allow.degeneration must be logical (T or F)")}
  if (!is.logical(return.nonrm.species)){stop("return.nonrm.species must be logical (T or F)")}
  if (!is.logical(byarea)){stop("byarea must be logical (T or F)")}
  
  ### M dimensions
  nr <- dim(M)[1] # Number of rows
  nc <- dim(M)[2] # Number of collums
  
  ### Constructing a array with r rows, c columns and 2 slices. This array represents the matrix area structure
  if (byarea){
    Matrix.area <- array(0, dim = c(nr, nc, 2))
    for (rr in 1:nr){
      for (cc in 1:nc){
        Matrix.area[rr,cc,1] <- R.partitions[rr]
        Matrix.area[rr,cc,2] <- C.partitions[cc]
      }
    }
    
  } else { # byarea==F
    ## Assigning all rows and columns to the same partition in order to run the code bellow
    Matrix.area <- array(1, dim = c(nr, nc, 2))
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
    
  }
  
  ### Null model simulation
  
  NullMatrices <- list() # list where the null matrices will be storage
  length(NullMatrices) <- N #assigning the number of null matrices to be saved in NullMatrices 
  
  ## Drawing interaction in each null matrix 
  for (nn in 1:N){
    
    R.part <- sort(unique(as.vector(Matrix.area[,, 1])))
    C.part <- sort(unique(as.vector(Matrix.area[,, 2])))
    finalmat <- matrix(NA, nr, nc)
    
    for (R.p in R.part){
      for (C.p in C.part){
        
        M.a <- as.matrix(M[R.partitions == R.p, C.partitions == C.p])
        Pij.a <- Pij.Prob[R.partitions == R.p, C.partitions == C.p]
        
        r.a <- dim(M.a)[1]
        c.a <- dim(M.a)[2]
        
        P.a <- P1.a <- Pij.a
        finalmat.a <- matrix(0, r.a, c.a)
        
        if(allow.degeneration == FALSE & R.p == C.p){
          
          ## Ensuring that the dimensions of the null matrix will be the same of the original matrix
          
          D.int.finalmat.a <- 0 # The number of rows + columns occupied of the null matrix 
          while (D.int.finalmat.a < sum(dim(M.a))) { # While the dimensions of the null matrix was smaller then the original matrix, keep going
            sel <- sample(1:length(M.a), 1, prob = P.a) # Sample an cell of M.a with probability P.a
            finalmat.a[sel] <- 1 
            P.a[outer(as.numeric(rowSums(finalmat.a) > 0), as.numeric(colSums(finalmat.a)>0)) == 1] <- 0 #probability of cell with both rows and column non empty are zero
            D.int.finalmat.a <- sum(rowSums(finalmat.a) > 0) + sum(colSums(finalmat.a) > 0) # Setting the new number of dimensions occupied
          }
          # When the number of occupied dimensions of the null matrix was the same as the original matrix, continue
        }
        
        conn.remain <- sum(M.a > 0) - sum(finalmat.a > 0) # The number of cells remaining to be occupied to mantain the original connectance
        
        if (conn.remain > 0) {
          if (connectance == T){
            if (length(which(finalmat.a == 0)) == 1) {
              add <- which(finalmat.a == 0)
            } else {
              add <- sample(which(finalmat.a == 0), conn.remain, 
                            prob = P1.a[finalmat.a == 0], replace = FALSE)
            }
          }else {
            add <- sample(1:length(finalmat.a), conn.remain, 
                          prob = P1.a, replace = TRUE)
          }
          for (add1 in add){
            finalmat.a[add1] <- finalmat.a[add1] + 1
          }
        }
        
        ### Checking if there are still interactions to be drawn. If applicable, draw.
        int.remain <- (sum(M.a) - sum(finalmat.a))
        if (int.remain > 0) {
          if (length(which(finalmat.a > 0)) == 1) {
            add <- rep(which(finalmat.a > 0), int.remain)
          } else {
            add <- sample(which(finalmat.a > 0), int.remain, prob = P1.a[which(finalmat.a >0)], replace = TRUE)
          }
          finalmat.a[as.numeric(names(table(add)))] <- finalmat.a[as.numeric(names(table(add)))] + (table(add)) 
        }
        
        finalmat[R.partitions == R.p, C.partitions == C.p] <- finalmat.a
      }
    }
    
    # Saving outputs
    R2keep <- which(rowSums(finalmat) != 0)
    C2keep <- which(colSums(finalmat) != 0)
    finalmat2 <- finalmat[R2keep,C2keep]
    if (return.nonrm.species){
      NullMatrices[[nn]] = list(NullMatrix = finalmat2, R.Kept = R2keep, C.Kept = C2keep)
    } else { # return.nonrm.species==F:
      NullMatrices[[nn]] = finalmat2
    }
    if (print.null){print(nn)}
  }
  return(NullMatrices = NullMatrices)
}

# helper function, not exported or documented in detail:

PosteriorProb <- function(web, R.partitions = NULL, C.partitions = NULL, Prior.Pij = "degreeprob", conditional.level="modules"){
  #### Compute probabilities of interaction in a network with a modular structure.
  #### Authors: Gabriel M. Felix, Rafael B. P. Pinheiro, and Marco A. R. Mello.
  #### Ecological Synthesis Lab (SintECO)
  
  M <- web
  
  if (conditional.level == "matrix"){
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
  }
  
  if (conditional.level == "modules" & (is.null(R.partitions) | is.null(C.partitions))) stop("When using conditional.level='modules', R- and C-partitions must be given, based on modularity. See example.")
  
  # Test of assumptions:
  if (!is.matrix(M)){stop("M is not a matrix.")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated. There are rows and/or columns without interactions in the matrix. Remove them before proceeding.")}
  if (!is.numeric(R.partitions) | !is.numeric(C.partitions)) {stop("Partitions are not numeric.")}
  if (length(R.partitions) != nrow(M) | length(C.partitions) != ncol(M)) {stop("Partitions and matrix dimensions have different sizes.")}
  if (!(conditional.level %in% c("matrix","modules","areas"))) {stop("conditional.level should be 'matrix','modules' or 'areas'.")}
  if (!(Prior.Pij %in% c("degreeprob","equiprobable","degreeprob.byarea"))) {stop("Pij.probs should be 'equiprobable' or 'degreeprob' or 'degreeprob.byarea.")}
  
  # M dimensions
  nr <- dim(M)[1] # Number of rows
  nc <- dim(M)[2] # Number of columns
  #array()
  
  # Making an array with r rows, c columns, and 3 slices. This array represents the modular structure.
  # The first slice informs if a given cell M(rc) is within (1) or outside (0) a module.
  # The second slice informs to which module the species in the row (r) of a given cell M(rc) belongs.
  # The third slice informs to which module the species in the column (c) of a given cell M(rc) belongs .
  
  Matrix.mod <- array(0, dim = c(nr, nc, 3))
  for (rr in 1:nr){
    for (cc in 1:nc){
      Matrix.mod[rr, cc, 1] <- ifelse(R.partitions[rr] == C.partitions[cc], 1, 0)
      Matrix.mod[rr, cc, 2] <- R.partitions[rr]
      Matrix.mod[rr, cc, 3] <- C.partitions[cc]
    }
  }
  
  # Defining a priori Pij probabilities.
  if (Prior.Pij == "equiprobable"){
    Pi <- rep(1 / nr, times=nr) 
    Pj <- rep(1 / nc, times=nc)
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }
  if (Prior.Pij == "degreeprob"){
    Pi <- rowSums(M) / sum(rowSums(M)) 
    Pj <- colSums(M) / sum(colSums(M))
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }
  if(Prior.Pij == "degreeprob.byarea"){
    Prior.Pij.species <- M
    RMod <- sort(unique(R.partitions))
    CMod <- sort(unique(C.partitions))
    for (rr in RMod){
      for (cc in CMod){
        M.rr.cc <- matrix(M[R.partitions == rr,C.partitions == cc], sum(1*(R.partitions == rr)), sum(1*(C.partitions == cc)))
        Pi.rr.cc <- rowSums(M.rr.cc) / sum(rowSums(M.rr.cc)) 
        Pj.rr.cc <- colSums(M.rr.cc) / sum(colSums(M.rr.cc))
        Prior.Pij.species[R.partitions == rr, C.partitions == cc] <- tcrossprod(Pi.rr.cc, Pj.rr.cc)
      }
    }
  }
  
  # Defining conditional probabilities by area based on species degrees and connectance by area.
  
  if (conditional.level == "matrix"){
    Post.Pij <- Prior.Pij.species
  } else {
    Prior.Pij.area <- matrix(NA, nr, nc)
    Cond.Pij.area <- matrix(NA, nr, nc)
    if (conditional.level == "modules"){
      WMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 1])
      OMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 0])
      Prior.Pij.area[Matrix.mod[,,1] == 1] <- WMod.prior
      Prior.Pij.area[Matrix.mod[,,1] == 0] <- OMod.prior
      
      WMod.cond <- sum(M[Matrix.mod[,,1] == 1]) / sum(M)
      OMod.cond <- sum(M[Matrix.mod[,,1] == 0]) / sum(M)
      Cond.Pij.area[Matrix.mod[,,1] == 1] <- WMod.cond
      Cond.Pij.area[Matrix.mod[,,1] == 0] <- OMod.cond
    }
    if (conditional.level == "areas"){
      RMod <- sort(unique(R.partitions))
      CMod <- sort(unique(C.partitions))
      for (rr in RMod){
        for (cc in CMod){
          WArea.prior <- sum(Prior.Pij.species[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc])
          Prior.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.prior
          
          WArea.cond <- sum(M[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc]) / sum(M)
          Cond.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.cond
        }
      } 
    }
    
    # Adjusting the prior Pij prob by conditional probabilities. 
    
    Post.Pij <- Prior.Pij.species * (Cond.Pij.area / Prior.Pij.area)
  }
  
  return(Post.Pij = Post.Pij)
}


# #Compute modularity
# Mod <- computeModules(Safariland)
# 
# #Recover the partitions
# Part <- module2constraints(Mod)
# row.Part <- Part[1:nrow(Safariland)]
# col.Part <- Part[(nrow(Safariland)+1):(nrow(Safariland)+ncol(Safariland))]
# 
# 
# #Pij <- PosteriorProb(M = Safariland, 
# #                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
# #                     Prior.Pij = "degreeprob", #Choose the null model
# #                     conditional.level = "modules") #Choose the kind of constraints
# 
# #Check what those probabilities look like
# #Pij
# 
# #Generate randomized networks with the null model of your choice, considering the interaction probabilities calculated before. 
# nulls <- restrictednull(web = Safariland, 
#                        R.partitions = row.Part, C.partitions = col.Part)
# # nulls <- restrictednull(web = Safariland) # should return error!
# 
# #Calculate the same nestedness metric for all randomized networks
# null <- sapply(nulls, nest.smdm, constraints = Part, weighted = TRUE, decreasing = "abund")
# (WNODA.null <- unlist(null[1,])) # WNODArow
# (WNODAsm.null <- unlist(null[2,])) # WNODAcol
# (WNODAdm.null <- unlist(null[3,])) # WNODAmatrix
# # observed values:
# nest.smdm(Safariland, weighted = TRUE, decreasing = "abund")

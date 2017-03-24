# LPA_wb_plus.R
# Label propagation algorithm for weighted bipartite networks that finds modularity.
# Contains the LPAwb+ and the DIRTLPAwb+ algorithms
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

###############################################

LPA_wb_plus <- function(MATRIX, initialmoduleguess=NA) {
  
  #Make sure the smallest matrix dimension represent the red labels by making them the rows (If matrix is transposed here, will be transposed back at the end)
  flipped = 0
  if(dim(MATRIX)[1] > dim(MATRIX)[2]) {
    MATRIX = t(MATRIX)
    flipped = 1
  }
  
  Matsum = sum(MATRIX)
  col_marginals = colSums(MATRIX)
  row_marginals = rowSums(MATRIX)
  BMatrix = BarbersMatrix(MATRIX)
  
  #initiliase labels
  bluelabels = rep(NA,dim(MATRIX)[2])
  
  if (is.na(initialmoduleguess))
    redlabels = 1:dim(MATRIX)[1]
  else
    redlabels = sample(1:(initialmoduleguess+1),dim(MATRIX)[1],replace=TRUE)
  
  #Run Phase 1: Locally update lables to maximise Qb
  outlist = StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  #Run Phase 2: Connect divisions from top-down if it improves Qb, then run
  #phase 1 again. Repeat till Qb cannot be improved.
  outlist = StageTwo_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,Qb_now)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  if(flipped==1) { #If matrix was flipped, swap row and column labels
    holder = redlabels
    redlabels = bluelabels
    bluelabels = holder
  }
  
  return(list(Row_labels=redlabels, Col_labels=bluelabels, modularity=Qb_now))
}

###############################################

DIRT_LPA_wb_plus <- function(MATRIX,mini=4,reps=10) {
  A=LPA_wb_plus(MATRIX)
  
  mods=length(unique(A[[1]]))
  
  if((mods-mini) > 0) {
    for(aa in mini:mods) {
      for(bb in 1:reps) {
        B=LPA_wb_plus(MATRIX,aa)
        if(B[[3]]>A[[3]])
          A=B
      }
    }
  }
  
  return(list(Row_labels=A[[1]], Col_labels=A[[2]], modularity=A[[3]]))
}




###############################################

BarbersMatrix <- function(MATRIX) {
  return(MATRIX - (cbind(rowSums(MATRIX))%*%rbind(colSums(MATRIX)))/sum(MATRIX))
}

###############################################

WEIGHTEDMODULARITY <- function(BMatrix,Matsum,redlabels,bluelabels) {
  #see equation 8
  holdsum = 0
  
  for (rr in 1:length(redlabels)) {
    for (cc in 1:length(bluelabels)) {
      kroneckerdelta = redlabels[rr] == bluelabels[cc]
      holdsum = holdsum + BMatrix[rr,cc] * kroneckerdelta
    }
  }
  return(holdsum/Matsum)
}

###############################################
TRACE <- function(MATRIX) { return(sum(diag(MATRIX))) }

###############################################

WEIGHTEDMODULARITY2 <- function(BMatrix,Matsum,redlabels,bluelabels) {
  #see equation 9
  UNIred = unique(redlabels)
  Lred = length(UNIred)
  UNIblu = unique(bluelabels)
  Lblu = length(UNIblu)
  LABELMAT1 = matrix(0,Lred,length(redlabels))
  LABELMAT2 = matrix(0,length(bluelabels),Lblu)
  
  for(aa in 1:length(redlabels))
    LABELMAT1[which(UNIred == redlabels[aa]),aa] = 1
  
  for(aa in 1:length(bluelabels))
    LABELMAT2[aa,which(UNIblu == bluelabels[aa])] = 1
  
  return(TRACE(LABELMAT1 %*% BMatrix %*% LABELMAT2)/Matsum)
  
}

###############################################

DIVISION <- function(redlabels,bluelabels) {
  
  divisionsFound <- intersect(redlabels,bluelabels)
  
  return(divisionsFound)
}

###############################################

StageTwo_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels, Qb_now) {
  
  divisionsFound = DIVISION(redlabels,bluelabels)
  NUMdiv = length(divisionsFound)
  IterateFlag = 1
  while(IterateFlag == 1) {
    CombinedDivisionsThisTime = 0
    if(NUMdiv > 1) {
      for(div1check in 1:(NUMdiv-1)) {
        Mod1 = divisionsFound[div1check]
        for(div2check in (div1check+1):NUMdiv) {
          CHECK_RED = redlabels
          CHECK_RED[which(redlabels==Mod1)] = divisionsFound[div2check]
          CHECK_BLUE = bluelabels
          CHECK_BLUE[which(bluelabels==Mod1)] = divisionsFound[div2check]
          
          
          QQ = WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED,CHECK_BLUE)			
          if(QQ > Qb_now) { #If a division will improve modularity - find the best way to do this
            FoundBetter = 0
            for(aa in 1:NUMdiv) {
              CHECK_RED2 = redlabels
              CHECK_RED2[which(redlabels==divisionsFound[aa])] = Mod1
              CHECK_BLUE2 = bluelabels
              CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = Mod1
              if(WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                FoundBetter = 1
              }
              CHECK_RED2 = redlabels
              CHECK_RED2[which(redlabels==divisionsFound[aa])] = divisionsFound[div2check]
              CHECK_BLUE2 = bluelabels
              CHECK_BLUE2[which(bluelabels==divisionsFound[aa])] = divisionsFound[div2check]
              if(WEIGHTEDMODULARITY2(BMatrix,Matsum,CHECK_RED2,CHECK_BLUE2) > QQ) {
                FoundBetter = 1
              }
            }
            if(FoundBetter == 0) { #If no better configuration found - JOIN.
              redlabels = CHECK_RED
              bluelabels = CHECK_BLUE
              CombinedDivisionsThisTime = CombinedDivisionsThisTime + 1
            }
          }
        }
      }
      if(CombinedDivisionsThisTime == 0) {#If no divisions were joined move on
        IterateFlag = 0
      }
    }
    else {
      IterateFlag = 0		
    }
    
    outlist=StageOne_LPAwbdash(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels)  ##
    redlabels = outlist[[1]]
    bluelabels = outlist[[2]]
    Qb_now = outlist[[3]]
    divisionsFound = DIVISION(redlabels,bluelabels)
    NUMdiv = length(divisionsFound)
  }
  
  return(list(redlabels, bluelabels, Qb_now))
}


###############################################

StageOne_LPAwbdash <- function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels) {
  #Create storage containers for total marginals attached to each red(row)
  #label and blue(column) label
  
  BLUELABELLENGTH=length(bluelabels)
  REDLABELLENGTH=length(redlabels)
  TotalRedDegrees = rep(NA,max(redlabels))
  TotalBlueDegrees = rep(NA,max(BLUELABELLENGTH,REDLABELLENGTH))
  
  #Fill up these containers according to current labels
  #Red
  for(aa in 1:REDLABELLENGTH) {
    if(is.na(TotalRedDegrees[redlabels[aa]])) {
      TotalRedDegrees[redlabels[aa]] = row_marginals[aa]
    }
    else {
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa]
    }
  }
  
  #Blue
  if(sum(is.na(bluelabels)) != BLUELABELLENGTH) { #occurs first time through as blue nodes unlabelled
    for(bb in 1:BLUELABELLENGTH) {
      if(is.na(TotalBlueDegrees[bluelabels[bb]])) {
        TotalBlueDegrees[bluelabels[bb]] = col_marginals[bb]
      }
      else {
        TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb]
      }
    }
  }
  else {
    TotalBlueDegrees = rep(0,max(BLUELABELLENGTH,REDLABELLENGTH))
  }
  
  
  #locally maximise modularity!
  outlist = LOCALMAXIMISATION(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees)
  redlabels = outlist[[1]]
  bluelabels = outlist[[2]]
  Qb_now = outlist[[3]]
  
  return(list(redlabels, bluelabels, Qb_now))
  
}



###############################################

LOCALMAXIMISATION <-  function(row_marginals,col_marginals,MATRIX,BMatrix,Matsum,redlabels,bluelabels,TotalRedDegrees,TotalBlueDegrees) {
  
  #Find score for current partition
  QbAfter = WEIGHTEDMODULARITY2(BMatrix,Matsum,redlabels,bluelabels)
  
  if(is.na(QbAfter)) { QbAfter = -999 }
  
  IterateFlag = 1
  while(IterateFlag == 1) {
    #Save old information
    QbBefore = QbAfter
    old_redlabels = redlabels
    old_bluelabels = bluelabels
    old_TRD = TotalRedDegrees
    old_TBD = TotalBlueDegrees
    
    #Update Blue Nodes using red node information (see equation 10)
    bluelabelchoices = unique(redlabels)
    
    for(bb in 1:length(bluelabels)) {
      if(is.na(bluelabels[bb]) == FALSE) {
        TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] - col_marginals[bb] 
      }
      changebluelabeltest = rep(NA,length(bluelabelchoices))
      
      for(ww in 1:length(bluelabelchoices)) {
        changebluelabeltest[ww] = sum( (redlabels == bluelabelchoices[ww]) * MATRIX[,bb])  -  col_marginals[bb]*TotalRedDegrees[bluelabelchoices[ww]]/Matsum 
      }
      
      #assign new label based on maximisation of above condition  
      
      labels = which(changebluelabeltest == max(changebluelabeltest,na.rm =TRUE))
      newlabelindex = labels[sample(1:length(labels),1)]
      bluelabels[bb] = bluelabelchoices[newlabelindex[1]]
      if(bluelabels[bb] > length(TotalBlueDegrees)) {
        TotalBlueDegrees[bluelabels[bb]] = 0
      }
      
      #Update total marginals on new labelling
      TotalBlueDegrees[bluelabels[bb]] = TotalBlueDegrees[bluelabels[bb]] + col_marginals[bb]
    }
    
    #Now update red node labels based on blue node information (see equation 10)
    redlabelchoices = unique(bluelabels)
    
    
    for(aa in 1:length(redlabels)) {
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] - row_marginals[aa]
      changeredlabeltest = rep(NA,length(redlabelchoices))
      
      for(ww in 1:length(redlabelchoices)) {
        changeredlabeltest[ww] = sum( (bluelabels == redlabelchoices[ww]) * MATRIX[aa,])  -  row_marginals[aa]*TotalBlueDegrees[redlabelchoices[ww]]/Matsum  
      }
      
      #assign new label based on maximisation of above condition
      labels = which(changeredlabeltest == max(changeredlabeltest,na.rm = TRUE))
      newlabelindex = labels[sample(1:length(labels),1)]
      redlabels[aa] = redlabelchoices[newlabelindex[1]]
      
      if(redlabels[aa] > length(TotalRedDegrees)) {
        TotalRedDegrees[redlabels[aa]] = 0
      }
      TotalRedDegrees[redlabels[aa]] = TotalRedDegrees[redlabels[aa]] + row_marginals[aa]
    }
    
    
    #Find the new modularity score based on node label updates.
    QbAfter = WEIGHTEDMODULARITY(BMatrix,Matsum,redlabels,bluelabels)
    
    #If this modularity is not as good as previous stop iterating and
    #use that previous best information
    
    if(QbAfter <= QbBefore) {
      redlabels = old_redlabels
      bluelabels = old_bluelabels
      TotalRedDegrees = old_TRD
      TotalBlueDegrees = old_TBD
      IterateFlag = 0
    }
    
  }
  
  Qb_now = QbAfter
  
  
  return(list(redlabels, bluelabels, Qb_now))
}


# convert2moduleWeb.R
# Use output from LPAwb+ to create a moduleWeb class object that can interface with the bipartite library in R
# Author :  Stephen Beckett ( https://github.com/sjbeckett/weighted-modularity-LPAwbPLUS )
# MIT License

#Specifically this can be used in order to make use of the plotModuleWeb function in bipartite to visualise detected modular structure
convert2moduleWeb <- function(MATRIX, MODINFO){
  #NET - network biadjacency matrix
  #MODINFO - output from running LPA_wb_plus or Exhaustive_LPA_wb_plus , contains
  
  #if names are NULL - assign names to network
  if (is.null(rownames(MATRIX)))
    rownames(MATRIX) = 1:dim(MATRIX)[1]
  
  if (is.null(colnames(MATRIX)))
    colnames(MATRIX) = 1:dim(MATRIX)[2]
  
  #module sorting
  ROW_IX <- order(MODINFO$Row_labels)
  COL_IX <- order(MODINFO$Col_labels)
  ROWS <- MODINFO$Row_labels[ROW_IX]
  COLS <- MODINFO$Col_labels[COL_IX]
  MODS <- unique(ROWS)
  LMod <- length(MODS)
  
  #Creating the modules matrix
  #as LPAwb+ does not search at multiple depths can fix the first two columns
  Col1 <- c(0,rep(1,LMod)) #depth
  Col2 <- c(1,rep(0,LMod-1),1) # module markers
  Vals <- 1:(length(ROWS)+length(COLS))
  #create matrix store, with LMod rows - where within module indices are shown and others are zero.
  store <- Vals
  for(bb in MODS) {
    Ix1 <- MODINFO$Row_labels!=bb
    Ix2 <- MODINFO$Col_labels!=bb
    Ix <- c(Ix1,Ix2)
    New <- Vals
    New[Ix] <- 0
    store <- rbind(store,New)
  }
  
  modules <- cbind(Col1,Col2,store)
  rownames(modules) <- NULL
  colnames(modules) <- NULL
  
  #assigning values to moduleWeb object	
  
  ModwebObj <- new("moduleWeb")
  ModwebObj@originalWeb <- MATRIX
  ModwebObj@moduleWeb <- MATRIX[ROW_IX,COL_IX]
  ModwebObj@orderA <- ROW_IX
  ModwebObj@orderB <- COL_IX
  ModwebObj@modules <- modules
  ModwebObj@likelihood <- MODINFO$modularity
  
  return(ModwebObj)
}

## 
# resLPA <- LPA_wb_plus(Safariland) # find labels and weighted modularity using LPAwb+
# computeModules(Safariland) # higher value (0.429 vs. 0.427)
# resDirtLPA <- DIRT_LPA_wb_plus(Safariland) # find labels and weighted modularity using DIRTLPAwb+
# better than computeModules (0.430)
# mod <- convert2moduleWeb(Safariland, resLPA)
# plotModuleWeb(mod)

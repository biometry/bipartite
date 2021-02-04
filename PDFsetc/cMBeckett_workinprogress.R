# This function actually prepares the recursive computation of the modules and returns an object of class "moduleWeb"
cMBeckett = function(web, depth, nrOfModule, ytop, xleft, ybottom, xright, prev_orderA, prev_orderB, modules, deepCompute, delete, steps, tolerance, experimental) {
  
  result = list();
  
  webName = paste("web-", depth, "-", nrOfModule, sep="");
  
  #### Start change for Beckett ####
  
  # write web to file
  #web2edges(web[ytop:ybottom, xleft:xright], webName=webName);
  #argv = c("identifyModules", "-filename", paste(webName, ".pairs", sep=""), "-steps", round(steps), "-tolerance", as.double(tolerance)); # round instead of "as.integer(steps)" because the latter works only up to 1E9!!
  #if(experimental) {
  #  argv = append(argv, c("-method", "Strauss"));
  #}
  #argv = as.character(argv);
  #argc = as.integer(length(argv));
  
  #### Here DormannStrauss calls the C-dynamic link library, returning files!
  
  #.C("identifyModules", argc, argv, PACKAGE="bipartite");
  ## because of unresolved issues, the dll needs to be unloaded/reloaded after each run.
  ## this does not work under Linux (see ?dyn.unload)
  
  ##	LIBS <- .dynLibs()
  ##	bipLIB <- which(unlist(sapply(LIBS, function(x) x[1])) == "bipartite")
  ##	IMpath <- 	LIBS[[bipLIB]][[2]] # absolute path on the system to the dll!!!
  ###library.dynam.unload("bipartite", libpath=IMpath)
  ###library.dynam("bipartite", package="bipartite", lib.loc=find.package("bipartite"))
  ##	dyn.unload(IMpath)
  ##	dyn.load(IMpath)
  
  #### Here the files are read in again!
  
  ## read in data from result files
  # data = readModuleData(webName, deleteOriginalFiles=delete); #delete.edgesOriginalFiles=delete);
  
  #orderAFile	= data[[1]];
  #orderBFile	= data[[2]];
  #modulesFile	= data[[3]];
  #likelihood	= as.numeric(data[[4]]);
  
  
  
  mod <- if (forceLPA) LPA_wb_plus(web) else  DIRT_LPA_wb_plus(web)
  data <- convert2moduleWeb(web,  mod)  
  
  orderAFile <- data@orderA
  orderBFile <- data@orderB
  modulesFile <- data@modules
  likelihood <- data@likelihood
  
  ##### End Change for Beckett #####
  
  n_a		= dim(web)[1];
  n_b		= dim(web)[2];
  n		= n_a + n_b;
  
  # Permutation of the graph and therefore actualization of the permutation vectors is necessary
  # if more than one module are suggested
  if (likelihood >= 0 && nrow(modulesFile) > 1) {
    
    # Actualization of the permutation vectors
    tempA				= c(1:dim(web)[1]);
    tempB				= c(1:dim(web)[2]);
    tempA[ytop:ybottom]	= tempA[ytop:ybottom][orderAFile];
    tempB[xleft:xright]	= tempB[xleft:xright][-(orderBFile-length(orderBFile))+1]; # CFD
    orderAFile			= prev_orderA[tempA];
    orderBFile			= prev_orderB[tempB];
    
    result[[1]] = as.matrix(web[tempA, tempB, drop=FALSE])
    result[[2]] = orderAFile;
    result[[3]] = orderBFile;
    result[[4]] = NA;
    
    # The matrix M containing the information about the identified modules is formatted in the following way:
    # Each row i contains the information about one certain module while
    # M[i,1] represents the nesting depth of the module,
    # M[i,2] = 1 iff the module is the last one to plot within its nesting module and 0 else,
    # and the M[i,j] = j-offset_M iff node j-offset_M is part of the module and 0 else.
    offset_M			= 2;
    offset_modulesFile	= 1;	
    nrOfModules			= nrow(modulesFile);
    M					= matrix(0, nrOfModules, (sum(dim(web))+offset_M));
    
    if (experimental) {
      
      M[,1]	= modulesFile[,1];
      M[1,2]	= 1;
      
      for(i in nrOfModules:1) {
        
        if(i > 1) {
          indicesOfNextRowsWithSameDepth = rev(which(M[1:(i - 1),1] == M[i,1]));
          
          if(	length(indicesOfNextRowsWithSameDepth) == 0
              ||
              (indicesOfNextRowsWithSameDepth[1] + 1 < i && length(which(M[(indicesOfNextRowsWithSameDepth[1] + 1):(i - 1),1] < M[i,1])) > 0)) {
            M[i,2] = 1;
          }
        }
        
        M[i, (1 + offset_M):(ncol(modulesFile)+1)] = modulesFile[i, (1 + offset_modulesFile):ncol(modulesFile)];
      }
      
      M_temp		= matrix(0, nrow(M), ncol(M));
      maxDepth	= max(M[,1]);
      rowCounter	= 1;
      
      for(i in 0:maxDepth) {
        indicesOfRowsWithSameDepth		= rev(which(M[,1] == i));
        nrOfIndicesOfRowsWithSameDepth	= length(indicesOfRowsWithSameDepth);
        M_temp[rowCounter:(rowCounter + nrOfIndicesOfRowsWithSameDepth - 1),] = M[indicesOfRowsWithSameDepth,];
        rowCounter = rowCounter + nrOfIndicesOfRowsWithSameDepth;
      }
      
      M = M_temp;
    }
    else {
      
      M[, 1]								= depth;
      M[nrOfModules, offset_M]			= 1;
      
      colvals = append(prev_orderA[ytop:ybottom], prev_orderB[xleft:xright]+n_a);
      rowvals = which(modulesFile[, offset_M:ncol(modulesFile)] > 0) %% nrow(modulesFile);
      rowvals[rowvals == 0] = nrow(modulesFile);
      if (length(colvals) == length(rowvals)){
        for(i in 1:length(colvals)) {
          M[rowvals[i], colvals[i]+offset_M] = colvals[i];
        }
      }
    }
    
    modulesFile = M;
    
    if(experimental) {
      result[[4]]	= M;
    }
    else {
      result[[4]]	= rbind(modules, M);
    }
    
    # Computation of potential modules nested within the ones found until now
    if(deepCompute) {
      order = append(orderAFile, (orderBFile+n_a));
      
      # Apply the recursive function cM(...) to each module
      for (i in 1:nrow(modulesFile)) {
        mod = modulesFile[i, (offset_M+1):ncol(modulesFile)]
        mod = mod[order]
        
        # Calculate the coordinates of the part of the web we are looking at
        j = n_a + 1;
        while (mod[j] == 0) { j = j+1 }			# calculate x-coordinate of left lower corner of module border
        xleft_new = j - n_a;
        
        j = n_a;
        while(mod[j] == 0) { j = j-1; }			# calculate y-coordinate of left lower corner of module border
        ybottom_new = j;
        
        j = n;
        while(mod[j] == 0) { j = j-1; }			# calculate x-coordinate of right upper corner of module border
        xright_new = j - n_a;
        
        j = 1;
        while(mod[j] == 0) { j = j+1; }			# calculate y-coordinate of right upper corner of module border
        ytop_new = j;
        
        # An invocation of cM(...) is necessary only if there is the possibility to find more than one submodule the current module consists of
        if((ybottom_new - ytop_new)+1 > 1 && (xright_new - xleft_new)+1 > 1) {
          
          print(paste("Recursive invocation (depth: ", depth+1, ", module nr. ", i, ")", sep=""));
          
          result = cMBeckett(web=result[[1]], depth+1, nrOfModule=i, ytop=ytop_new, xleft=xleft_new, 
                      ybottom=ybottom_new, xright=xright_new, prev_orderA=result[[2]], prev_orderB=result[[3]], 
                      modules=result[[4]], deepCompute, delete, steps, tolerance, experimental);
          
          # Make sure that all computed modules have modularity >= 0
          if(result[[5]] < 0) {
            likelihood = result[[5]];
          }
        }
      }
    }
  }
  else {
    result[[1]] = web;
    result[[2]] = prev_orderA;
    result[[3]] = prev_orderB;
    result[[4]] = modules;
  }
  
  result[[5]] = likelihood;
  
  return(result)
} 


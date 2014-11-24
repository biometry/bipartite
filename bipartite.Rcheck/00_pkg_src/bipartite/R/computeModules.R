setClass("moduleWeb", representation(originalWeb="matrix", moduleWeb="matrix", orderA="vector", orderB="vector", modules="matrix", likelihood="numeric"));

computeModules = function(web, deep=FALSE, deleteOriginalFiles=TRUE, steps=1000000, tolerance=1e-10, experimental=FALSE) {

	if (deep && experimental) {
		print("ERROR: Flag deep has to be FALSE if flag experimental is set to TRUE.");
		NULL;
	} else {

		web <- as.matrix(empty(web, count=TRUE)) # to get rid of empty columns and rows
		if (any(attr(web, "empty")) > 0) warning("Some empty columns or rows were deleted.")
		
		cat("This function is SLOW! \nMonitor your system's activity; depending on steps-setting this may take minutes to hours!\n\n")
		
		result		<-  cM(web, depth=1, nrOfModule=1, ytop=1, xleft=1, ybottom=dim(web)[1], 
						  xright=dim(web)[2], prev_orderA=c(1:dim(web)[1]), prev_orderB=c(1:dim(web)[2]), 
						  modules=matrix(c(0, 1, c(1:sum(dim(web)))), 1), 
						  deepCompute=deep, delete=deleteOriginalFiles, steps=steps, tolerance=tolerance, experimental=experimental) #try()
	#	if (inherits(result, "try-error")) {
	#		print(result)
	#		break;
	#	}
		result[[4]]	= result[[4]][order(result[[4]][,1]),];

		# Make sure the modularity is non-negative
		if (result[[5]] >= 0) {
			new("moduleWeb", originalWeb=web, moduleWeb=as.matrix(result[[1]]), orderA=result[[2]], orderB=result[[3]], modules=result[[4]], likelihood=result[[5]]);
		}
		else {
			if(deep) {
				print("ERROR: At least one of the computed modularity values is smaller than 0. Please increase the number of steps in order to reach a valid result.");
			}
			else {
				print("ERROR: The computed modularity value is smaller than 0. Please increase the number of steps in order to reach a valid result.");
			}
			NULL;
		}
	}
}


# This function actually prepares the recursive computation of the modules and returns an object of class "moduleWeb"
cM = function(web, depth, nrOfModule, ytop, xleft, ybottom, xright, prev_orderA, prev_orderB, modules, deepCompute, delete, steps, tolerance, experimental) {

	result = list();

	webName = paste("web-", depth, "-", nrOfModule, sep="");

	# write web to file
	web2edges(web[ytop:ybottom, xleft:xright], webName=webName);
	argv = c("identifyModules", "-filename", paste(webName, ".pairs", sep=""), "-steps", round(steps), "-tolerance", as.double(tolerance)); # round instead of "as.integer(steps)" because the latter works only up to 1E9!!
	if(experimental) {
		argv = append(argv, c("-method", "Strauss"));
	}
	argv = as.character(argv);
	argc = as.integer(length(argv));

	.C("identifyModules", argc, argv, PACKAGE="bipartite");
## because of unresolved issues, the dll needs to be unloaded/reloaded after each run.
## this does not work under Linux (see ?dyn.unload)

#	LIBS <- .dynLibs()
#	bipLIB <- which(unlist(sapply(LIBS, function(x) x[1])) == "bipartite")
#	IMpath <- 	LIBS[[bipLIB]][[2]] # absolute path on the system to the dll!!!
	##library.dynam.unload("bipartite", libpath=IMpath)
	##library.dynam("bipartite", package="bipartite", lib.loc=find.package("bipartite"))
#	dyn.unload(IMpath)
#	dyn.load(IMpath)

	# read in data from result files
	data = readModuleData(webName, deleteOriginalFiles=delete); #delete.edgesOriginalFiles=delete);

	n_a		= dim(web)[1];
	n_b		= dim(web)[2];
	n		= n_a + n_b;

	orderAFile	= data[[1]];
	orderBFile	= data[[2]];
	modulesFile	= data[[3]];
	likelihood	= as.numeric(data[[4]]);
	
	# Permutation of the graph and therefore actualization of the permutation vectors is necessary
	# if more than one module are suggested
	if (likelihood >= 0 && nrow(modulesFile) > 1) {

		# Actualization of the permutation vectors
		tempA				= c(1:dim(web)[1]);
		tempB				= c(1:dim(web)[2]);
		tempA[ytop:ybottom]	= tempA[ytop:ybottom][orderAFile];
		tempB[xleft:xright]	= tempB[xleft:xright][orderBFile-length(orderAFile)];
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

		if(experimental) {
		
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
			if(length(colvals) == length(rowvals)) {
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
			for(i in 1:nrow(modulesFile)) {
				mod = modulesFile[i, (offset_M+1):ncol(modulesFile)];
				mod = mod[order];

				# Calculate the coordinates of the part of the web we are looking at
				j = n_a+1;
				while(mod[j] == 0) { j = j+1; }			# calculate x-coordinate of left lower corner of module border
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

					result = cM(web=result[[1]], depth+1, nrOfModule=i, ytop=ytop_new, xleft=xleft_new, 
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

	result

}


# This function reads in the data written down to the hard drive disk by the invoked C code and optionally triggers the deletion of the appropriate files
readModuleData = function(webName=NULL, deleteOriginalFiles=TRUE) {
	orderAFile = as.vector(as.matrix(read.table(paste(webName, ".ordA", sep=""), header=TRUE, sep="\t")));
	orderBFile = as.vector(as.matrix(read.table(paste(webName, ".ordB", sep=""), header=TRUE, sep="\t")));
	modulesFile = as.matrix(read.table(paste(webName, ".mod", sep=""), header=TRUE, sep="\t"));
	infoFile = strsplit(readLines(paste(webName, ".info", sep=""), n=-1)[26], ": ")[[1]][2]  # reads in only the likelihood!
	#read.table(paste(webName, ".info", sep=""), header=TRUE, sep="\t");

	if(deleteOriginalFiles) deleteModuleData(webName);

	result = list();
	result[[1]] = orderAFile;
	result[[2]] = orderBFile;
	result[[3]] = modulesFile;
	result[[4]] = infoFile;

	result;
}


# This function deletes the files written down to the hard drive disk by the invoked C code
deleteModuleData = function(webName=NULL) {
	unlink(paste(webName, ".pairs", sep=""));
	unlink(paste(webName, "-names.lut", sep=""));
	unlink(paste(webName, ".ordA", sep=""));
	unlink(paste(webName, ".ordB", sep=""));
	unlink(paste(webName, ".mod", sep=""));
	unlink(paste(webName, ".info", sep=""));
	unlink(paste(webName, ".den", sep=""));
	unlink(paste(webName, ".lut", sep=""));
}

# This function is ALSO in drawModules.R!!
# Auxiliary function checking whether the passed object is an object of class "moduleWeb" and contains correctly formatted information
isCorrectModuleWebObject = function(moduleWebObject) {

	if (!is(moduleWebObject, "moduleWeb")) {
		warning("Object of wrong class.");
		FALSE;
	}
	else if(dim(slot(moduleWebObject, "originalWeb")) == 0 ||  dim(slot(moduleWebObject, "moduleWeb")) != dim(slot(moduleWebObject, "originalWeb")) || dim(slot(moduleWebObject, "modules")) == 0) {
		warning("Object corrupt.");
		FALSE;
	}
	else if(min(slot(moduleWebObject, "originalWeb")) < 0 || min(slot(moduleWebObject, "moduleWeb")) < 0) {
		warning("entries of matrix have to be greater than or equal to 0.");
		FALSE;
	}
	else {
		TRUE;
	}
}
######----- ########

listModuleInformation = function(moduleWebObject) {

	if(isCorrectModuleWebObject(moduleWebObject)) {

		result	= list();

		web	= slot(moduleWebObject, "originalWeb");
		modules	= slot(moduleWebObject, "modules");

		n_a	= nrow(web);
		n_b	= ncol(web);

		offset	= 2;

		for(depth in unique(modules[,1])) {
			result[[depth+1]] = list();

			counter = 1;

			for(i in 1:nrow(modules)) {
				if(modules[i,1] == depth) {
					result[[depth+1]][[counter]]		= list();
					result[[depth+1]][[counter]][[1]]	= rownames(web)[modules[i,(offset+1):(n_a+offset)][modules[i,(offset+1):(n_a+offset)] > 0]];
					result[[depth+1]][[counter]][[2]]	= colnames(web)[(modules[i,(n_a+offset+1):(n_a+n_b+offset)][modules[i,(n_a+offset+1):(n_a+n_b+offset)] > 0])-n_a];

					counter = counter + 1;
				}
			}
		}

		result;
	}
}


printoutModuleInformation = function(moduleWebObject) {

	if(isCorrectModuleWebObject(moduleWebObject)) {

		modules	= slot(moduleWebObject, "modules");

		listOfModules = listModuleInformation(moduleWebObject);

		linebreak = "\n";

		a = linebreak;

		for(depth in unique(modules[,1])) {
			for(i in 1:length(listOfModules[[depth+1]])) {
				a = paste(a, "Depth: ", depth, linebreak, linebreak, "Nr of module: ", i, linebreak, linebreak, "Rownames: ", linebreak, sep="");
				for(j in 1:length(listOfModules[[depth+1]][[i]][[1]])) {
					a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[1]][j], sep=""), sep=linebreak);
				}
				a = paste(a, linebreak, linebreak, "Colnames: ", linebreak, sep="");
				for(j in 1:length(listOfModules[[depth+1]][[i]][[2]])) {
					a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[2]][j], sep=""), sep=linebreak);
				}
				a = paste(a, linebreak, linebreak, "__________________________", linebreak, linebreak, sep="");
			}
			a = paste(a, "__________________________", linebreak, linebreak, sep="");
		}

		cat(a);
	}
}

# web <- Safariland
# mod <- computeModules(web, steps=1E9)
# dir()
# read.delim("web-1-1.pairs") # works fine
# read.delim("web-1-1-names.lut") # works fine
# mod <- computeModules(web, steps=1E9) # does not work (on first run)

# Auxiliary function extending the actual web by empty rows and columns for plotting the rectangles representing modules
prepareWebForPlottingModules = function(moduleWebObject, fromDepth=0, upToDepth=-1) {

	web	= slot(moduleWebObject, "moduleWeb");
	orderA	= slot(moduleWebObject, "orderA");
	orderB	= slot(moduleWebObject, "orderB");
	modules	= slot(moduleWebObject, "modules");

	if(upToDepth < 0 || upToDepth > max(modules[,1])) upToDepth = max(modules[,1]);

	n_a		= nrow(web);
	n_b		= ncol(web);
	n		= n_a + n_b;

	order		= append(orderA, (orderB+n_a));

	offset		= 2;

	resultWeb = web;

	modules = modules[order(modules[,1], decreasing=TRUE),];

	for(i in 1:nrow(modules)) {
		if(fromDepth <= modules[i,1] && modules[i,1] <= upToDepth) {
			mod = modules[i,c((offset+1):ncol(modules))];
			mod = mod[order];

			j = 1
			while(mod[j] == 0) { j = j+1; }
			resultWeb = addEmptyRowToMatrix(resultWeb, which(rownames(resultWeb) == rownames(web)[j]));

			j = n_a+1;
			while(mod[j] == 0) { j = j+1; }
			j = j - n_a;
			resultWeb = addEmptyColToMatrix(resultWeb, which(colnames(resultWeb) == colnames(web)[j]));

			# if current module is last submodule within its nesting module add an empty row and an empty column
			# also at its right lower coordinate
			if(modules[i,2]) {
				j = n_a;
				while(mod[j] == 0) { j = j-1; }
				resultWeb = addEmptyRowToMatrix(resultWeb, which(rownames(resultWeb) == rownames(web)[j])+1);

				j = n;
				while(mod[j] == 0) { j = j-1; }
				resultWeb = addEmptyColToMatrix(resultWeb, which(colnames(resultWeb) == colnames(web)[j-n_a])+1);
			}
		}
	}

	resultWeb;

}


addEmptyRowToMatrix = function(matrix, x) {
	if (is.na(x)) { # added by CFD 23-Sep-2011
		warning("Error in addEmptyRowToMatrix: counter j is NA.")
		return(matrix)
	}
	if(x == 1) {
		rbind(0, matrix);
	}
	else if(x > 1 && x < nrow(matrix)) {
		name = rownames(matrix)[x-1];
		matrix = rbind(matrix[1:(x-1),], 0, matrix[x:nrow(matrix),]);
		rownames(matrix)[x-1] = name;
		matrix;
	}
	else if(x == nrow(matrix)) {
		name = rownames(matrix)[nrow(matrix)];
		matrix = rbind(matrix[1:(x-1),], 0, matrix[x:nrow(matrix),]);
		rownames(matrix)[nrow(matrix)] = name;
		matrix;
	}
	else if(x == nrow(matrix)+1) {
		rbind(matrix, 0);
	}
}


addEmptyColToMatrix = function(matrix, x) {
	if (is.na(x)) { # added by CFD 23-Sep-2011
		warning("Error in addEmptyRowToMatrix: counter j is NA.")
		return(matrix)
	}
	if(x == 1) {
		cbind(0, matrix);
	}
	else if(x > 1 && x < ncol(matrix)) {
		name = colnames(matrix)[x-1];
		matrix = cbind(matrix[,1:(x-1)], 0, matrix[,x:ncol(matrix)]);
		colnames(matrix)[x-1] = name;
		matrix;
	}
	else if(x == ncol(matrix)) {
		name = colnames(matrix)[ncol(matrix)];
		matrix = cbind(matrix[,1:(x-1)], 0, matrix[,x:ncol(matrix)]);
		colnames(matrix)[ncol(matrix)] = name;
		matrix;
	}
	else if(x == ncol(matrix)+1) {
		cbind(matrix, 0);
	}
}
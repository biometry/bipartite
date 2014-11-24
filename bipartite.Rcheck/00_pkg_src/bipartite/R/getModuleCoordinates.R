# This function computes and returns the coordinates of the corners of the rectangles representing modules
getModuleCoordinates = function(moduleWebObject, fromDepth=0, upToDepth=-1) {

	if (!is(moduleWebObject, "moduleWeb")) { warning("Object of wrong class."); 
  }	else {

		web		= slot(moduleWebObject, "originalWeb");
		moduleWeb	= slot(moduleWebObject, "moduleWeb");
		# moduleWeb contains empty rows and columns since function prepareWebForPlottingModules(...) has been applied to moduleWebObject before
		emptyModuleWeb	= empty(moduleWeb);
		orderA		= slot(moduleWebObject, "orderA");
		orderB		= slot(moduleWebObject, "orderB");
		modules		= slot(moduleWebObject, "modules");

		if(upToDepth < 0 || upToDepth > max(modules[,1])) upToDepth = max(modules[,1]);

		n_a		= nrow(web);
		n_b		= ncol(web);
		n		= n_a + n_b;

		foundModules	= matrix(0, nrow(modules), 5);

		order		= append(orderA, (orderB+n_a));
		modules		= modules[order(modules[,1], decreasing=TRUE),];

		offset		= 2;

		for(i in 1:nrow(modules)) {
			if (fromDepth <= modules[i,1] && modules[i,1] <= upToDepth) {
				mod = modules[i,c((offset+1):ncol(modules))];			# i.th row of "modules" without two first elements (containing information about depth and whether current module is the last submodule of its depth to plot within its nesting module)
				mod = mod[order];

				j = n_a+1;
				while (mod[j] == 0) { j = j+1; }					# calculate x-coordinate of left lower corner of module border
				x_left = which(colnames(moduleWeb) == colnames(emptyModuleWeb)[j - n_a]) - 1;

				j = n_a;
				while(mod[j] == 0) { j = j-1; }					# calculate y-coordinate of left lower corner of module border
				y_bottom = nrow(moduleWeb) - which(rownames(moduleWeb) == rownames(emptyModuleWeb)[j]);

				j = n;
				while(mod[j] == 0) { j = j-1; }					# calculate x-coordinate of right upper corner of module border
				x_right = which(colnames(moduleWeb) == colnames(emptyModuleWeb)[j - n_a]);

				j = 1;
				while(mod[j] == 0) { j = j+1; }					# calculate y-coordinate of right upper corner of module border
				y_top = nrow(moduleWeb) - which(rownames(moduleWeb) == rownames(emptyModuleWeb)[j]) + 1;

				# Since, possibly, there are already rectangles around some of the nodes comprised by the rectangle spanned by the computed coordinates,
				# we have to find the next empty or column, respectively
				while (moduleWeb[(-y_top + nrow(moduleWeb)), x_left] < -modules[i,1]) {
					x_left	= x_left - 1;
					y_top		= y_top + 1;
				}

				while (moduleWeb[(-y_bottom + nrow(moduleWeb)), x_right] < -modules[i,1]) {
					y_bottom	= y_bottom - 1;
					x_right		= x_right + 1;
				}

				# The rectangles representing the modules are plotted onto the empty rows and colums
				# centerOffset is used for making each line of the rectangles run exactly through the
				# middle of the appropriate empty row or column, respectively
				centerOffset = 0.5;

				foundModules[i,1] = modules[i,1];
				foundModules[i,2] = x_left - centerOffset;
				foundModules[i,3] = y_bottom - centerOffset;
				foundModules[i,4] = x_right + centerOffset;
				foundModules[i,5] = y_top + centerOffset;

				# Set entries of moduleWeb over which the current rectangle will be plotted to negative values
				# in order to be able to compute the correct coordinates of potential supermodules

				moduleWeb[(-y_top + nrow(moduleWeb)):(-y_bottom + nrow(moduleWeb)), x_right]		= -modules[i,1];
				moduleWeb[(-y_top + nrow(moduleWeb)):(-y_bottom + nrow(moduleWeb)), x_left]		= -modules[i,1];
				moduleWeb[(-y_bottom + nrow(moduleWeb)), x_left:x_right]				= -modules[i,1];
				moduleWeb[(-y_top + nrow(moduleWeb)), x_left:x_right]					= -modules[i,1];
			}
		}
	}

	foundModules;
}

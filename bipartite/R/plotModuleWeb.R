plotModuleWeb <- function (moduleWebObject, plotModules=TRUE, rank=FALSE, weighted=TRUE, displayAlabels=TRUE, 
	displayBlabels=TRUE, labsize=1, xlabel="", ylabel="", square.border="white", fromDepth=0, upToDepth=-1) {

	if (isCorrectModuleWebObject(moduleWebObject)) {

		if (plotModules) {
			web	<- prepareWebForPlottingModules(moduleWebObject, fromDepth, upToDepth)
			moduleWebObject	<- new("moduleWeb", originalWeb=slot(moduleWebObject, "originalWeb"), 
				moduleWeb=web, orderA=slot(moduleWebObject, "orderA"), orderB=slot(moduleWebObject, "orderB"), 
				modules=slot(moduleWebObject, "modules"))
		}
		else {
			web	<- slot(moduleWebObject, "originalWeb")
		}

		n_a	= nrow(web);
		n_b	= ncol(web);
		n	= n_a + n_b;

		maxMarginBottom = max(strwidth(colnames(web), units = "inches"));
		maxMarginLeft = max(strwidth(rownames(web), units = "inches"));
		par(mai = c(maxMarginBottom, maxMarginLeft, 0, 0));

		# clear plot screen
		plot(1, type = "n", axes = FALSE, xlim = c(0, n_b), ylim = c(0, n_a), asp = 1, xlab = xlabel, ylab = ylabel);

		# displayA labels
		Alabels = 0;
		if (displayAlabels && !is.null(rownames(web))) {
			for (i in 1:n_a) {
				s = rownames(web)[n_a - i + 1];
				if (!is.null(s)) { Alabels[i] = s; }
				else Alabels[i] = "";
       			}
			axis(2, (1:n_a) - 0.5, labels = Alabels, tick = FALSE, mgp = c(0, 0, 0), las = 2, cex.axis = labsize);
   		}

		# displayB labels
		Blabels = 0
		if (displayBlabels && !is.null(colnames(web))) {
			for (j in 1:n_b) {
				s = colnames(web)[j]
				if (!is.null(s)) { Blabels[j] = s; }
				else Blabels[i] = "";
			}
			axis(1, (1:n_b) - 0.5, labels = Blabels, tick = FALSE, mgp = c(0, 0, 0), las = 2, cex.axis = labsize);
		}

		nl	= length(unique(rank(web))) - 1;
		lev	= as.numeric(names(table(web)));

		# set color of squares
		# and draw them
		for (i in 1:n_a) {
			for (j in 1:n_b) {
				if (!weighted) {
					if(web[n_a - i + 1, j] > 0) {
						squareColor = rgb(0,0,1);				# set color of ij.th square to blue if interaction exists
					}
					else { squareColor = rgb(1,1,1); }			# else set to white
				}
				else {
					if (rank) {
						red  = green = (1 - (which(lev == web[n_a - i + 1, j]) - 1)/nl);
						blue = 1;
					}
					else {
						red  = green = 1 - web[n_a - i + 1, j]/max(web);
						blue = 1;		# set color of ij.th square to blue with color intensity proportional to number of interaction
					}

					squareColor = rgb(red, green, blue);
				}
				rect(j - 1, i - 1, j, i, col = squareColor, border=square.border);
			}
		}

		if (plotModules) {
			foundModules = getModuleCoordinates(moduleWebObject, fromDepth, upToDepth);
			drawModules(foundModules);
		}
	}
}
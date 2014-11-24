# This function plots rectangles representing modules onto the web plotted by function plotModuleWeb(...).
# The coordinates are previously computed by function getModuleCoordinates(...).
drawModules = function(foundModules) {

	nrOfFoundModules = nrow(foundModules);

	colors = c("black", "red", "blue", "green", "yellow", "purple", "orange", "grey");

	for(i in 1:nrOfFoundModules) {
		rect(foundModules[i,2], foundModules[i,3], foundModules[i,4], foundModules[i,5], border=colors[(foundModules[i,1] %% length(colors)+1)]);	# plot modules by drawing a border around the corresponding squares
		rect(-100, foundModules[i,3], 0, foundModules[i,3], border=colors[(foundModules[i,1] %% length(colors)+1)]);					# Auxiliary line left
		rect(foundModules[i,2], -50, foundModules[i,2], 0, border=colors[(foundModules[i,1] %% length(colors)+1)]);					# Auxiliary line bottom
		rect(-100, foundModules[i,5], 0, foundModules[i,5], border=colors[(foundModules[i,1] %% length(colors)+1)]);					# Auxiliary line left 2
		rect(foundModules[i,4], -50, foundModules[i,4], 0, border=colors[(foundModules[i,1] %% length(colors)+1)]);					# Auxiliary line bottom 2
	}
}

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
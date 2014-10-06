plotPAC <- function(web, scaling=1, plot.scale=1, fill.col=rgb(.2,.2,.2,.5), arrow.col=rgb(.5,.5,.5,.5), outby=1, text=TRUE, circles=FALSE, radius=1){
	# function to draw a circular PAC-plot, as in Morris et al. 2005
	# PAC is the "Potential for Apparent Competition and is computed using the function with the same name in bipartite
	# by default, this function yields a plot for the lower trophic level
	# author: Carsten Dormann, 07 Sept 2009
	#
	# web 		a community matrix with two trophic levels
	
	toCartesian <- function (t1, rP) {
		# I stole this function from the package fisheyeR (sorry, but it was not worth including it as dependent only for three lines of code)
    	x1 = rP * cos(t1)
    	y1 = rP * sin(t1)
    	return(cbind.data.frame(x = x1, y = y1))
	}
	
	toPolar <- function (xy){
		# same source as toCartesian ...
		# nicked and vectorised
		xy <- t(as.matrix(xy))
		t1 = atan2(xy[,2], xy[,1])
   		rP = sqrt(xy[,1]^2 + xy[,2]^2)
	    return(c(t1 = t1, rP = rP))
	}

	pointsoncircle <- function(N){
		# helper function
		# computes positions of equidistant points (i.e. higher trophic level species) on a circle
		
		rhos <- seq(0, 2*pi, length=N+1)
		out <- as.matrix(toCartesian(rhos, 1)[-(N+1),2:1])
		colnames(out) <- c("x", "y")
		out
	}
	
	coords <- pointsoncircle(NROW(web))
	rs <- rowSums(web)
		
	# plot position and size of species:
	par(mar=c(0,0,0,0)+.1)
	plot(coords, cex=sqrt(rs)*0.75*scaling, xlab="", ylab="", axes=FALSE, xlim=c(-1, 1)*1.25*plot.scale, ylim=c(-1, 1)*1.25*plot.scale,asp=1)
	

	# compute PACs:
	PV <- PAC(web)
	
	# plot self-loop (i.e. diagonals) as filling:
	D <- diag(PV)
	points(coords, cex=sqrt(rs)*0.75*scaling*D, pch=16, col=fill.col)


	# draw PAC-triangles (rectangles?):
	for (i in (1:NROW(PV))[order(rs)]){
		for (j in (1:NROW(PV))[order(rs)]){
			if (i <= j) next # dAB and dBA are drawn simultaneously
			
			arrow.direction <- toPolar(coords[j,] - coords[i,])[1] #arrow from j to i
			orthog <- arrow.direction + pi/2
			
			# a cex=1 is 0.05 units diameter
			# to scale the absolute width to cex-equivalents, we need to multiply with 0.05:
			width.i <- PV[j, i]/2*0.025 *sqrt(rs[i]) *0.75*scaling # /2 because the width goes in both directions later
			width.j <- PV[i, j]/2*0.025 *sqrt(rs[j]) *0.75*scaling
			
			upper.i <- coords[i,] + toCartesian(orthog, width.i)
			lower.i <- coords[i,] - toCartesian(orthog, width.i)
			upper.j <- coords[j,] + toCartesian(orthog, width.j)
			lower.j <- coords[j,] - toCartesian(orthog, width.j)
			polygon(rbind(upper.i, lower.i, lower.j, upper.j), col=arrow.col, border=NA) #from j to i
			
		}
	}
	
	if (text) text(coords*1.25*outby, as.character(1:NROW(web)))
	if (circles) symbols(coords*1.25*outby, circles=rep(0.07*radius, NROW(web)), add=TRUE, inches=FALSE)
	
}
#-- plotwebr = a revised new implementation of plotweb, for easier making changes ----
# by Jochen Fründ, May2021
# a quick and dirty attempt for a new plotweb version "engine"; many things taken from previous plotweb

# for development only
# web <- testweb

plotwebr <- function(web, 
  space.perc = c(10,10), # space % between boxes (lower, higher)
  scale.spacing = 0.01, # factor for increasing the space.perc with increasing number of species (keep at 0 if setting a custom spacing with space.perc!!)
  # basic color stuff: better allow vectors for high, low, interactions [with these following the higher/lower sp!]
  col.boxes = c("darkgreen","grey10"), # name the elements?, maybe call it col.species?
  col.int = "grey80",  # maybe change back to original name "col.interaction"? or partial matching
  col.add_abun = c("green","green"),
  box.height = c(0.1,0.1), # low, high; in proportion of vertical size of one network 
  text.rot = c(30,30), # low, high;
  abbr.lab = c(NA, 6), # low, high; if not NA, using abbreviate with minlength = value; cf. lablength in old plotweb
  sequence = NULL,
  add = FALSE,
  x.lim = c(0,1),
  y.lim = c(0,1),
  # y.lim = c(0, 1 + 1+0.3), # for two plots
  mar = c(4,3,4,3),  # if NULL, preserve user setting of par
  add_abun.low = NULL, 
  add_abun.high = NULL,
  plot.add_abun = FALSE, # FALSE is a good default for plot2webs, but should be TRUE for single use of plotwebr
  rescale.boxwidth = TRUE,
  method = "cca",
  empty = TRUE
  )
{
  # setting margins # make this settable for user! (old plotweb overwrites!)
  if (is.null(mar)) {mar <- par()$mar} else {par(mar = mar)}  
  
  # initialize plot
  if (add==FALSE){
    yshift <- 0
    # make the yshift more flexible!! (for plot2webs)
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i") # original plotweb code, but option plot.axes removed
    # could use par(xaxs = "i", yaxs = "i") for making the plotting region follow xlim/ylim exactly
  } else {
    yshift <- 1+0.3
  }
  
  #-- preparatory calculations --
  nr <- nrow(web)
  nc <- ncol(web)
  if (is.null(add_abun.low)) add_abun.low <- rep(0,nr)
  if (is.null(names(add_abun.low))) names(add_abun.low) <- rownames(web)
  if (is.null(add_abun.high)) add_abun.high <- rep(0,nc)
  if (is.null(names(add_abun.high))) names(add_abun.high) <- colnames(web)
  
  # rearrangement of web: now outsourced!
  web <- sortweb2(web, sequence=sequence, empty=empty, sort.order=method)
  
  # also re-sort abundances!
  add_abun.low <- add_abun.low[rownames(web)]
  add_abun.high <- add_abun.high[colnames(web)]
  
  # labels etc (truncation only after re-sorting abuns)
  labels.low <- rownames(web)
  if (!is.na(abbr.lab[1])) {labels.low <- abbreviate(labels.low, minlength=abbr.lab[1])} 
  labels.high <- colnames(web)
  if (!is.na(abbr.lab[2])) {labels.high <- abbreviate(labels.high, minlength=abbr.lab[2])}

  # aesthetics preparations
  srt.low <- text.rot[1]
  srt.high <- ifelse(length(text.rot)==1, text.rot[1], text.rot[2])
  col.low <- col.boxes[1]
  col.high <- ifelse(length(col.boxes)==1, col.boxes[1], col.boxes[2]) # allow similar functionality elsewhere!
  col.add_abun.low <- col.add_abun[1]
  col.add_abun.high <- ifelse(length(col.add_abun)==1, col.add_abun[1], col.add_abun[2]) # allow similar functionality elsewhere!
  #! if color is a vector, don't forget to also resort it! (unless strictly calling by name?)
  
  # web-stuff
  websum <- sum(web)
  nr <- nrow(web) # recalculating species numbers here, as it might have changed with empty
  nc <- ncol(web)
  prop.low <- rowSums(web) / (websum + sum(add_abun.low))
  prop.high <- colSums(web) / (websum + sum(add_abun.high))
  prop.add.low <- add_abun.low / (websum + sum(add_abun.low))
  prop.add.high <- add_abun.high / (websum + sum(add_abun.high))
  # with add_abun, size of one guild might have to be rescaled for parallel shape of link edge (as in plotweb), e.g. in host-parasitoid web, but maybe not e.g. with arrow
  if (rescale.boxwidth){
    rescale.low <- min(1, (websum + sum(add_abun.low)) / (websum + sum(add_abun.high)))
    rescale.high <- min(1, (websum + sum(add_abun.high)) / (websum + sum(add_abun.low)))
  }
  
  # calculate box-spacing automatically!
    # for better plots of large webs, total space is not fixed %, but increases with No. species (of guild with more species)
  space.low <- space.perc[1] / (100*(nr-1)) * (1 + max(nr,nc)*scale.spacing)
  space.high <- space.perc[2] / (100*(nc-1)) * (1 + max(nr,nc)*scale.spacing)

  
  #-- lower boxes --
  coord.low.xl <- c(0, cumsum(prop.low[-nr])) + c(0, cumsum(rep(space.low, nr-1))) + c(0, cumsum(prop.add.low[-nr]))
  coord.low.xr <- cumsum(prop.low) + c(0, cumsum(rep(space.low, nr-1))) + c(0, cumsum(prop.add.low[-nr]))
  coord.addlow.xl <- coord.low.xr
  coord.addlow.xr <- coord.low.xr + prop.add.low
  center <- function(x){x + (1-rescale.low)/2} # a specific function for centering in case the x-spread of boxes is below 1
  # rescale all coords to maximum of 1 (or possibly lower if add_abun given)
  coord.low.xl <- center(coord.low.xl * rescale.low / coord.addlow.xr[nr])
  coord.low.xr <- center(coord.low.xr * rescale.low / coord.addlow.xr[nr])
  coord.addlow.xl <- center(coord.addlow.xl * rescale.low / coord.addlow.xr[nr])
  coord.addlow.xr <- center(coord.addlow.xr * rescale.low / coord.addlow.xr[nr])
  # draw the boxes
  rect(coord.low.xl, 0 + yshift,  coord.low.xr, box.height[1] + yshift, col=col.low)
  # optional plotting of additional abundance boxes (not just whitespace)
  if (plot.add_abun){
    rect(coord.addlow.xl, 0 + yshift,  coord.addlow.xr, box.height[1] + yshift, col=col.add_abun.low)
  }
  
  # lower labels
  #! hoffset-Zeug einbauen (original function, strwidth)
  text(labels=labels.low, x=(coord.low.xl + coord.low.xr)/2, y=-0.01 + yshift, cex=0.6, adj=1, srt=srt.low, xpd=TRUE)
  
  
  #-- higher boxes --
  coord.high.xl <- c(0, cumsum(prop.high[-nc])) + c(0, cumsum(rep(space.high, nc-1))) + c(0, cumsum(prop.add.high[-nc]))
  coord.high.xr <- cumsum(prop.high) + c(0, cumsum(rep(space.high, nc-1))) + c(0, cumsum(prop.add.high[-nc]))
  coord.addhigh.xl <- coord.high.xr
  coord.addhigh.xr <- coord.high.xr + prop.add.high
  center <- function(x){x + (1-rescale.high)/2} # a specific function for centering in case the x-spread of boxes is below 1
  # rescale all coords to maximum of 1 (or possibly higher if add_abun given)
  coord.high.xl <- center(coord.high.xl * rescale.high / coord.addhigh.xr[nc])
  coord.high.xr <- center(coord.high.xr * rescale.high / coord.addhigh.xr[nc])
  coord.addhigh.xl <- center(coord.addhigh.xl * rescale.high / coord.addhigh.xr[nc])
  coord.addhigh.xr <- center(coord.addhigh.xr * rescale.high / coord.addhigh.xr[nc])
  # draw the boxes
  rect(coord.high.xl, 1 + yshift,  coord.high.xr, 1-box.height[1] + yshift, col=col.high)
  # optional plotting of additional abundance boxes (not just whitespace)
  if (plot.add_abun){
    rect(coord.addhigh.xl, 1 + yshift,  coord.addhigh.xr, 1-box.height[1] + yshift, col=col.add_abun.high)
  }
  
  # higher labels
  #! hoffset-Zeug einbauen (original function, strwidth)
  text(labels=labels.high, x=(coord.high.xl + coord.high.xr)/2, y=1.01 + yshift, cex=0.6, adj=0, srt=srt.high, xpd=TRUE)
  
  
  #-- interactions --
  web.df <- data.frame(row=rep(1:nr, nc), col=rep(1:nc, each=nr), weight=c(web))
  # web.df <- web.df[order(-web.df$weight),] # messes up positions, moved to loop
  web.df <- web.df[web.df$weight>0,]
  # XYcoords <- as.matrix(web.df[, 1:2])  # I guess not needed
  web.df[, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")] <- NA # x-coordinates of interactions: tl=topleft, etc
  
  # low coordinates for interactions (in order of the web.df)
  # for i in lower species
  for (i in unique(web.df$row)){
    # i <- 3
    links.i <- web.df[web.df$row==i, ]
    relpos <- cumsum(links.i$weight) / sum(links.i$weight)
    coords.int.low <- (coord.low.xl[i] + relpos*(coord.low.xr[i] - coord.low.xl[i]))
    web.df[web.df$row==i, "xcoord.bl"] <- c(coord.low.xl[i], coords.int.low[-nrow(links.i)])
    web.df[web.df$row==i, "xcoord.br"] <- c(coords.int.low)
    # abline(v=coords.int.high, col="red") # for checking only
  }
  
  # high coordinates for interactions (in order of the web.df)
  # for j in higher species
  for (j in unique(web.df$col)){
    # j <- 3
    links.j <- web.df[web.df$col==j, ]
    relpos <- cumsum(links.j$weight) / sum(links.j$weight)
    coords.int.high <- (coord.high.xl[j] + relpos*(coord.high.xr[j] - coord.high.xl[j]))
    web.df[web.df$col==j, "xcoord.tl"] <- c(coord.high.xl[j], coords.int.high[-nrow(links.j)])
    web.df[web.df$col==j, "xcoord.tr"] <- c(coords.int.high)
    # abline(v=coords.int.high, col="red") # for checking only
  }
  
  # loop through interactions
  for (linki in order(-web.df$weight)){
   polygon(web.df[linki, 4:7], y=c(1-box.height[2], 1-box.height[2], box.height[1], box.height[1]) + yshift, col=col.int)
  }  
}

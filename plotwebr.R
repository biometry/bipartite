#-- plotwebr = a revised new implementation of plotweb, for easier making changes ----
# by Jochen Fründ, May2021
# a quick and dirty attempt for a new plotweb version "engine"; many things taken from previous plotweb

# for development only
# web <- testweb

# Note: earlier arguments are more likely to be changed by user

plotwebr <- function(web,
  arrow = "no", # display type of connection between upper and lower boxes, options are down, up, both, down.center, up.center, both.center and no, default is no, which is a polygonal connection between boxes.
  method = "cca", # the sorting method, passed to sortweb2
  sequence = NULL,
  empty = TRUE,
  add = FALSE,
  abun.low = NULL,
  abun.high = NULL,
  add_abun.low = NULL, 
  add_abun.high = NULL,
  col.boxes = c("darkgreen","grey10"), # can be a single value, or 2 values, or a list with two elements (each being a value or vector over species of the low/high level, resp., in order! of the web); archived: name the elements? NO, maybe call it col.species? NO
  col.int = "grey80",  # maybe change back to original name "col.interaction"? or partial matching; can be a single value or a matrix; to assign interaction color by higher species: e.g. matrix(c("red",rep("grey80",ncol(web)-1)),byrow=T,nrow=nrow(web),ncol=ncol(web))
  col.add_abun = c("green","green"),
  # need support for border colors!
  border.boxes = c("black","black"),  # currently used also for add_abun
  border.int = "black", # can be a single color value or a matrix (see col.int)
  box.height = c(0.1,0.1), # low, high; in proportion of vertical size of one network
  text.rot = c(30,30), # low, high; should be between 0 and 180, otherwise positioning (adj) may fail
  abbr.lab = c("s3", "s3"), # low, high; NA does not abbreviate; a number cuts label to this number of character (like lablength in old plotweb); a number preceded by "s" (default) tries good species abbreviation by taking this number of letter from genus and species name (if those are separated by " ", "." or "_")
  x.lim = c(0,1),
  y.lim = c(0,1),
  # y.lim = c(0, 1 + 1+0.3), # for two plots
  mar = c(4,3,4,3),  # if NULL, preserve user setting of par
  plot.boxes=c(TRUE, TRUE),  # low, high; maybe set to FALSE for combined plots (e.g. multitrophic)
  plot.labels=c(TRUE, TRUE),
  plot.add_abun = FALSE, # FALSE is a good default for plot2webs, but should be TRUE for single use of plotwebr
  rescale.boxwidth = "choose",  # either TRUE, FALSE, or "choose", which uses FALSE for given abun.low/abun.high, and TRUE otherwise
  cex.lab = c(0.6,0.6),
  space.perc = c(10,10), # space % between boxes (lower, higher)
  space.scaling = 0.01 # factor for increasing the space.perc with increasing number of species (keep at 0 if setting a custom spacing with space.perc!!)
  # basic color stuff: better allow vectors for high, low, interactions [with these following the higher/lower sp!]
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

  # prepare abundance vectors:
  # give rownames & colnames if missing
  rownames(web) <- rownames(web, do.NULL = FALSE, prefix="r")
  colnames(web) <- colnames(web, do.NULL = FALSE, prefix="c")
  # create 0 add_abun if not provided
  if (is.null(add_abun.low)) add_abun.low <- rep(0, nrow(web))
  if (is.null(add_abun.high)) add_abun.high <- rep(0, ncol(web))
  # naming abundance vectors (-> now supporting unnamed input of correct order!)
  if (is.null(names(add_abun.low))) names(add_abun.low) <- rownames(web)
  if (is.null(names(add_abun.high))) names(add_abun.high) <- colnames(web)
  if (is.null(names(abun.low))) names(add_abun.low) <- rownames(web)
  if (is.null(names(abun.high))) names(add_abun.high) <- colnames(web)
  if (is.list(col.boxes)){
    if (length(col.boxes[[1]])>1 & is.null(names(col.boxes[[1]]))){
      names(col.boxes[[1]]) <- rownames(web)
    }
    if (length(col.boxes[[2]])>1 & is.null(names(col.boxes[[2]]))){
      names(col.boxes[[2]]) <- colnames(web)
    }
  }
  if (is.list(col.add_abun)){
    if (length(col.add_abun[[1]])>1 & is.null(names(col.add_abun[[1]]))){
      names(col.add_abun[[1]]) <- names(add_abun.low)
    }
    if (length(col.add_abun[[2]])>1 & is.null(names(col.add_abun[[2]]))){
      names(col.add_abun[[2]]) <- names(add_abun.high)
    }
  }
  if (is.list(border.boxes)){
    if (length(border.boxes[[1]])>1 & is.null(names(border.boxes[[1]]))){
      names(border.boxes[[1]]) <- rownames(web)
    }
    if (length(border.boxes[[2]])>1 & is.null(names(border.boxes[[2]]))){
      names(border.boxes[[2]]) <- colnames(web)
    }
  }
  if (is.matrix(col.int) & is.null(dimnames(col.int))){
    dimnames(col.int) <- dimnames(web)
  }
  if (is.matrix(border.int) & is.null(dimnames(border.int))){
    dimnames(border.int) <- dimnames(web)
  }
  
  # rearrangement of web: now outsourced!
  web <- sortweb2(web, sequence=sequence, empty=empty, sort.order=method)
 
  # preparing and re-sorting abundances according to sorted web
    #! still to do for colors etc!
  nr <- nrow(web)
  nc <- ncol(web)
  add_abun.low <- add_abun.low[rownames(web)]
  add_abun.high <- add_abun.high[colnames(web)]
  if (!is.null(abun.low)){abun.low <- abun.low[rownames(web)]}
  if (!is.null(abun.high)){abun.high <- abun.high[colnames(web)]}
  
  # labels etc (truncation only after re-sorting abuns)
  abbr.sp <- function(x, nlett=2){
    # a function to abbreviate species names (move out of plotwebr??)
    if (substr(nlett,1,1)=="s"){
      nlett <- as.numeric(substr(nlett,2,3))
      splitter <- c(" ", ".", "_")[which.max(sapply(c(" ", ".", "_"), function(mysplit)sum(grepl(mysplit, x, fixed=TRUE))))]
      step1 <- strsplit(x, split=splitter, fixed=TRUE)
      step2 <- lapply(step1, FUN=function(x){substr(x,1,nlett)})
      return(sapply(step2, paste, collapse="."))
    } else {
      return(substr(x, 0, as.numeric(nlett)))
    }
  }
  # here is where expression-style labels could be used instead
  labels.low <- rownames(web)
  if (!is.na(abbr.lab[1])) {
    # labels.low <- abbreviate(labels.low, minlength=abbr.lab[1])
    labels.low <- abbr.sp(labels.low, nlett=abbr.lab[1])
  } 
  labels.high <- colnames(web)
  if (!is.na(abbr.lab[2])) {
    # labels.high <- abbreviate(labels.high, minlength=abbr.lab[2])
    labels.high <- abbr.sp(labels.high, nlett=abbr.lab[2])
  }

  # aesthetics preparations (incl. resorting colors)
  srt.low <- text.rot[1]
  srt.high <- ifelse(length(text.rot)==1, text.rot[1], text.rot[2])
  col.low <- unlist(col.boxes[1])
  col.high <- unlist(ifelse(length(col.boxes)==1, col.boxes[1], col.boxes[2])) # allow similar functionality elsewhere!
  if (length(col.low)>1) col.low <- col.low[rownames(web)]
  if (length(col.high)>1) col.high <- col.high[colnames(web)]
  border.low <- border.boxes[1]
  border.high <- ifelse(length(border.boxes)==1, border.boxes[1], border.boxes[2]) 
  if (length(border.low)>1) border.low <- border.low[rownames(web)]
  if (length(border.high)>1) border.high <- border.high[colnames(web)]
  col.add_abun.low <- unlist(col.add_abun[1])
  col.add_abun.high <- unlist(ifelse(length(col.add_abun)==1, col.add_abun[1], col.add_abun[2]))
  if (length(col.add_abun.low)>1) col.add_abun.low <- col.add_abun.low[rownames(web)]
  if (length(col.add_abun.high)>1) col.add_abun.high <- col.add_abun.high[colnames(web)]
  if (is.matrix(col.int)) {
    col.int <- col.int[rownames(web), colnames(web)]
  }
  if (is.matrix(border.int)) {
    border.int <- border.int[rownames(web), colnames(web)]
  }
  
  
  # preparations for calculating coordinates
  if (rescale.boxwidth=="choose"){
        rescale.boxwidth <- is.null(abun.low) & is.null(abun.high)
  }
  if (is.null(abun.low)){abun.low <- rowSums(web)}
  prop.low <- abun.low / (sum(abun.low) + sum(add_abun.low))
  prop.add.low <- add_abun.low / (sum(abun.low) + sum(add_abun.low))
  if (is.null(abun.high)){abun.high <- colSums(web)}
  prop.high <- abun.high / (sum(abun.high) + sum(add_abun.high))
  prop.add.high <- add_abun.high / (sum(abun.high) + sum(add_abun.high))
  # with add_abun, size of one guild might have to be rescaled for parallel shape of link edge (as in plotweb), e.g. in host-parasitoid web, but maybe not e.g. with arrow
  if (rescale.boxwidth){
    rescale.low <- min(1, (sum(abun.low) + sum(add_abun.low)) / (sum(abun.high) + sum(add_abun.high)))
    rescale.high <- min(1, (sum(abun.high) + sum(add_abun.high)) / (sum(abun.low) + sum(add_abun.low)))
  } else {
    rescale.low <- 1
    rescale.high <- 1
  }
  
  # calculate box-spacing automatically!
    # for better plots of large webs, total space is not fixed %, but increases with No. species (of guild with more species)
  space.low <- space.perc[1] / (100*(nr-1)) * (1 + max(nr,nc)*space.scaling)
  space.high <- space.perc[2] / (100*(nc-1)) * (1 + max(nr,nc)*space.scaling)

  
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
  if (plot.boxes[1]){
    rect(coord.low.xl, 0 + yshift,  coord.low.xr, box.height[1] + yshift, col=col.low, border=border.low)
  }
  # optional plotting of additional abundance boxes (not just whitespace)
  if (plot.add_abun){
    rect(coord.addlow.xl, 0 + yshift,  coord.addlow.xr, box.height[1] + yshift, col=col.add_abun.low, border=border.low)
  }
  
  # lower labels
  if (plot.add_abun) {
    # maybe I should always use this choice, for good labels in plot2webs?
    coord.lab.low <- (coord.low.xl + coord.addlow.xr)/2
  } else {
     coord.lab.low <- (coord.low.xl + coord.low.xr)/2   
  }
  labshift <- numeric(nr)
  labelheight <- strheight(labels.low, cex=cex.lab[2])
  # shifting labels up to avoid overwriting (only with hoirzontal labels)
  if (nr > 1 & (srt.low==0 | srt.low==180)){
    labelwidth <- strwidth(labels.low, cex=cex.lab[2])
    label.xl <- coord.lab.low - labelwidth/2
    label.xr <- coord.lab.low + labelwidth/2
    for (j in 2:nr){
      labshift[j] <- 0
      while (suppressWarnings(max(label.xr[labshift==labshift[j] & (1:nr)<j]) >= label.xl[j])){
        labshift[j] <- labshift[j] + 1
      }
    }
  }  
  adj.low <- c(ifelse(srt.low %in% c(0,180), 0.5, 1), ifelse(srt.low>90,-0.5,1))
  # plotting the labels:
  if (plot.labels[1]){
    text(labels=labels.low, x=coord.lab.low, y=-0.01 + yshift - labshift*labelheight, cex=cex.lab[1], adj=adj.low, srt=srt.low, xpd=TRUE)
  }
  
  
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
  if (plot.boxes[2]){
    rect(coord.high.xl, 1 + yshift,  coord.high.xr, 1-box.height[1] + yshift, col=col.high, border=border.high)
  }
  # optional plotting of additional abundance boxes (not just whitespace)
  if (plot.add_abun){
    rect(coord.addhigh.xl, 1 + yshift,  coord.addhigh.xr, 1-box.height[1] + yshift, col=col.add_abun.high, border=border.high)
  }
  
  # higher labels
  if (plot.add_abun) {
    # maybe I should always use this choice, for good labels in plot2webs?
    coord.lab.high <- (coord.high.xl + coord.addhigh.xr)/2
  } else {
     coord.lab.high <- (coord.high.xl + coord.high.xr)/2   
  }
  labshift <- numeric(nc)
  labelheight <- strheight(labels.high, cex=cex.lab[2])
  # shifting labels up to avoid overwriting (only with hoirzontal labels)
  if (nc > 1 & (srt.high==0 | srt.high==180)){
    labelwidth <- strwidth(labels.high, cex=cex.lab[2])
    label.xl <- coord.lab.high - labelwidth/2
    label.xr <- coord.lab.high + labelwidth/2
    for (j in 2:nc){
      labshift[j] <- 0
      while (suppressWarnings(max(label.xr[labshift==labshift[j] & (1:nc)<j]) >= label.xl[j])){
        labshift[j] <- labshift[j] + 1
      }
    }
  }  
  adj.high <- c(ifelse(srt.high %in% c(0,180), 0.5, 0), 0.5*sin(srt.high*pi/180))
  # plotting the labels:
  if (plot.labels[2]){
    text(labels=labels.high, x=coord.lab.high, y=1.01 + yshift + labshift*labelheight, cex=cex.lab[2], adj=adj.high, srt=srt.high, xpd=TRUE)
  }
  
  
  #-- interactions --
  web.df <- data.frame(row=rep(1:nr, nc), col=rep(1:nc, each=nr), weight=c(web), col.int = c(col.int), border.int = c(border.int))
  web.df <- web.df[web.df$weight>0,]
  web.df[, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")] <- NA # x-coordinates of interactions: tl=topleft, etc
  
  # low coordinates for interactions (in order of the web.df)
  for (i in unique(web.df$row)){   # for i in lower species
    # i <- 3
    links.i <- web.df[web.df$row==i, ]
    relpos <- cumsum(links.i$weight) / sum(links.i$weight)
    coords.int.low <- (coord.low.xl[i] + relpos*(coord.low.xr[i] - coord.low.xl[i]))
    web.df[web.df$row==i, "xcoord.bl"] <- c(coord.low.xl[i], coords.int.low[-nrow(links.i)])
    web.df[web.df$row==i, "xcoord.br"] <- c(coords.int.low)
    if (arrow %in% c("down.center", "both.center")){
      web.df[web.df$row==i, "xcoord.bl"] <- web.df[web.df$row==i, "xcoord.br"] <- mean(c(coord.low.xl[i], coord.low.xr[i]))
    }    
  }
  if (arrow %in% c("down","both")){
    web.df[, "xcoord.bl"] <- web.df[, "xcoord.br"] <- rowMeans(web.df[, c("xcoord.bl", "xcoord.br")])
  }
  
  # high coordinates for interactions (in order of the web.df)
  for (j in unique(web.df$col)){   # for j in higher species
    # j <- 3
    links.j <- web.df[web.df$col==j, ]
    relpos <- cumsum(links.j$weight) / sum(links.j$weight)
    coords.int.high <- (coord.high.xl[j] + relpos*(coord.high.xr[j] - coord.high.xl[j]))
    web.df[web.df$col==j, "xcoord.tl"] <- c(coord.high.xl[j], coords.int.high[-nrow(links.j)])
    web.df[web.df$col==j, "xcoord.tr"] <- c(coords.int.high)
    if (arrow %in% c("up.center", "both.center")){
      web.df[web.df$col==j, "xcoord.tl"] <- web.df[web.df$col==j, "xcoord.tr"] <- mean(c(coord.high.xl[j], coord.high.xr[j]))
    }  
  }
  if (arrow %in% c("up","both")){
    web.df[, "xcoord.tl"] <- web.df[, "xcoord.tr"] <- rowMeans(web.df[, c("xcoord.tl", "xcoord.tr")])
  }
  
  # loop through interactions
  for (linki in order(-web.df$weight)){
   polygon(web.df[linki, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")], y=c(1-box.height[2], 1-box.height[2], box.height[1], box.height[1]) + yshift, col=web.df$col.int[linki], border=web.df$border.int[linki])
  }  
}

#!! currently not working!

#-- plot2webs: 2-web alignment function (replaces plotweb2 with very different approach) -------
# by Jochen Fründ, May2021
  #- input: webarray or 2 webs

# for development only
# web1 <- testweb1
# web2 <- testweb2


plot2webs <- function(web1, web2,
                      flipweb2=TRUE,
                      sorting="pooledweb",
                      # optional arguments for sortweb2
                      sequence=NULL,
                      col.boxes=c("grey10","darkgreen"),
                      abbr.lab=c(NA,6),
                      ydist.add = 0.2,
                      plot.abun_other = FALSE # should a species' excess abundance in the other web of the pair be plotted as additional abundance or left empty? (used as setting for plot.add_abun in the call to plotwebr; note that externally provided "additional abundances" are not supported in plot2webs)
                      ){
  # sort the webs together!
    # allow webs2array to be called in here; currently assuming the webarray is already given 
  if (length(dim(web1))==3){
    mywebarray <- web1
    # should empty() the array! but not currently implemented?!?
    web1 <- mywebarray[,,1]
    web2 <- mywebarray[,,2]
  } else {
    mywebarray <- webs2array(list(web1, web2))
    web1 <- mywebarray[,,1] # now with matched names
    web2 <- mywebarray[,,2] # now with matched names
    # Note: if original webs are emptied, non-interacting species are missing from web1 / web2 (different from case with web1 being a webarray)
    #--> let's see if this causes problems, or gives nice options
    # can always reassign webs from array, as in option before the else
  }
  
  # three cases of sorting: given sequence OR cca of pooledweb OR by names of web1, with species unique to web2 added afterwards
  if (sorting=="pooledweb"){
    if (!is.null(sequence)){warning("given sequence is overridden by sorting='pooledweb'")}
    #? how to deal with empty here?
    sorted.aggweb <- sortweb2(apply(mywebarray, 1:2, sum))
    sequence <- dimnames(sorted.aggweb)
    names(sequence) <- c("seq.low","seq.high")
  } else { # sorts by given sequence
    if (is.null(sequence)){ # not sure if this works (below)
      sequence <- list(seq.low=unique(c(rownames(web1),rownames(web2))),
                    seq.high=unique(c(colnames(web1),colnames(web2))))
    } 
  }
  
  # calculate the "additional abundances" as the excess abundances in one vs. other web
  #! currently only works without given external abundances!! (should support independent abuns, but not given additional abuns)
  add_abun1.low <- rowSums(web2) - rowSums(web1)
    add_abun1.low[add_abun1.low < 0] <- 0
  add_abun1.high <- colSums(web2) - colSums(web1)
    add_abun1.high[add_abun1.high < 0] <- 0
  add_abun2.low <- rowSums(web1) - rowSums(web2)
    add_abun2.low[add_abun2.low < 0] <- 0    
  add_abun2.high <- colSums(web1) - colSums(web2)
    add_abun2.high[add_abun2.high < 0] <- 0
  
  # plot the webs
  plotwebr(web1, y.lim = c(0, 1 + 1+ydist.add), empty=FALSE, method="sequence", sequence=sequence, col.boxes = col.boxes, abbr.lab = abbr.lab, add_abun.low = add_abun1.low, add_abun.high = add_abun1.high, plot.add_abun = plot.abun_other)
  if (flipweb2){ # often nice for web alignment [multilayer, centering the shared species], but not for multitrophic webs!
    web2 <- t(web2)
    # flipping the sequence: need to switch the vectors, but keep the names
    seq2 <- sequence
    seq2[2:1] <- sequence
    plotwebr(web2, add=TRUE, empty=FALSE, method="sequence", sequence=seq2, col.boxes = col.boxes[2:1], abbr.lab = abbr.lab[2:1], add_abun.low = add_abun2.high, add_abun.high = add_abun2.low, plot.add_abun = plot.abun_other, ydist.add=ydist.add, plot.labels = c(FALSE, TRUE))
    # beware of changes in col.boxes!
  }
  
  #! implement option without flipping web2 (for multitrophic webs, and also the starting point for series of more than 2 webs)
  # plotwebr(web2, add=TRUE, col.boxes = col.boxes, abbr.lab = c(6,NA), empty=FALSE, method="normal") # use abbr.lab = c(0,NA) to omit the lower labels in upper plot
}

# example:
# plot2webs(Safariland, t(vazarr), sorting="given")
# plot2webs(webs2array(Safariland, vazarr))

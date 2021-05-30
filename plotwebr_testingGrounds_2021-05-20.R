#-- the main file for development of the plotweb revisions ----
# by Jochen Fründ, May2021

library(bipartite)

# three new functions (helper, main, meta)

#*** To Do sortweb2 ***
#! change name to sortweb (also checking dependencies)
  # check if problems occur because new default sorting of sortweb is not longer "dec"!
  # then use sortweb in plotwebr and plot2webs
#! make helpfile
  #- note that the sequence now must have names seq.low and seq.high, for consistency with plotweb! (was seq.lower and seq.higher)
  #- partial matching now allowed, e.g. "dec" or "decreasing"
#! check through code, comment, and extract to do from comments
source("sortweb2.R")

#*** To Do plotwebr ***
#! check through code, comment, and extract to do from comments
#! color customization per species (as in plotweb), and smarter for interactions (e.g. based on the species)
#! cleanup arguments / sort them better, etc
#! add all customization options from plotweb
#! make helpfile
  #- note that !is.null(sequence) now changes method (sort.order in sortweb) to "seq"
#! check if it works with many different webs, including those with 1sp in a guild, 0cols, etc
source("plotwebr.R")  

#*** To Do plot2webs ***
#! check through code, comment, and extract to do from comments
#! fix it to work provisionally (with array or 2 webs)
#- with or without labels, space in between or not, etc.
#- add ... functionality for (some) plotweb arguments
#- (vertical? ggf einfach über gutes text.rot)
#! should I rather flipweb1 instead of web2? (i.e. is it more common to center higher or lower??)
#! make helpfile
source("plot2webs.R")

#*** To Do all functions ***
#* rename argument "sequence"? conflicts with base function; maybe use specseq?
#* maybe visweb could be taken next, calling sortweb and making arguments more consistent with plotweb?


# before installing latest bipartite, call the fixed webs2array function:
source("bipartite/R/webs2array.R")


#-- developmental testing grounds ----
testweb1 <- Safariland[1:4,1:5]
testweb2 <- Safariland[c(1:3,5), c(1:3, 7:9)]

testweb1 <- memmott1999#[1:2,]

testweb1 <- Safariland
testweb2 <- vazarr

# plotweb(testweb1)
plotwebr(testweb1, col.boxes=c("darkgreen", "purple"))
plotwebr(testweb2, empty=TRUE)

plotwebr(testweb1, col.boxes=list(c("darkgreen", "purple"),"orange"))

# a flashy version
plotwebr(testweb1, col.boxes=list(c("darkgreen", "purple", "turquoise", "goldenrod"),"orange"), border.boxes=NA, border.int="pink", method="inc")

plotwebr(testweb1, col.boxes=list(c("darkgreen", "purple", "turquoise", "goldenrod"),"orange"), border.boxes=NA, border.int="pink", method="inc")


webcolors <- matrix(c("grey80", "red", "grey80", "grey80", "grey80"),byrow=T,nrow=nrow(testweb1),ncol=ncol(testweb1))
plotwebr(testweb1, col.boxes=list("darkgreen", c("purple", "red",  "purple", "purple",  "purple")), col.int = webcolors, border.int = webcolors )


# testweb1 <- testweb1 * (sum(testweb2) / sum(testweb1))

plot2webs(webs2array(testweb1, testweb2))

plot2webs(testweb1, testweb2) # not yet working


source("plotwebr.R")

#-- playground -------

#locator(1)

testweb <- Safariland[1:4,1:5]
testweb <- Safariland
testweb <- memmott1999
plotwebr(testweb)
plotwebr(testweb, y.lim = c(0, 1 + 1+0.3))
testweb <- t(vazarr)
plotwebr(testweb, add=T)

plotwebr(testweb, text.rot=45)
plotwebr(testweb, text.rot=0)
plotwebr(testweb, text.rot=90)
plotwebr(testweb, text.rot=180)
plotwebr(testweb, text.rot=135)
plotwebr(testweb, text.rot=90, cex.lab=c(0.6,0.4))
plotwebr(testweb, text.rot=c(315,45))

plotwebr(testweb, text.rot=0, abbr.lab = c(3,3), space.scaling = 0.04)
plotweb(testweb, text.rot=0, high.lablength = 3)

plotwebr(testweb, text.rot=0, abbr.lab = c("s2","s2"), space.scaling = 0.04)

plotwebr(testweb, text.rot=0, abbr.lab = c("s2","s2"), space.scaling = 0.04, cex.lab=c(0.6,0.5))

# for Safariland
plotwebr(testweb, add_abun.low = c(500,1:8), add_abun.high=seq(10,270, by=10), plot.add_abun = T)

# for subset of Safari
plotwebr(testweb, add_abun.low = c(500,50,100,200), add_abun.high=seq(10,50, by=10), plot.add_abun = T, col.add_abun = list(c("white","lightblue","white","yellow"), rep("red",5)))

plotweb(testweb)

testweb <- Safariland[1:4,1:5]
abun.high <- 5:1
names(abun.high) <- colnames(testweb)
abun.low <- rep(100,4)
names(abun.low) <- rownames(testweb)
plotweb(testweb, empty = FALSE)
plotweb(testweb, empty = FALSE, high.abun=abun.high, low.abun=abun.low)
plotweb(testweb, empty = FALSE, high.abun=abun.high, low.abun=abun.low, abuns.type="independent")

# hidden additional abuns
plotwebr(testweb, empty = FALSE, add_abun.high=abun.high, add_abun.low=abun.low)
# additional abuns
plotwebr(testweb, empty = FALSE, add_abun.high=abun.high, add_abun.low=abun.low, plot.add_abun = TRUE)
# independent abuns
plotwebr(testweb, empty = FALSE, abun.high=abun.high, abun.low=abun.low, plot.add_abun = TRUE)
plotwebr(testweb, empty = FALSE, abun.high=abun.high, abun.low=abun.low, plot.add_abun = TRUE, arrow="down")
plotwebr(testweb, empty = FALSE, abun.high=abun.high, abun.low=abun.low, plot.add_abun = TRUE, arrow="both")
plotwebr(testweb, empty = FALSE, abun.high=abun.high, abun.low=abun.low, plot.add_abun = TRUE, arrow="down.center")

plotwebr(testweb, empty = FALSE, abun.high=abun.high, abun.low=abun.low, plot.add_abun = TRUE, arrow="down.center", add_abun.high=rep(10,5))

plotwebr(testweb, plot.boxes=c(T,F), plot.labels=c(T,F), border.boxes=NA, border.int=NA)

#-- Example code from old plotweb help (edited for plotwebr dvlp) ---------
data(Safariland)
plotweb(Safariland)

# shorter labels
plotweb(Safariland, high.lablength=3, low.lablength=0, arrow="down")

# centered triangles for displaying interacions
plotwebr(Safariland, text.rot=90, arrow="down.center", col.int="wheat2",
	y.lim=c(-1,2.5))

#orginal sequence, up arrows and different box width
plotweb(Safariland, method="normal", arrow="up", y.width.low=0.3, low.lablength=4)
# interactions as lines
plotweb(Safariland, arrow="both", y.width.low=0.05, text.rot=90, col.high="blue", 
	col.low="green")

# add an abundance vector for lower trophic species 
low.abun = round(runif(dim(Safariland)[1],1,400)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction="purple", 
	y.width.low=0.05, y.width.high=0.05)
# plotwebr: 
plotwebr(Safariland, text.rot=90, col.int="purple")

high.abun = round(runif(dim(Safariland)[2],1,400)) #create
names(high.abun) <- colnames(Safariland)
plotwebr(Safariland, text.rot=90, add_abun.high=high.abun, col.int="purple",plot.add_abun = TRUE)
plotwebr(Safariland, text.rot=90, add_abun.high=high.abun, col.int="purple",plot.add_abun = FALSE)

plotwebr(Safariland, text.rot=90, add_abun.low=low.abun, col.int="purple",plot.add_abun = TRUE)
plotwebr(Safariland, text.rot=90, add_abun.low=low.abun, col.int="purple",plot.add_abun = TRUE, rescale.boxwidth = FALSE)
plotwebr(Safariland, text.rot=90, add_abun.low=low.abun, col.int="purple",plot.add_abun = TRUE, rescale.boxwidth = FALSE, space.perc = c(10,30))


plotweb(Safariland, text.rot=90, low.abun=low.abun, col.interaction ="red", 
	bor.col.interaction="red", arrow="down")


# now vectors for all colours can be given, to mark certain species or 
# interactions. Colour vectors are recycled if not of appropriate length
plotweb(Safariland,col.high=c("orange","green"))
plotweb(Safariland,col.low=c("orange","green"),col.high=c("white","grey","purple"),
	text.high.col=c("blue","red"), col.interaction=c("red",rep("green",26),rep("brown",242)),
	bor.col.interaction=c(rep("green",26),rep("brown",242)),method="normal", 
	text.rot=90, low.lablength=10, high.lablength=5)


#example one (tritrophic)
plotweb(Safariland,y.width.low=0.1, y.width.high=0.05,method="normal", 
	y.lim=c(0,3), arrow="up", adj.high=c(0.5,1.5), col.high="orange",
	high.lablength=3,high.lab.dis=0)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal",
	add=TRUE,low.y=1.5,high.y=2.5, col.low="green", text.low.col="red", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1),low.plot=FALSE)

#example two (4 trophic with abundance)
low.abun = round(runif(dim(Safariland)[1],1,40)) #create
names(low.abun) <- rownames(Safariland)
plotweb(Safariland, text.rot=90, high.abun=low.abun, col.interaction="purple", 
	y.lim=c(0,4.5), high.lablength=0, arrow="up", method="normal", 
	y.width.high=0.05)

plotweb(t(Safariland), y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=1.7,high.y=2.7, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4, 
	high.lablength=0)

plotweb(Safariland,y.width.low=0.05, y.width.high=0.1, method="normal", 
	add=TRUE, low.y=2.95, high.y=3.95, col.low="green", text.low.col="black", 
	low.lab.dis=0, arrow="down", adj.low=c(0.5,1.1), low.lablength=4)

# now some examples with the abuns.type-option:
plotweb(Safariland, abuns.type='independent',arrow="down.center")
plotweb(Safariland, abuns.type='additional',arrow="down.center")



#-- the main file for development of the plotweb revisions ----
# by Jochen Fründ, May2021

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
#! independent abundances
  #- abundance options ("independent" AND "additional", with "additional" possible to suppress [for alignment], and position left or right)
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
plotwebr(testweb1)
plotwebr(testweb2, empty=TRUE)

# testweb1 <- testweb1 * (sum(testweb2) / sum(testweb1))

plot2webs(webs2array(testweb1, testweb2))

plot2webs(testweb1, testweb2) # not yet working




#-- playground -------

#locator(1)

testweb <- Safariland[1:4,1:5]
testweb <- Safariland
testweb <- memmott1999
plotwebr(testweb)
plotwebr(testweb, y.lim = c(0, 1 + 1+0.3))
testweb <- t(vazarr)
plotwebr(testweb, add=T)

plotwebr(testweb, add_abun.low = c(500,1:8), add_abun.high=seq(100,2700, by=100))

plotweb(testweb)
abun.high <- 5:1
names(abun.high) <- colnames(testweb)
abun.low <- rep(1,4)
names(abun.low) <- rownames(testweb)
plotweb(testweb, empty = FALSE)
plotweb(testweb, empty = FALSE, high.abun=abun.high, low.abun=abun.low)
plotweb(testweb, empty = FALSE, high.abun=abun.high, low.abun=abun.low, abuns.type="independent")



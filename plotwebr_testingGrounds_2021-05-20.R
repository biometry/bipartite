#-- the main file for development of the plotweb revisions ----
# by Jochen Fründ, May2021

# three new functions (helper, main, meta)

#*** To Do sortweb2 ***
#! check through code, comment, and extract to do from comments
#! merge with sortweb !  (also allow sort.order inc, dec)
#! give warnings when method is changed to "normal" or cca overridden by sequence
  # simplify the code (issue1, issue2, issue3 together)
#! likely simplify cca-method, always using the compartment-loop even with 1 compartment
#! advanced sorting functionality
#! make helpfile
source("sortweb2.R")

#*** To Do plotwebr ***
#! check through code, comment, and extract to do from comments
#! color customization per species (as in plotweb), and smarter for interactions (e.g. based on the species)
#! independent abundances
  #- abundance options ("independent" AND "additional", with "additional" possible to suppress [for alignment], and position left or right)
#! cleanup arguments / sort them better, etc
#! add all customization options from plotweb
#! make helpfile
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

# before installing lates bipartite:
source("bipartite/R/webs2array.R")


#-- developmental testing grounds ----
testweb1 <- Safariland[1:4,1:5]
testweb2 <- Safariland[c(1:3,5), c(1:3, 7:9)]

testweb1 <- Safariland
testweb2 <- vazarr

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



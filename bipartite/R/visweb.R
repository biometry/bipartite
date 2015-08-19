visweb <- function (web, type = "nested", prednames = TRUE, preynames = TRUE, 
    labsize = 1, plotsize = 12, square = "interaction", text = "no", 
    frame = NULL, textsize = 1, textcol = "red", pred.lablength = NULL, 
    prey.lablength = NULL, clear = TRUE, xlabel = "", ylabel = "", 
    boxes = TRUE, circles = FALSE, circle.col = "black", circle.min = 0.2, 
    circle.max = 2, outerbox.border = "white", outerbox.col = "white", 
    box.border = "black", box.col = "black", def.col="blue", max.digits=4, NA.col="red") 
{
	backweb <- as.matrix(web)
    if (!is.matrix(web)) {
        
        web <- as.matrix(web)
        warning("Object converted to matrix.")
    }
    if (length(as.numeric(names(table(web))))!= length(def.col) & substr(square,1,1)=="d") {
        warning("Defined colors not of equal length to number of levels of interactions")
    }
    if (type != "diagonal" && is.null(frame)) 
        frame <- FALSE
    if (clear == TRUE) 
        web <- empty(web)
    if (type == "diagonal") {
        #web <- empty(web)
        ca <- cca(web)
        web <- web[order(summary(ca)$sites[, 1], decreasing = TRUE), 
            order(summary(ca)$species[, 1], decreasing = TRUE)]
    }
    if (type == "nested") {
        #web <- empty(web)
        web <- web[order(rowSums(web), decreasing = TRUE), order(colSums(web), 
            decreasing = TRUE)]
    }
    n.pred <- ncol(web)
    n.prey <- nrow(web)
    plotsize = plotsize/2.54
    mcol = max(web)
    if (n.pred > n.prey) {
        wx <- plotsize
        wy <- (plotsize)/n.pred * n.prey
    } else {
        wy <- plotsize
        wx <- (plotsize)/n.prey * n.pred
    }
    m.predsize = max(strwidth(colnames(web), units = "inches"))
    m.preysize = max(strwidth(rownames(web), units = "inches"))
    cellsize = wx/n.pred
    if (substr(text, 1, 1) == "i") {
        s <- as.character(round(max(web),max.digits))
    } else {s = "A"}
    lettersize = strwidth(s, units = "inches")
    clratio = cellsize/lettersize
    #mm.pred <- max(m.predsize, m.preysize)
  
    #pin = c(wx, wy),
    par(mai = c(m.predsize*0.4 * labsize * clratio+0.2, m.preysize*0.4 * labsize * clratio+0.2, 0, 0))
      
    plot(1, type = "n", axes = FALSE, xlim = c(0, n.pred), ylim = c(0, 
        n.prey), asp = 1, xlab = xlabel, ylab = ylabel)
    rect(0, 0, n.pred, n.prey, col = outerbox.col, border = outerbox.border)
    pnl <- 0
    if (prednames && !is.null(colnames(web))) {
        for (ii in 1:n.pred) {
            s <- colnames(web)[ii]
            if (!is.null(s)) {
                pop <- regexpr("_", s)
                if (pop > 0) {
                  cn <- paste(c(substr(s, 1, min(pop - 1, 9)), 
                    "\n", substr(s, pop + 1, pop + 9)), sep = "", 
                    collapse = "")
                } else {
                  if (!is.null(pred.lablength)) {
                    cn <- substr(s, 1, pred.lablength)
                  }
                  else {
                    cn <- s
                  }
                }
                pnl[ii] <- cn
            }
        }
        axis(1, (1:n.pred) - 0.5, labels = pnl, tick = FALSE, 
            mgp = c(0, 0, 0), las = 2, cex.axis = 0.4 * labsize * 
                clratio)
    }
    ynl <- 0
    if (preynames && !is.null(rownames(web))) {
        for (i in 1:n.prey) {
            s <- rownames(web)[n.prey - i + 1]
            if (!is.null(s)) {
                pop <- regexpr("_", s)
                if (pop > 0) {
                  cr <- paste(c(substr(s, 1, min(pop - 1, 9)), 
                    "\n", substr(s, pop + 1, pop + 9)), sep = "", 
                    collapse = "")
                } else {
                  if (!is.null(prey.lablength)) {
                    cr <- substr(s, 1, prey.lablength)
                  } else {
                    cr <- s
                  }
                }
                ynl[i] <- cr
            }
        }
        axis(2, (1:n.prey) - 0.5, labels = ynl, tick = FALSE, 
            mgp = c(0, 0, 0), las = 2, cex.axis = 0.4 * labsize * 
                clratio)
    }
    cross = function(web, start, comp) {
        n.r = nrow(web)
        n.c = ncol(web)
        r = start[1]
        c = start[2]
        web[r, c] = -comp
        for (i in 1:n.r) {
            if (web[i, c] > 0) web <- cross(web, c(i, c), comp)
        }
        for (i in 1:n.c) {
            if (web[r, i] > 0) web <- cross(web, c(r, i), comp)
        }
        return(web)
    }
    comp = 1
    w <- web
    while (max(w) > 0) {
        start = which(w == max(w), arr.ind = TRUE)[1, ]
        w <- cross(w, start, comp)
        comp = comp + 1
    }
    w <- abs(w)
    nl <- max(length(unique(rank(web))) - 1,1)
    lev <- as.numeric(names(table(web)))
    max.web <- max(web)
    for (i in 1:n.prey) {
        for (ii in 1:n.pred) {
            if (substr(square, 1, 1) == "c") {
                c = 1 - (w[n.prey - i + 1, ii]) * (1/(max(w)))
            } else {
              if (substr(square, 1, 1) == "i") {
                c = (1 - (which( round(lev,6) == round(web[n.prey - i + 1, ii],6)) - 1)/nl) 
              } else {
                if (substr(square, 1, 1) == "b") {
                  c = floor((1 - web[n.prey - i + 1, ii]/mcol))
                } else { 
                  c = 1
                }
              }
            }
            if (substr(square, 1, 1) == "b" & c == 0) {
              c = box.col
            } else {
              c = gray(c)
            }
            
            if (substr(square, 1, 1) == "d") c=def.col[ which(round(lev,6) == round(web[n.prey - i + 1, ii],6))]
            
            if (is.na(backweb[n.prey-i+1,ii])) c <- NA.col
            if (circles == TRUE) 
                c = outerbox.col
            if (boxes == TRUE) 
                rect(ii - 1, i - 1, ii, i, col = c, border = box.border)
            if (circles == TRUE && web[n.prey - i + 1, ii] >  0) {
                points(ii - 0.5, i - 0.5, pch = 21, col = circle.col, 
                  bg = circle.col, cex = circle.min + 
                  (circle.max * web[n.prey - i + 1, ii]/max.web))
            }
            tc <- ""
            if (substr(text, 1, 1) == "i") {
                if (web[n.prey - i + 1, ii] > 0) 
                  tc <- round(web[n.prey - i + 1, ii],max.digits)
            } else {
              if (substr(text, 1, 1) == "c") {
                if (w[n.prey - i + 1, ii] > 0) {
                  tc <- LETTERS[w[n.prey - i + 1, ii]]
                }
              }
            }
            text(ii - 0.5, i - 0.5, tc, col = textcol, cex = 1 * 
                textsize * clratio * 0.5, adj = c(0.5, 0.5))
        }
    }
    if (is.null(frame)) 
        frame = TRUE
    if (frame) {
        for (i in 1:max(abs(w))) {
            squares <- which(w == i, arr.ind = TRUE)
            minsqx <- min(squares[, 2]) - 1
            maxsqx <- max(squares[, 2])
            minsqy <- min(squares[, 1]) - 1
            maxsqy <- max(squares[, 1])
            rect(minsqx, n.prey - maxsqy, maxsqx, n.prey - minsqy, 
                lwd = 3)
        }
    }
}


#x<-matrix(c(1,0,0,1,
#0,1,0,0,
#1,1,1,0,
#0,0,1,1,
#1,0,0,0,
#1,1,1,1),nrow=6,ncol=4)
# 
#colnames(x)<-c("A","B","C","D")
#rownames(x)<-c("Aaron","Betty","Charles","Dina","Elias","Fiona")
#visweb(x,type="none",circles=TRUE, clear=F)
#
#

#m <- matrix(3,3,4)


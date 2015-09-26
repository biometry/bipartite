`plotweb2` <-
function(web,  web2, method = "cca", empty = FALSE, labsize = 1, ybig = 1,
    y_width = 0.1, spacing = 0.05, arrow="no", col.interaction="grey80",
    col.pred = "grey10", col.prey="grey10", lab.space=1,
    lablength = NULL, sequence=NULL,low.abun=NULL,low.abun.col="green",
    high.abun=NULL, high.abun.col="red",
    method2 = "cca", empty2 = TRUE, spacing2 = 0.05, arrow2="no",
    col.interaction2="grey80", col.pred2 = "grey30", col.prey2="grey20",
    lablength2 = NULL, sequence.pred2=NULL,low.abun2=NULL,low.abun.col2="green",
    high.abun2=NULL, high.abun.col2="red")



{
  if (empty) web <- empty(web) else method <- "normal"
  web<-as.matrix(web) # to convert data.frames into matrix: needed for cumsum

  meths <- c("normal", "cca")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) stop("Choose plot-method: normal/cca.\n")
  if (meths.match==2)
  {
    ca <- cca(web)
    web <- web[order(summary(ca)$sites[,1], decreasing=TRUE), order(summary(ca)$species[,1], decreasing=TRUE)]

  }

        if (!is.null(sequence)) {
            cs <- sequence$seq.pred %in% colnames(web)
            rs <- sequence$seq.prey %in% rownames(web)
            web <- web[sequence$seq.prey[rs], sequence$seq.pred[cs]]
        }

   websum <- sum(web)
   difff <- diffh <-0 # if no abundances are set leave plotsize as before

        ###rearrange if lowfreq is set  # lowfreq is a named vector!!!
        if (!is.null(low.abun)) {

        lowfreq = rowSums(web)

        dummy <- lowfreq

        for (i in 1:length(low.abun) )
        {
        ind <- which(names(low.abun)[i] == names(dummy))
        
        lowfreq[ind] <- lowfreq[ind]+low.abun[i]

        }
        #websum <- sum(lowfreq)
        difff = (lowfreq-rowSums(web))/websum
        }

        ###rearrange if highfreq is set  # lowfreq is a named vector!!!
        if (!is.null(high.abun)) {

        highfreq = colSums(web)

        dummy <- highfreq

        for (i in 1:length(high.abun) )
        {
        ind <- which(names(high.abun)[i] == names(dummy))

        highfreq[ind] <- highfreq[ind]+high.abun[i]

        }
        #websum <- sum(highfreq)
        diffh = (highfreq-colSums(web))/websum
        }
        if (is.null(high.abun)) pred_prop <- colSums(web)/websum else pred_prop <- highfreq/websum
        if (is.null(low.abun)) prey_prop <- rowSums(web)/websum else prey_prop <- lowfreq/websum
        n.pred <- length(pred_prop)
        n.prey <- length(prey_prop)
        pred_x <- 0
        pred_xold <- -1
        pred_versatz <- 0
        pred_y <- 1.5
        prey_x <- 0
        prey_xold <- -1
        prey_versatz <- 0
        prey_y <- 0.5
        if (length(colnames(web)) == 0)
            colnames(web) <- colnames(web, do.NULL = FALSE)
        if (length(rownames(web)) == 0)
            rownames(web) <- rownames(web, do.NULL = FALSE)
        if (!is.null(lablength))
            colnames(web) <- substr(colnames(web), 1, lablength)
        if (!is.null(lablength))
            rownames(web) <- substr(rownames(web), 1, lablength)

        par(mai = c(0.2, 0.2, 0.2, 0.2))
        pred_spacing =  (n.prey - 1)/(n.pred - 1)
        prey_spacing =   (n.pred - 1)/(n.prey -1)
        pred_spacing <- pred_spacing*spacing
        prey_spacing <- prey_spacing*spacing

        if (n.pred>n.prey) prey_spacing <- pred_spacing*(n.pred-1)/(n.prey-1) else pred_spacing<- prey_spacing*(n.prey-1)/(n.pred-1)
        
        if (!is.null(low.abun)) pred_spacing <- pred_spacing+sum(difff)/n.pred
        
        if (!is.null(high.abun)) prey_spacing <- prey_spacing+sum(diffh)/n.prey

        
        wleft = 0
        wright = (max(n.pred, n.prey)) * min(prey_spacing,pred_spacing) +1+max(sum(diffh),sum(difff))

        wup <- 2.6 + y_width + lab.space * 0.05 #we need some space for labels
        wdown <- 0.4 - y_width - lab.space * 0.05 #we need some space for labels

        plot(0, type = "n", xlim = range(wleft, wright), ylim = range(wdown/ybig,
            wup * ybig), axes = FALSE, xlab = "", ylab = "")
        pred_x = 0
        hoffset <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(colnames(web)[1], cex = 0.6)
        for (i in 1:n.pred) {
            rect(pred_x, pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = col.pred)
            #### coloured boxes at the end if highfreq is given
            if (!is.null(high.abun))
              {
              rect(pred_x + pred_prop[i]-diffh[i], pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = high.abun.col)
              }

            breite <- strwidth(colnames(web)[i], cex = 0.6 *
                labsize)
            links <- pred_x + pred_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- pred_x + pred_prop[i]/2 + breite/2
                hoffset <- 0
            }
            text(pred_x + pred_prop[i]/2, pred_y + y_width +
                hoehe + hoffset, colnames(web)[i], cex = 0.6 *
                labsize, offset = 0)
            pred_x <- pred_x + pred_prop[i] + pred_spacing
        }
        prey_x <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(rownames(web)[1], cex = 0.6)
        hoffset <- hoehe
        for (i in 1:n.prey) {
            rect(prey_x, prey_y - y_width, prey_x + prey_prop[i],
                prey_y + y_width, col = col.prey)
            #### coloured boxes at the end if lowfreq is given
            if (!is.null(low.abun))
              {
              rect(prey_x + prey_prop[i]-difff[i], prey_y - y_width, prey_x + prey_prop[i],
                prey_y + y_width, col = low.abun.col)
              }
            breite <- strwidth(rownames(web)[i], cex = 0.6 *
                labsize)
            links <- prey_x + prey_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- prey_x + prey_prop[i]/2 + breite/2
                hoffset <- hoehe
            }
            text(prey_x + prey_prop[i]/2, prey_y - y_width -
                hoffset, rownames(web)[i], cex = 0.6 * labsize,
                offset = 0)
            prey_x <- prey_x + prey_prop[i] + prey_spacing
        }
       pred_x <- 0
        zwischenweb <- web
        XYcoords <- matrix(ncol = 2, nrow = length(zwischenweb))
        for (i in 1:length(zwischenweb)) {
            XYcoords[i, ] <- which(zwischenweb == max(zwischenweb),
                arr.ind = TRUE)[1, ]
            zwischenweb[XYcoords[i, 1], XYcoords[i, 2]] <- -1
        }
        y1 <- pred_y - y_width
        y2 <- y1
        y3 <- prey_y + y_width
        y4 <- y3
        for (p in 1:sum(web > 0)) {
            i <- XYcoords[p, 1]
            j <- XYcoords[p, 2]
            if (j == 1 & i == 1)
                x1 <- 0
            else x1 <- (j - 1) * pred_spacing + cumsum(web)[(j -
                1) * nrow(web) + (i - 1)]/websum
            if (!is.null(high.abun) && j>1) x1 <- x1 +cumsum(diffh)[j-1]
            x2 <- x1 + web[i, j]/websum
            if (arrow=="up" || arrow=="both") {x2<-(x1+x2)/2; x1<-x2}
            tweb <- t(web)
            if (j == 1 & i == 1)
                x3 <- 0
            else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i -
                1) * nrow(tweb) + (j - 1)]/websum
            if (!is.null(low.abun) && i>1) x3 <- x3 +cumsum(difff)[i-1]
            x4 <- x3 + tweb[j, i]/websum
            if (arrow=="down" || arrow=="both") {x4<-(x3+x4)/2; x3<-x4}
            polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction)
        }



######web2 ######
###################


  web2<-as.matrix(web2) # to convert data.frames into matrix: needed for cumsum




####idee

#erstmal das netz umbauen, da leere Felder auch geplottet werden muessen (komplett unparasitiert!!!!!)




         for (i in 1:dim(web)[2] )
             {
             dn <- dimnames(web)[[2]][i]
             if (is.na(match(dn,dimnames(web2)[[1]])))
                {
                dummy <- matrix(rep(0,dim(web2)[2]),nrow=1)
                rownames(dummy) <- dn
                web2<-rbind(web2,dummy)

                }

             }
###umsortieren....nach der ordnung von web
         web2<-web2[order(dimnames(web2)[[1]]),]  #erst alphabetisch
         web2<-web2[rank(dimnames(web)[[2]]),] # dann die reihenfolge von web


         #websum <- sum(web2)
         difff <- diffh <-0 # if no abundances are set leave plotsize as before

        ###rearrange if lowfreq is set  # lowfreq is a named vector!!!

        dummy <- colSums(web)
        lowfreq = rowSums(web2)
        ### sortieren !!!!
        for (i in 1:length(dummy))
            {
            dummy[i] <- dummy[i] -lowfreq[which(names(lowfreq)==names(dummy[i]))]

            }



        low.abun2 <- dummy
        lowfreq=lowfreq+low.abun2

        difff = low.abun2/websum






        ###rearrange if highfreq is set  # lowfreq is a named vector!!!
        if (!is.null(high.abun2)) {

        highfreq = colSums(web2)

        dummy <- highfreq

        for (i in 1:length(high.abun2) )
        {
        ind <- which(names(high.abun2)[i] == names(dummy))

        highfreq[ind] <- highfreq[ind]+high.abun2[i]

        }
        #websum <- sum(highfreq)
        diffh = (highfreq-colSums(web2))/websum
        }






        if (is.null(high.abun2)) pred_prop <- colSums(web2)/websum else pred_prop <- highfreq/websum
        if (is.null(low.abun2)) prey_prop <- rowSums(web2)/websum else prey_prop <- lowfreq/websum
        n.pred <- length(pred_prop)
        n.prey <- length(prey_prop)
        pred_x <- 0
        pred_xold <- -1
        pred_versatz <- 0
        pred_y <- 2.5
        prey_x <- 0
        prey_xold <- -1
        prey_versatz <- 0
        prey_y <- 1.5
        if (length(colnames(web2)) == 0)
            colnames(web2) <- colnames(web2, do.NULL = FALSE)
        if (length(rownames(web2)) == 0)
            rownames(web2) <- rownames(web2, do.NULL = FALSE)
        if (!is.null(lablength2))
            colnames(web2) <- substr(colnames(web2), 1, lablength2)
        if (!is.null(lablength2))
            rownames(web2) <- substr(rownames(web2), 1, lablength2)

        prey_spacing =  pred_spacing
        pred_spacing =  (n.prey - 1)/(n.pred - 1)

        pred_spacing <- pred_spacing*spacing2
        #prey_spacing <- prey_spacing*spacing2

        if (n.pred<n.prey)  pred_spacing<- prey_spacing*(n.prey-1)/(n.pred-1) # else prey_spacing <- pred_spacing*(n.pred-1)/(n.prey-1)

        if (!is.null(low.abun2)) pred_spacing <- pred_spacing+sum(difff)/n.pred

        #if (!is.null(high.abun2)) prey_spacing <- prey_spacing+sum(diffh)/n.prey

        pred_x = 0
        hoffset <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(colnames(web2)[1], cex = 0.6)
        for (i in 1:n.pred) {
            rect(pred_x, pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = col.pred2)
            #### coloured boxes at the end if highfreq is given
            if (!is.null(high.abun2))
              {
              rect(pred_x + pred_prop[i]-diffh[i], pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = high.abun.col2)
              }

            breite <- strwidth(colnames(web2)[i], cex = 0.6 *
                labsize)
            links <- pred_x + pred_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- pred_x + pred_prop[i]/2 + breite/2
                hoffset <- 0
            }
            text(pred_x + pred_prop[i]/2, pred_y + y_width +
                hoehe + hoffset, colnames(web2)[i], cex = 0.6 *
                labsize, offset = 0)
            pred_x <- pred_x + pred_prop[i] + pred_spacing
        }
        prey_x <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(rownames(web2)[1], cex = 0.6)
        hoffset <- hoehe
        for (i in 1:n.prey) {
            rect(prey_x, prey_y - y_width, prey_x + prey_prop[i],
                prey_y + y_width, col = col.prey2)
            #### coloured boxes at the end if lowfreq is given
            if (!is.null(low.abun2))
              {
              rect(prey_x + prey_prop[i]-difff[i], prey_y , prey_x + prey_prop[i],
                prey_y + y_width, col = low.abun.col2)
              }
            breite <- strwidth(rownames(web)[i], cex = 0.6 *
                labsize)
            links <- prey_x + prey_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- prey_x + prey_prop[i]/2 + breite/2
                hoffset <- hoehe
            }
            text(prey_x + prey_prop[i]/2, prey_y - y_width -
                hoffset, rownames(web2)[i], cex = 0.6 * labsize,
                offset = 0)
            prey_x <- prey_x + prey_prop[i] + prey_spacing
        }
       pred_x <- 0
        zwischenweb <- web2
        XYcoords <- matrix(ncol = 2, nrow = length(zwischenweb))
        for (i in 1:length(zwischenweb)) {
            XYcoords[i, ] <- which(zwischenweb == max(zwischenweb),
                arr.ind = TRUE)[1, ]
            zwischenweb[XYcoords[i, 1], XYcoords[i, 2]] <- -1
        }
        y1 <- pred_y - y_width
        y2 <- y1
        y3 <- prey_y + y_width
        y4 <- y3
        if (sum(web2>0)) {
        for (p in 1:sum(web2 > 0)) {
            i <- XYcoords[p, 1]
            j <- XYcoords[p, 2]
            if (j == 1 & i == 1)
                x1 <- 0
            else x1 <- (j - 1) * pred_spacing + cumsum(web2)[(j -
                1) * nrow(web2) + (i - 1)]/websum
            if (!is.null(high.abun2) && j>1) x1 <- x1 +cumsum(diffh)[j-1]
            x2 <- x1 + web2[i, j]/websum
            if (arrow=="up" || arrow=="both") {x2<-(x1+x2)/2; x1<-x2}
            tweb <- t(web2)
            if (j == 1 & i == 1)
                x3 <- 0
            else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i -
                1) * nrow(tweb) + (j - 1)]/websum
            if (!is.null(low.abun2) && i>1) x3 <- x3 +cumsum(difff)[i-1]
            x4 <- x3 + tweb[j, i]/websum
            if (arrow=="down" || arrow=="both") {x4<-(x3+x4)/2; x3<-x4}
            polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = col.interaction2)
        }
        }


}


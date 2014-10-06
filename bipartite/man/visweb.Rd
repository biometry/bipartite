\encoding{UTF-8}
\name{visweb}
\alias{visweb}

\title{ Plotting function to visualize a bipartite food web }
\description{ This function draws a food web as a grid using a matrix. Colnames and rownames are used as labels and entries in the matrix are visualized by text and colours. It can also be used to plot bipartite webs in the style of Vázquez et al. (2009).

}

\usage{
visweb(web, type="nested",  prednames=TRUE, preynames=TRUE, labsize=1,
  plotsize=12, square="interaction", text="no", frame=NULL, textsize=1,
  textcol="red", pred.lablength=NULL, prey.lablength=NULL, clear=TRUE,
  xlabel="", ylabel="", boxes=TRUE, circles=FALSE, circle.col="black",
  circle.min=0.2, circle.max=2, outerbox.border="white", 
  outerbox.col="white", box.border="black", box.col="black", def.col="blue", 
  max.digits=4, NA.col="red")
}

\arguments{
  \item{web}{ A matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey. }
  \item{type}{ type changes the sorting of rows and coloumns of the web and can be\cr
            \bold{nested  }: (sorted by row/colSums)\cr
            \bold{diagonal}: (highest number of interactions appear along the diagonal, good for showing compartments)\cr
            \bold{none    }: (as is)\cr
  }
  \item{prednames}{ labels can be switched of by prednames=F}
  \item{preynames}{ labels can be switched of by preynames=F}
  \item{labsize}{ factor for size of labels }
  \item{plotsize}{  size of plot (length of width or height), depending on the dimension of the web in cm,  default is 12 cm.}
  \item{square}{ square is used to indicate number of interactions  and belonging to compartments by coloured grid cells\cr
            \bold{interaction}: (level of grey indicates the number of intercation, white means no interaction)\cr
            \bold{compartment}: (level of grey indicates belonging to the same compartment)\cr
            \bold{black      }: (black grid cells if number of interaction is bigger than one)\cr
            \bold{none       }: (no coloured squares are drawn)\cr
  }
  \item{text}{ number of interactions or belonging are plotted into each grid cell\cr 
            \bold{interaction}: (number of interactions are drawn)\cr
            \bold{compartment}: (belonging to same compartment indicated by capitel Letters)\cr
            \bold{none       }: (no text is shown)\cr
  }
  \item{frame}{ a frame is drawn around each compartment (frame=TRUE), best used with type=\dQuote{diagonal}  }
  \item{textsize}{ factor for size of text in squares, default=1 }
  \item{textcol}{  color of text in grid cells, default =\dQuote{red}}
  \item{pred.lablength}{length of predators (upper) labels that should be displayed}
  \item{prey.lablength}{length of prey (lower) labels that should be displayed}
  \item{clear}{ delete species with no interactions (compulsory done for \dQuote{nested} and \dQuote{diagonal}) }
  \item{xlabel}{label on the x-axis, make sure prey.lablength is set accordingly, default is empty}
    \item{ylabel}{label on the y-axis, make sure pred.lablength is set accordingly, default is empty}
  \item{boxes}{logical, if boxes should be drawn. Default is set to TRUE}
  \item{circles}{logical, if circles in a Vazquez et al. style should be drawn. Default is set to FALSE, size and colours of circles and background can be set by the following arguments}
  \item{circle.col}{Colour of circles, works only if circles=TRUE}
  \item{circle.min}{minimal size of circles, use to rescale circles appropriately, default is 0.2}
  \item{circle.max}{maximal size of circle, , use to rescale circles appropriately, default is 2}
  \item{outerbox.border}{Colour of outerbox border if option circles=TRUE}
  \item{outerbox.col}{Colour of background if option circles=TRUE}
  \item{box.border}{Colour of border of boxes, if option square=\dQuote{b} is set}
  \item{box.col}{Colour of boxes, if option square=\dQuote{b} is set}
  \item{def.col}{A user defined color vector for all ranks (starting from low to high) of occurring values in the network}  
  \item{max.digits}{defaults to 4.}
  \item{NA.col}{Which colour should be used for missing data (NAs)? Defaults to red.}
}

\value{
  A plot window with appropriate size according to the dimensions of the web.
}

\note{If labels don't fit, resize window by hand!}

\references{Vázquez, P.D., Chacoff, N.,P. and  Cagnolo, L. (2009) Evaluating multiple determinants of the structure of plant-animal mutualistic networks. \emph{Ecology}}

\author{ Bernd Gruber }

\seealso{ For a different plot on food web see \code{\link{plotweb}}. }

\examples{

data(Safariland)
visweb(Safariland)
visweb(Safariland, type="diagonal", square="compartment", text="none", 
	frame=TRUE)
visweb(Safariland, type="nested", text="compartment")
visweb(Safariland, circles=TRUE,  boxes=FALSE,  labsize=1, circle.max=3, 
	text="no")
visweb(Safariland,square="b",box.col="green",box.border="red")

#define your colors here,length has to be the numbers of different entries
cols <-0:(length(table(Safariland))-1) 
visweb(Safariland, square="defined", def.col=cols) 
}
\keyword{ package }


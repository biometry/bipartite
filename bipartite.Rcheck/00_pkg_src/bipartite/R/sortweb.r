sortweb<-function(web, sort.order="dec", sequence=NULL){
# sort a matrix by the following option
#      inc (increasing col and row totals)
#      dec (decreasing col and row totals)
#      seq (by a given sequence, which is a list of names for seq.pred and seq.prey)
if (sort.order=="inc") web <- web[order(rowSums(web),decreasing=FALSE), order(colSums(web),decreasing=FALSE) ]

if (sort.order=="dec") web <- web[order(rowSums(web),decreasing=TRUE), order(colSums(web),decreasing=TRUE) ]

if (sort.order=="seq"){
       if (is.null(colnames(web)) | is.null(rownames(web))) stop("This web is unnamed:can't sort it by name ...")
       if (is.null(sequence)) stop("Please give sequences as list (see ?sortweb).")
       web <- web[sequence$seq.lower, sequence$seq.higher]
   }
   
  return(web)
}


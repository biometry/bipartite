`empty` <-
function(web, count=FALSE){
    # simple helper function to do two things:
    # 1. clean up a web from empty columns/rows
    # 2. return how many cols and rows had to be removed
    # together with extinction, it counts how many species were affected by the
    # removal of a species, because these rows/cols are now empty
    # C.F. Dormann, 6 Mar 2007

    web[is.na(web)] <- 0

    if (NCOL(web)==1 | NROW(web)==1)
    {
        if (NCOL(web)==1 & NROW(web) != 1) {nr <- sum(web>0); nc <- 1}
        if (NROW(web)==1 & NCOL(web) != 1) {nc <- sum(web>0); nr <- 1}
	      if (NROW(web)==1 & NCOL(web) == 1) {nr <- 1; nc <- 1}
        out <- web[1:nr, 1:nc, drop=FALSE]
        if (count) attr(out, "empty") <- c("empty rows"=NROW(web)-nr, "empty columns"=NCOL(web)-nc)
        return(out)
    }

    cempty <- which(colSums(web)==0)
    rempty <- which(rowSums(web)==0)

    cind <- if (length(cempty)==0) 1:NCOL(web) else (1:NCOL(web))[-cempty]
    rind <- if (length(rempty)==0) 1:NROW(web) else (1:NROW(web))[-rempty]

#    if (min(dim(as.matrix(web[rind, cind])))==1)
#    { #problem "solved" here: R turns a one-column matrix into a vector, loosing thus
#      #the information if it was the rows or the columns that were zero
#      out <- matrix(web[rind, cind], nrow=length(rind), ncol=length(cind))
#    } else
    out <- web[rind, cind, drop=FALSE]

    if (count) attr(out, "empty") <- c("empty rows"=length(rempty), "empty columns"=length(cempty))

    return(out)

}


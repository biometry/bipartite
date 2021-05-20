#-- sortweb2: adjustments to matrix, outsorced from plotweb and adapted -------
# by Jochen Fründ, May2021  (but largely copied code from Bernd Gruber's plotweb)

# for development only
# web <- testweb

sortweb2 <- function(web, method="cca", empty=TRUE, sequence=NULL){
  if (empty) {web <- empty(web)} else {method <- "normal"}  # cannot use cca with 0-rows/cols (issue1)
  web <- as.matrix(web) # to convert data.frames into matrix: needed for cumsum
  
  # give rownames & colnames if missing
  colnames(web) <- colnames(web, do.NULL = FALSE)
  rownames(web) <- rownames(web, do.NULL = FALSE)
  
  # the following should be combined in one place with 0-rows/cols and "same value in each cell" conditions, which all share that cca isn't possible (but should only trigger if cca was asked for!) (issue1, issue2, issue3)
  if (NROW(web) == 1 | NCOL(web) ==1) {
  	# sequence <- NULL # I guess this is wrong and not needed?
  	method <- "normal" # cannot use cca with 1 species in a guild
  }

  meths <- c("normal", "cca")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) stop("Choose sorting-method: normal/cca.\n")
  if (length(unique(as.vector(web))) == 1) { # (issue3)
    meths.match <- 1
    warning("CCA-sorting does not work with same value in each cell. Uses method='normal' instead.")
  }
  
  # simplest sorting: as is (i.e. meths.match==1 and sequence=NULL)
  row.seq <- 1:nrow(web)
  col.seq <- 1:ncol(web)
  
  # option "cca" = the web is re-arranged by ordination (& separating compartments)
  if (meths.match==2){
    # Problem: cca sometimes doesn't get the compartments right!
    # Solution: Function "compart" returns a matrix with links assigned to compartments
    # So, we need to extract the compartments there and put them in sequence,
    # order by cca only within compartments
    co <- compart(web)
    if (co$n.compart>1){ # do the arrangement for each compartment separately
      row.seq <- NULL
      col.seq <- NULL
      for (m in 1:co$n.compart){
        comp.member <- which(abs(co$cweb)==m, arr.ind=TRUE)
        rs <- unique(comp.member[,1])
        cs <- unique(comp.member[,2])
        if (length(rs) < 3 | length(cs) < 3){
          row.seq <- c(row.seq, rs)
          col.seq <- c(col.seq, cs)
        } else { # works fine for webs with only one compartment
          ca <- cca(web[rs, cs])
          row.seq <- c(row.seq, rs[order(summary(ca)$sites[,1], decreasing=TRUE)])
          col.seq <- c(col.seq, cs[order(summary(ca)$species[,1], decreasing=TRUE)])
        }
      }
    } else {  # webs with 1 compartment; this may not be needed, as the code above likely works also for these
      ca <- cca(web)
      row.seq <- order(summary(ca)$sites[,1], decreasing=TRUE)
      col.seq <- order(summary(ca)$species[,1], decreasing=TRUE)
    }
  } # end for meths.match==2 condition; start of the "normal" plotting

  if (!is.null(sequence)) {
    if (meths.match==2) {warning("cca sorting is not used with given sequence")}
    row.seq <- sequence$seq.low[sequence$seq.low %in% rownames(web)]
    col.seq <- sequence$seq.high[sequence$seq.high %in% colnames(web)]
  }
  
  web <- web[row.seq, col.seq, drop=FALSE]  # now once for all methods
  return(web)
}

# Examples:
# sortweb2(testweb)
# sortweb2(testweb, sequence=sequence)
# sortweb2(testweb, empty=FALSE)

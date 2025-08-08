#-- sortweb2: adjustments to matrix, outsorced from plotweb and adapted -------
# by Jochen Fründ, May2021  (but largely copied code from Bernd Gruber's plotweb)


sortweb <- function(web,
                     sort.order = "dec",
                     empty = TRUE,
                     sequence = NULL) {
  if (empty) {
    web <- empty(web)
  }
  web <- as.matrix(web) # to convert data.frames into matrix: needed for cumsum

  # give rownames & colnames if missing
  colnames(web) <- colnames(web, do.NULL = FALSE)
  rownames(web) <- rownames(web, do.NULL = FALSE)

  # choose sort.order (method in plotwebr)
  methods <- c("normal", "ca", "sequence", "decreasing", "increasing")
  method.matched <- methods[pmatch(sort.order, methods)]
  if (is.na(method.matched)) stop("Choose one of the available sorting-methods (sort.order in sortweb)!\n")

  if (method.matched == "ca") {
    # all cases where ca wouldn't work now combined in one place
    if (NROW(web) == 1 | NCOL(web) == 1 | length(unique(as.vector(web))) == 1) {
      method.matched <- "normal"
    }
    if (any(rowSums(web) == 0, colSums(web) == 0)) {
      method.matched <- "normal"
      warning("cannot use ca with 0-rows/cols, using sort.order='normal' instead")
    }
  }

  # simplest sorting as start; for sort.order="normal", no other sorting will be applied
  row.seq <- 1:nrow(web)
  col.seq <- 1:ncol(web)

  # option "cca" = the web is re-arranged by ordination (& separating compartments)
  if (method.matched == "ca") {
    # Problem: ca sometimes doesn't get the compartments right!
    # Solution: Function "compart" returns a matrix with links assigned to compartments
    # So, we need to extract the compartments there and put them in sequence, sort by ca within
    co <- compart(web)
    # do the arrangement for each compartment separately
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
        ca <- ca(web[rs, cs])
        row.seq <- c(row.seq, rs[order(scores(ca)$sites[,1], decreasing=TRUE)])
        col.seq <- c(col.seq, cs[order(scores(ca)$species[,1], decreasing=TRUE)])
      }
    }
  } # end of ca method

  if (method.matched == "increasing") {
    web <- web[order(rowSums(web),decreasing=FALSE), order(colSums(web),decreasing=FALSE)]
  }

  if (method.matched == "decreasing") {
    web <- web[order(rowSums(web),decreasing=TRUE), order(colSums(web),decreasing=TRUE)]
  } 

  if (!is.null(sequence) & method.matched != "sequence") {
    # warning("giving a sorting sequence overrides other sorting options")
  }

  if (method.matched == "sequence"){ 
    if (is.null(sequence) | !is.list(sequence)) stop("Please give sequences as properly formatted list (see ?sortweb).")
    if (length(sequence) == 2){
      if (is.null(names(sequence))) {
        # for compatibility with old sortweb and for convenience, accept other or missing names
        names(sequence) <- c("seq.low", "seq.high")
      }
      if (all(names(sequence) == c("seq.lower", "seq.higher"))) {
        # for compatibility with old sortweb and for convenience, accept other or missing names
        names(sequence) <- c("seq.low", "seq.high")
      }
    }
    row.seq <- sequence$seq.low[sequence$seq.low %in% rownames(web)]
    col.seq <- sequence$seq.high[sequence$seq.high %in% colnames(web)]
  }

  web <- web[row.seq, col.seq, drop=FALSE]  # now called once for all methods (sort.orders)
  return(web)
}

# Examples:
# sortweb2(testweb)
# sortweb2(testweb, sequence=list(rownames(testweb)[3:2], colnames(testweb)[1:5]))
# sortweb2(testweb, empty=FALSE)

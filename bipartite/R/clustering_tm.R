`clustering_tm` <-
function(net,subsample=1,seed=NULL){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet)) {
    if(ncol(net)==3) {
      net <- as.tnet(net, type="weighted two-mode tnet")
    } else {
      net <- as.tnet(net, type="binary two-mode tnet")
    }
  }
  if(attributes(net)$tnet!="binary two-mode tnet" & attributes(net)$tnet!="weighted two-mode tnet")
    stop("Network not loaded properly")
  if(!is.null(seed))
    set.seed(as.integer(seed))

  weighted <- attributes(net)$tnet == "weighted two-mode tnet"

  # 1-paths i1+p1
  paths <- net
  if(subsample!=1) {
    if(subsample<1) {
      index <- sample.int(nrow(paths), round(nrow(paths)*subsample))
    } else {
      index <- sample.int(nrow(paths), as.integer(subsample))
    }
    index <- index[order(index)]
    paths <- paths[index,]
  }

  # ---------------------------------------------------------------------------
  # Incidence lookup for the 6-cycle test.
  # A 4-path i1-p1-i2-p2-i3 lies on a 6-cycle iff i1 and i3 share a primary node
  # other than p1 and p2. The original code answered that per row with
  #   ct <- c(net.list[[as.character(i1)]], net.list[[as.character(i3)]])
  #   any(duplicated(ct[ct != p1 & ct != p2]))
  # i.e. two list lookups, two string conversions and a duplicated() call for
  # every one of the (millions of) 4-paths. The same answer is one lookup in the
  # shared-partner matrix: i1~p1 and i3~p2 hold by construction, so the number of
  # shared partners excluding p1 and p2 is  S[i1,i3] - B[i3,p1] - B[i1,p2].
  # ---------------------------------------------------------------------------
  inode <- net[, 1]; pnode <- net[, 2]
  p.levels <- sort(unique(pnode))
  nI <- max(inode)
  B <- matrix(0L, nI, length(p.levels))
  B[cbind(inode, match(pnode, p.levels))] <- 1L
  S <- tcrossprod(B)                       # shared primary nodes per pair of secondary nodes
  closes.cycle <- function(i1, p1, p2, i3)
    (S[cbind(i1, i3)] - B[cbind(i3, match(p1, p.levels))] - B[cbind(i1, match(p2, p.levels))]) > 0

  # 2-paths i1+p1+i2
  if (weighted) {
    dimnames(paths)[[2]] <- c("i1","p1","w1")
    dimnames(net)[[2]] <- c("i2","p1","w2")
  } else {
    dimnames(paths)[[2]] <- c("i1","p1")
    dimnames(net)[[2]] <- c("i2","p1")
  }
  paths <- merge(paths, net, sort=FALSE)
  paths <- paths[paths[,"i1"] != paths[,"i2"],]
  # 3-paths i1+p1+i2+p2
  dimnames(net)[[2]] <- if (weighted) c("i2","p2","w3") else c("i2","p2")
  paths <- merge(paths, net, sort=FALSE)
  paths <- paths[paths[,"p1"] != paths[,"p2"],]

  # ---------------------------------------------------------------------------
  # 4-paths i1+p1+i2+p2+i3.
  # Expanding all of them at once is what made this function allocate several GB
  # on medium-sized webs (the 4-path count grows roughly cubically in web size).
  # Since only column sums of the per-path weights are needed, the expansion is
  # done in chunks and the sums accumulated.
  # ---------------------------------------------------------------------------
  edges.by.p <- split(seq_len(nrow(net)), match(net[, 2], p.levels))
  i1 <- paths[,"i1"]; i2 <- paths[,"i2"]; p1 <- paths[,"p1"]; p2 <- paths[,"p2"]
  if (weighted) { w1 <- paths[,"w1"]; w2 <- paths[,"w2"]; w3 <- paths[,"w3"] }
  rm(paths)

  cols <- c("bi","am","gm","ma","mi")
  denominator <- if (weighted) setNames(numeric(5), cols) else 0
  numerator   <- if (weighted) setNames(numeric(5), cols) else 0

  n3 <- length(i1)
  chunk <- max(1L, as.integer(2e5))
  starts <- if (n3 == 0) integer(0) else seq(1L, n3, by=chunk)
  for (st in starts) {
    sel <- st:min(st + chunk - 1L, n3)
    lst <- edges.by.p[match(p2[sel], p.levels)]
    lst[vapply(lst, is.null, logical(1))] <- list(integer(0))
    reps <- lengths(lst)
    if (sum(reps) == 0) next
    src <- rep.int(sel, reps)                      # which 3-path each 4-path came from
    e   <- unlist(lst, use.names=FALSE)            # which edge supplies p2--i3
    i3  <- net[e, 1]
    keep <- i1[src] != i3 & i2[src] != i3
    if (!any(keep)) next
    src <- src[keep]; e <- e[keep]; i3 <- i3[keep]
    hit <- closes.cycle(i1[src], p1[src], p2[src], i3)
    if (weighted) {
      W1 <- w1[src]; W2 <- w2[src]; W3 <- w3[src]; W4 <- net[e, 3]
      pw <- cbind(bi = rep.int(1, length(src)),
                  am = (W1 + W2 + W3 + W4)/4,
                  gm = sqrt(sqrt(W1 * W2 * W3 * W4)),
                  ma = pmax(W1, W2, W3, W4),
                  mi = pmin(W1, W2, W3, W4))
      denominator <- denominator + colSums(pw)
      numerator   <- numerator   + colSums(pw[hit, , drop=FALSE])
    } else {
      denominator <- denominator + length(src)
      numerator   <- numerator   + sum(hit)
    }
  }

  # Fraction
  return(numerator/denominator)
}

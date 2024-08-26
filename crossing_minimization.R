w <- function(web, x, j, l) {
  return(sum(web[j:l, x]))
}

c_u_v <- function(web, u, v) {
  ret <- 0
  web[v, ] * w(web, u, )
}

min_3_wolf <- function(web) {
  x <- 0
}

coarse_part <- function(web) {
  n0 <- nrow(web)
  n1 <- ncol(web)
  P <- vector("list", length = n0)
  for (u in 1:n1) {
    leftsum <- 0
    rightsum <- w(web, u, 2, n0)
    for (r in 0:(n0 - 1)) {
      if (leftsum >= rightsum) {
        break
      }
      leftsum <- leftsum + web[r + 1, u]
      rightsum <- rightsum - web[r + 2, u]
    }
    P[[r+1]] <- c(P[[r+1]], u)
  }
  return(P)
}
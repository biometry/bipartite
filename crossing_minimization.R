# Based on O. A. Çakıroḡlu, C. Erten, Ö. Karataş, and M. Sözdinler,
# “Crossing minimization in weighted bipartite graphs,”
# Journal of Discrete Algorithms, vol. 7, no. 4. Elsevier BV, pp. 439–452, Dec. 2009. 
# doi: 10.1016/j.jda.2008.08.003

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

fine_part_complete <- function(web, P) {
  Pi <- vector("list", length = length(P))
  for (r in seq_along(P)) {
    if (is.null(P[[r]])) {
      Pi[[r]] <- NULL
    } else {
      Pi[[r]] <- fine_part(web, P[[r]], r)
    }
  }
  return(Pi)
}

fine_part <- function(web, P, r) {
  if (length(P) == 1) {
    return(P)
  }
  n0 <- nrow(web)
  n_p <- length(P)
  pi_r <- vector(mode = "integer", length = n_p)
  half <- ceiling(n_p/2)
  Pr_1 <- P[1:half]
  Pr_2 <- P[(half + 1):n_p]
  Pr_1 <- fine_part(web, Pr_1, r)
  Pr_2 <- fine_part(web, Pr_2, r)
  u <- Pr_1[1]
  u_t <- 1
  v <- Pr_2[1]
  v_t <- 1
  for (t in seq_along(P)) {
    if (u_t > length(Pr_1)) {
      pi_r[t] <- v
      v_t <- v_t + 1
      v <- Pr_2[v_t]
    } else if (v_t > length(Pr_2)) {
      pi_r[t] <- u
      u_t <- u_t + 1
      u <- Pr_1[u_t]
    } else if ((w(web, v, 1, r) *  w(web, u, r + 1, n0)) <= (w(web, u, 1, r) * w(web, v, r + 1, n0))) {
      pi_r[t] <- u
      u_t <- u_t + 1
      u <- Pr_1[u_t]
    } else {
      pi_r[t] <- v
      v_t <- v_t + 1
      v <- Pr_2[v_t]
    }
  }
  return(pi_r)
}
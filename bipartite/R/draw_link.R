draw_link <- function(x1, x2, y1, y2, y3, y4, col,
                      curved = FALSE, horizontal = FALSE,
                      alpha = 0.4) {
  if (curved) { # Curved (polygon) links
    x <- matrix(c(x1, x2))
    xm <- 0.5 * (x1 + x2)
    A <- cbind(x^3, x^2, x, 1)
    A <- rbind(A, c(3 * x1^2, 2 * x1, 1, 0))
    A <- rbind(A, c(6 * xm, 2, 0, 0))
    b1 <- matrix(c(y1, y2, 0, 0))
    b2 <- matrix(c(y3, y4, 0, 0))
    poly1 <- solve(A, b1)
    poly2 <- solve(A, b2)
    xx <- seq(x1, x2, length.out = 500)
    curve1 <- poly1[1] * xx^3 + poly1[2] * xx^2 + poly1[3] * xx + poly1[4]
    curve2 <- poly2[1] * xx^3 + poly2[2] * xx^2 + poly2[3] * xx + poly2[4]
    if (horizontal) {
      polygon(c(xx, rev(xx)), c(curve1, rev(curve2)),
              col = adjustcolor(col, alpha.f = alpha),
              border = adjustcolor(col, alpha.f = alpha))
    } else {
      polygon(c(curve1, rev(curve2)), c(xx, rev(xx)),
              col = adjustcolor(col, alpha.f = alpha),
              border = adjustcolor(col, alpha.f = alpha))
    }
  } else { # Straight lines
    if (horizontal) { 
      polygon(c(x1, x1, x2, x2), c(y1, y3, y4, y2),
              col = adjustcolor(col, alpha.f = alpha),
              border = adjustcolor(col, alpha.f = alpha))
    } else {
      polygon(c(y1, y3, y4, y2), c(x1, x1, x2, x2),
              col = adjustcolor(col, alpha.f = alpha),
              border = adjustcolor(col, alpha.f = alpha))
    }
  }
}
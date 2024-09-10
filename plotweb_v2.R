plotweb_v2 <- function(web,
                       upper_abundances = NULL,
                       lower_abundances = NULL,
                       add_upper_abundances = NULL,
                       add_lower_abundances = NULL,
                       scaling = "relative",
                       font = NULL,
                       family = NULL,
                       srt = 0,
                       upper_italic = FALSE,
                       lower_italic = FALSE,
                       text_size = 1,#"auto",
                       x_lim = c(0, 1),
                       y_lim = c(0, 1),
                       lower_color = "black",
                       lower_border = "same",
                       upper_color = "black",
                       upper_border = "same",
                       horizontal = FALSE,
                       abbr_names = FALSE,
                       link_color = "lower",
                       link_border = "same",
                       link_alpha = 0.5,
                       style = "line",
                       arrow = "no",
                       spacing = "auto",
                       plot_axes = FALSE,
                       mar = c(1, 1, 1, 1),
                       mai = NULL) {
  #web <- sortweb2(web, sequence = NULL, empty = TRUE, sort.order = "cca")

  # Set the figure border via the mai or mar argument
  if (!is.null(mar) & !is.null(mai)) {
    warning("Both mar and mai have been set to values other than NULL. ",
            "This leads to mai overriding mar.")
  }
  if (!is.null(mar)){
    par(mar = mar)
  }
  if (!is.null(mai)){
    par(mai = mai)
  }

  r_names <- rownames(web)
  c_names <- colnames(web)

  nr <- nrow(web)
  nc <- ncol(web)

  # Recycle the color vectors
  if (length(upper_color) < nc) {
    upper_color <- rep_len(upper_color, nc)
  }
  if (length(lower_color) < nr) {
    lower_color <- rep_len(lower_color, nr)
  }

  # Check whether the color vector is a named vector
  if (!is.null(names(upper_color))) {
    # Sort the upper color vector the same way as the web
    upper_color <- upper_color[c_names]
  }

  if (abbr_names == "numbers") {
    r_names <- (nc + 1):(nc + nr)
    c_names <- 1:nc
    r_full_names <- rownames(web)
    c_full_names <- colnames(web)
  } else if (abbr_names == "letters") {
    r_names <- letters[(nc + 1):(nc + nr)]
    c_names <- letters[1:nc]
    r_names <- paste("lower_", r_names)
    c_names <- paste("upper_", c_names)
    r_full_names <- rownames(web)
    c_full_names <- colnames(web)
  }

  if (text_size == "auto") {
    total_str_h <- sum(strheight(c_names, srt = srt))
    text_size <- 1 / total_str_h
  }

  c_m_t_width <- max(strwidth(c_names, srt = srt,
                              cex = text_size, units = "inches"))
  r_m_t_width <- max(strwidth(r_names, srt = srt,
                              cex = text_size, units = "inches"))

  mai <- par()$mai

  if (horizontal) {
    mai[2] <- mai[2] + c_m_t_width
    mai[4] <- mai[4] + r_m_t_width
  } else {
    mai[1] <- mai[1] + c_m_t_width
    mai[3] <- mai[3] + r_m_t_width
  }

  par(mai = mai)

  if (horizontal) {
    plot(0, type = "n",
         xlim = x_lim,
         ylim = c(0, 1),
         axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  } else {
    plot(0, type = "n",
         ylim = y_lim,
         xlim = c(0, 1),
         axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }

  str_h <- 1.1 * strheight(c_names[1], srt = srt, cex = text_size)

  if (length(spacing) == 2) {
    c_space <- spacing[1] / (nc - 1)
    r_space <- spacing[2] / (nr - 1)
  } else if (length(spacing) == 1) {
    if (is.numeric(spacing)) {
      spacing <- c(spacing, spacing)
      c_space <- spacing[1] / (nc - 1)
      r_space <- spacing[2] / (nr - 1)
    }else if (spacing == "auto")  {
      space <- max(nc, nr) * str_h
      if (space > 1) {
        warning(paste("Text size is too large for auto spacing.",
                      "Decrease size or increase figure height for better results."))
      }
      spacing <- c(space, space)
      r_space <- space / (nr - 1)
      c_space <- space / (nc - 1)
    }
  } else {
    stop("Invalid spacing option.")
  }


  if (is.null(upper_abundances)) {
    upper_abundances <- colSums(web)
  }
  if (!is.null(add_upper_abundances)) {
    # Sort the additional upper abundances according to the column-order in the web
    add_upper_abundances <- add_upper_abundances[colnames(web)]
    upper_abundances <- c(rbind(upper_abundances, add_upper_abundances))
    #c_space <- c(0, c_space)
    upper_color <- rep(upper_color, each = 2)
  }

  if (is.null(lower_abundances)) {
    lower_abundances <- rowSums(web)
  }
  if (!is.null(add_lower_abundances)) {
    # Sort the additional lower abundances according to the row-order in the web
    add_lower_abundances <- add_lower_abundances[rownames(web)]
    # Merge the abundances and additional abundances vector
    # by alternating indices
    lower_abundances <- c(rbind(lower_abundances, add_lower_abundances))
    #r_space <- c(0, r_space)
    lower_color <- rep(lower_color, each = 2)
  }

  if (scaling == "relative") {
    c_prop_sizes <- (1 - spacing[1]) * upper_abundances / sum(upper_abundances)
    r_prop_sizes <- (1 - spacing[2]) * lower_abundances / sum(lower_abundances)
  } else if (scaling == "absolute") {
    max_abundances <- max(sum(upper_abundances), sum(lower_abundances))
    c_prop_sizes <- (1 - spacing[1]) * upper_abundances / max_abundances
    r_prop_sizes <- (1 - spacing[2]) * lower_abundances / max_abundances
    c_space <- (1 - sum(c_prop_sizes)) / (nc - 1)
    r_space <- (1 - sum(r_prop_sizes)) / (nr - 1)
  }

  if (!is.null(add_upper_abundances)) {
    c_space <- c(0, c_space)
    nc <- nc * 2
  }
  if (!is.null(add_lower_abundances)) {
    r_space <- c(0, r_space)
    nr <- nr * 2
  }

  c_xl <- c(0, cumsum(c_prop_sizes[-nc] + c_space))
  c_xr <- c(c_prop_sizes[1],
            c_prop_sizes[1] + cumsum(c_prop_sizes[-1] + c_space))

  r_xl <- c(0, cumsum(r_prop_sizes[-nr] + r_space))
  r_xr <- c(r_prop_sizes[1],
            r_prop_sizes[1] + cumsum(r_prop_sizes[-1] + r_space))

  if (!is.null(add_upper_abundances)) {
    c_tx <- (c_xl[seq(1, nc, 2)] + c_xr[seq(2, nc, 2)]) / 2
  } else{
    c_tx <- (c_xl + c_xr) / 2
  }
  if (!is.null(add_lower_abundances)) {
    r_tx <- (r_xl[seq(1, nr, 2)] + r_xr[seq(2, nr, 2)]) / 2
  } else {
    r_tx <- (r_xl + r_xr) / 2
  }

  if (upper_border == "same") {
    upper_border <- upper_color
  }
  if (lower_border == "same") {
    lower_border <- lower_color
  }

  c_m_t_width <- max(strwidth(c_names, srt = srt, cex = text_size))
  r_m_t_width <- max(strwidth(r_names, srt = srt, cex = text_size))

  if (upper_italic) {
    c_names <- lapply(c_names, function(x) bquote(italic(.(x))))
    c_names <- as.expression(c_names)
  }
  if (lower_italic) {
    r_names <- lapply(r_names, function(x) bquote(italic(.(x))))
    r_names <- as.expression(r_names)
  }

  # Draw the boxes and labels either horizontal or vertical to each other.
  if (horizontal) {
    rect(0, c_xl, 0.1, c_xr, col = upper_color, border = upper_border)
    rect(0.9, r_xl, 1, r_xr, col = lower_color, border = lower_border)
    text(-0.01, c_tx, c_names, adj = c(1, 0.5), srt = srt, cex = text_size,
         xpd = TRUE, font = font, family = family)
    text(1.01, r_tx, r_names, adj = c(0, 0.5), srt = srt, cex = text_size,
         xpd = TRUE, font = font, family = family)
  } else {
    rect(r_xl, 0, r_xr, 0.1, col = lower_color, border = lower_border)
    rect(c_xl, 0.9, c_xr, 1, col = upper_color, border = upper_border)
    text(r_tx, -0.01, r_names, adj = c(1, 0.5), cex = text_size,
         xpd = TRUE, srt = 90 + srt, font = font, family = family)
    # text(r_tx, 1.0 + strwidth(r_names, srt = srt) / 2, r_names,
    #      adj = c(.5, 0.5), cex = text_size, xpd = TRUE, srt = 90 + srt,
    #      font = font, family = family)
    text(c_tx, 1.01, c_names,
         adj = c(0, 0.5), cex = text_size, xpd = TRUE, srt = 90 + srt,
         font = font, family = family)
  }

  # Interactions
  if (!is.null(add_lower_abundances) && !is.null(add_upper_abundances)) {
    web.df <- data.frame(row = rep(seq(1, nr, 2), nc/2),
                         col = rep(seq(1, nc, 2), each = nr/2),
                         weight = c(web))
  } else if (!is.null(add_lower_abundances)) {
    web.df <- data.frame(row = rep(seq(1, nr, 2), nc),
                         col = rep(1:nc, each = nr/2),
                         weight = c(web))
  } else if (!is.null(add_upper_abundances)) {
    web.df <- data.frame(row = rep(1:nr, nc/2),
                         col = rep(seq(1, nc, 2), each = nr),
                         weight = c(web))
  } else {
    web.df <- data.frame(row = rep(1:nr, nc),
                         col = rep(1:nc, each = nr),
                         weight = c(web))
  }
  web.df <- web.df[web.df$weight > 0, ]

  # x-coordinates of interactions: tl = topleft, etc
  web.df[, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")] <- NA

  # low coordinates for interactions (in order of the web.df)
  for (i in unique(web.df$row)) { # for i in lower species
    # i <- 3
    links.i <- web.df[web.df$row == i, ]
    relpos <- cumsum(links.i$weight) / sum(links.i$weight)
    coords.int.low <- (r_xl[i] + relpos * (r_xr[i] - r_xl[i]))
    web.df[web.df$row == i, "xcoord.bl"] <- c(r_xl[i], coords.int.low[-nrow(links.i)])
    web.df[web.df$row == i, "xcoord.br"] <- c(coords.int.low)    
    if (arrow %in% c("down.center", "both.center")) {
      web.df[web.df$row == i, "xcoord.bl"] <- web.df[web.df$row == i, "xcoord.br"] <- mean(c(r_xl[i], r_xr[i]))
    }
  }
  if (arrow %in% c("down", "both")) {
    web.df[, "xcoord.bl"] <- web.df[, "xcoord.br"] <- rowMeans(web.df[, c("xcoord.bl", "xcoord.br")])
  }

  # high coordinates for interactions (in order of the web.df)
  for (j in unique(web.df$col)) { # for j in higher species
    links.j <- web.df[web.df$col == j, ]
    relpos <- cumsum(links.j$weight) / sum(links.j$weight)
    coords.int.high <- (c_xl[j] + relpos * (c_xr[j] - c_xl[j]))
    web.df[web.df$col == j, "xcoord.tl"] <- c(c_xl[j], coords.int.high[-nrow(links.j)])
    web.df[web.df$col == j, "xcoord.tr"] <- c(coords.int.high)    
    if (arrow %in% c("up.center", "both.center")) {
      web.df[web.df$col == j, "xcoord.tl"] <- web.df[web.df$col == j, "xcoord.tr"] <- mean(c(c_xl[j], c_xr[j]))
    }
  }
  if (arrow %in% c("up", "both")) {
    web.df[, "xcoord.tl"] <- web.df[, "xcoord.tr"] <- rowMeans(web.df[, c("xcoord.tl", "xcoord.tr")])
  }
  for (linki in order(-web.df$weight)) {
    link <- web.df[linki, ]
    if (horizontal) {
      x1 <- 0.1
      x2 <- 0.9
    } else {
      x1 <- 0.9
      x2 <- 0.1
    }
    y1 <- link$xcoord.tl
    y2 <- link$xcoord.bl
    y3 <- link$xcoord.tr
    y4 <- link$xcoord.br
    if (link_color == "lower") {
      l_col <- lower_color[link$row]
    } else if (link_color == "upper") {
      l_col <- upper_color[link$col]
    } else {
      l_col <- link_color
    }
    draw_link(x1, x2, y1, y2, y3, y4, l_col,
              horizontal = horizontal,
              alpha = link_alpha,
              style = style)
  }
}
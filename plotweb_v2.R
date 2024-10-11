plotweb_v2 <- function(web,
                       sorting = "none",
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
                       text_size = 1,
                       box_size = 0.1,
                       x_lim = c(0, 1),
                       y_lim = c(0, 1),
                       lower_color = "black",
                       lower_border = "same",
                       lower_add_color = "red",
                       upper_color = "black",
                       upper_border = "same",
                       upper_add_color = "red",
                       text_color = "black",
                       horizontal = FALSE,
                       #abbr_names = FALSE,
                       link_color = "lower",
                       link_border = "same",
                       link_alpha = 0.5,
                       style = "line",
                       arrow = "no",
                       spacing = "auto",
                       plot_axes = FALSE,
                       add = FALSE,
                       mar = c(1, 1, 1, 1),
                       mai = NULL) {
  # Set the figure border via the mai or mar argument
  if (!is.null(mar) & !is.null(mai)) {
    warning("Both mar and mai have been set to values other than NULL. ",
            "This leads to mai overriding mar.")
  }

  stopifnot(is.logical(upper_italic),
            is.logical(lower_italic),
            is.logical(horizontal),
            is.logical(plot_axes),
            (is.numeric(text_size) && text_size > 0))

  if (!is.null(mar)) {
    par(mar = mar)
  }
  if (!is.null(mai)) {
    par(mai = mai)
  }

  if (!is.null(rownames(web))) {
    r_names <- rownames(web)
  } else {
    r_names <- 1:nrow(web)
  }
  if (!is.null(rownames(web))) {
    c_names <- colnames(web)
  } else {
    c_names <- 1:ncol(web)
  }

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

  # if (upper_add_color == "same") {
  #   upper_color <- rep(upper_color, each = 2)
  # } else if (!is.null(names(upper_add_color))) {
  #   stopifnot(length(upper_add_color) == nc)
  #   if (!setequal(names(upper_add_color), c_names)) {
  #     stop("Names of upper_add_color does not match upper species names.")
  #   }
  #   upper_add_color <- upper_add_color[colnames(web1)]
  #   upper_color <- c(rbind(upper_color, upper_add_color))
  # } else if (length(upper_add_color) < nc) {
  #   upper_add_color <- rep_len(upper_color, nc)
  #   upper_color <- c(rbind(upper_color, upper_add_color))
  # }
  
  # if (lower_add_color == "same") {
  #   lower_color <- rep(lower_color, each = 2)
  # } else if (!is.null(names(lower_add_color))) {
  #   stopifnot(length(lower_add_color) == nr)
  #   lower_color <- c(rbind(lower_color, lower_add_color))
  # } else {
  #   lower_color <- c(rbind(lower_color, lower_add_color))
  # }

  ## TODO: Implement this correctly and add an option to create a legend
  # if (abbr_names == "numbers") {
  #   r_names <- (nc + 1):(nc + nr)
  #   c_names <- 1:nc
  #   r_full_names <- rownames(web)
  #   c_full_names <- colnames(web)
  # } else if (abbr_names == "letters") {
  #   r_names <- letters[(nc + 1):(nc + nr)]
  #   c_names <- letters[1:nc]
  #   r_names <- paste("lower_", r_names)
  #   c_names <- paste("upper_", c_names)
  #   r_full_names <- rownames(web)
  #   c_full_names <- colnames(web)
  # }


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
    mai[3] <- mai[3] + c_m_t_width
    mai[1] <- mai[1] + r_m_t_width
  }

  par(mai = mai)

  if (horizontal) {
    if (add == FALSE) {
      plot(0, type = "n",
          xlim = x_lim,
          ylim = y_lim,
          axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    }
    space_size <- y_lim[2] - y_lim[1]
    space_start <- y_lim[1]
  } else {
    if (add == FALSE) {
      plot(0, type = "n",
          ylim = y_lim,
          xlim = x_lim,
          axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    }
    space_size <- x_lim[2] - x_lim[1]
    space_start <- x_lim[1]
  }

  print(space_size)

  # if (horizontal) {
  #   str_h <- 1.1 * strheight(c_names[1], srt = srt, cex = text_size)
  # } else {
  #   str_h <- 1.1 * strwidth(c_names[1], srt = srt, cex = text_size)
  # }

  theta <- srt * pi / 180
  print(theta)

  # Get the unrotated height and width of the text
  H <- strheight(c_names[1], srt = srt, cex = text_size)
  print(H)
  W <- strwidth(c_names[1], srt = srt, cex = text_size)
  print(W)

  # Calculate the effective rotated height and width
  str_h <- abs(H * cos(theta)) + abs(W * sin(theta))

  print(str_h)

  if (length(spacing) == 2) {
    c_space <- space_size * spacing[1] / (nc - 1)
    r_space <- space_size * spacing[2] / (nr - 1)
  } else if (length(spacing) == 1) {
    if (is.numeric(spacing)) {
      spacing <- c(spacing, spacing)
      c_space <- space_size * spacing[1] / (nc - 1)
      r_space <- space_size * spacing[2] / (nr - 1)
    } else if (spacing == "auto")  {
      space <- space_size * max(nc, nr) * str_h
      if (space > 1) {
        stop(paste("Text size is too large for auto spacing.",
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
    c_prop_sizes <- (space_size - spacing[1]) * upper_abundances / sum(upper_abundances)
    r_prop_sizes <- (space_size - spacing[2]) * lower_abundances / sum(lower_abundances)
  } else if (scaling == "absolute") {
    max_abundances <- max(sum(upper_abundances), sum(lower_abundances))
    c_prop_sizes <- (space_size - spacing[1]) * upper_abundances / max_abundances
    r_prop_sizes <- (space_size - spacing[2]) * lower_abundances / max_abundances
    c_space <- (space_size - sum(c_prop_sizes)) / (nc - 1)
    r_space <- (space_size - sum(r_prop_sizes)) / (nr - 1)
  }

  if (!is.null(add_upper_abundances)) {
    c_space <- c(0, c_space)
    nc <- nc * 2
  }
  if (!is.null(add_lower_abundances)) {
    r_space <- c(0, r_space)
    nr <- nr * 2
  }

  c_xl <- c(space_start, cumsum(c_prop_sizes[-nc] + c_space))
  c_xr <- c(c_prop_sizes[1],
            c_prop_sizes[1] + cumsum(c_prop_sizes[-1] + c_space))

  r_xl <- c(space_start, cumsum(r_prop_sizes[-nr] + r_space))
  r_xr <- c(r_prop_sizes[1],
            r_prop_sizes[1] + cumsum(r_prop_sizes[-1] + r_space))

  if (!is.null(add_upper_abundances)) {
    c_tx <- (c_xl[seq(1, nc, 2)] + c_xr[seq(2, nc, 2)]) / 2
  } else {
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

  if (length(box_size) == 1) {
    upper_box_size <- box_size
    lower_box_size <- box_size
  } else if (length(box_size) == 2) {
    upper_box_size <- box_size[1]
    lower_box_size <- box_size[2]
  }

  # Draw the boxes and labels either horizontal or vertical to each other.
  if (horizontal) {
    x_start <- x_lim[1]
    x_end <- x_lim[2]
    rect(x_start, c_xl, x_start + upper_box_size, c_xr,
         col = upper_color, border = upper_border)
    rect(x_end - lower_box_size, r_xl, x_end, r_xr,
         col = lower_color, border = lower_border)
    text(x_start - 0.01, c_tx, c_names, adj = c(1, 0.5), srt = srt,
         cex = text_size, xpd = TRUE, font = font, family = family)
    text(x_end + 0.01, r_tx, r_names, adj = c(0, 0.5), srt = srt,
         cex = text_size, xpd = TRUE, font = font, family = family)
  } else {
    y_start <- y_lim[1]
    y_end <- y_lim[2]
    rect(r_xl, y_start, r_xr, y_start + lower_box_size, col = lower_color, border = lower_border)
    rect(c_xl, y_end - upper_box_size, c_xr, y_end, col = upper_color, border = upper_border)
    text(r_tx, y_start - 0.01, r_names, adj = c(1, 0.5), cex = text_size,
         xpd = TRUE, srt = srt + 90, font = font, family = family)
    # text(r_tx, 1.0 + strwidth(r_names, srt = srt) / 2, r_names,
    #      adj = c(.5, 0.5), cex = text_size, xpd = TRUE, srt = 90 + srt,
    #      font = font, family = family)
    text(c_tx, y_end + 0.01, c_names,
         adj = c(0, 0.5), cex = text_size, xpd = TRUE, srt = srt + 90,
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
      x1 <- x_start + upper_box_size
      x2 <- x_end - lower_box_size
    } else {
      x1 <- y_end - lower_box_size#0.9
      x2 <- y_start + upper_box_size#0.1
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

  ## TODO: Implement legend plotting
  # if (legend == TRUE) {
  #   plot()
  # }
}
draw_link <- function(x1, x2, y1, y2, y3, y4, col,
                      curved = FALSE,
                      horizontal = FALSE,
                      alpha = 0.4,
                      border = "same") {
  if (border == "same") {
    border <- adjustcolor(col, alpha.f = alpha)
  }
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
              border = border)
    } else {
      polygon(c(curve1, rev(curve2)), c(xx, rev(xx)),
              col = adjustcolor(col, alpha.f = alpha),
              border = border)
    }
  } else { # Straight lines
    if (horizontal) { 
      polygon(c(x1, x1, x2, x2), c(y1, y3, y4, y2),
              col = adjustcolor(col, alpha.f = alpha),
              border = border)
    } else {
      polygon(c(y1, y3, y4, y2), c(x1, x1, x2, x2),
              col = adjustcolor(col, alpha.f = alpha),
              border = border)
    }
  }
}

plotweb <- function(web,
                    sorting = "normal",
                    empty = FALSE,
                    higher_abundances = NULL,
                    lower_abundances = NULL,
                    add_higher_abundances = NULL,
                    add_lower_abundances = NULL,
                    higher_labels = NULL,
                    lower_labels = NULL,
                    scaling = "relative",
                    font = NULL,
                    family = NULL,
                    srt = 0,
                    higher_italic = FALSE,
                    lower_italic = FALSE,
                    text_size = "auto",
                    spacing = 0.3,
                    box_size = 0.1,
                    x_lim = c(0, 1),
                    y_lim = c(0, 1),
                    lab_distance = 0.05,
                    lower_color = "black",
                    lower_border = "same",
                    lower_add_color = "red",
                    lower_text_color = "black",
                    higher_color = "black",
                    higher_border = "same",
                    higher_add_color = "red",
                    higher_text_color = "black",
                    horizontal = FALSE,
                    link_color = "higher",
                    link_border = "same",
                    link_alpha = 0.5,
                    curved_links = FALSE,
                    arrow = "no",
                    plot_axes = FALSE,
                    add = FALSE,
                    mar = c(1, 1, 1, 1),
                    mai = NULL) {
  # Set the figure border via the mai or mar argument
  if (!is.null(mar) & !is.null(mai)) {
    warning("Both mar and mai have been set to values other than NULL. ",
            "This leads to mai overriding mar.")
  }

  stopifnot(is.logical(higher_italic),
            is.logical(lower_italic),
            is.logical(horizontal),
            is.logical(plot_axes),
            is.numeric(box_size),
            is.numeric(lab_distance),
            is.logical(curved_links),
            is.logical(add),
            ((is.numeric(text_size) && text_size > 0) | text_size == "auto"))

  # Update the user defined margin in either rows or inches
  if (!is.null(mar)) {
    par(mar = mar)
  }
  if (!is.null(mai)) {
    par(mai = mai)
  }
  # Extract newly set margin in inches
  mai <- par()$mai

  # Sort the web according to the user defined method
  web <- sortweb(web, sort.order = sorting, empty = empty)

  # Extract the number of rows and columns
  nr <- nrow(web)
  nc <- ncol(web)

  # Extract the row and column names if existing
  # otherwise generate number sequence instead
  if (!is.null(rownames(web))) {
    r_names <- rownames(web)
  } else {
    r_names <- 1:nr
  }
  if (!is.null(rownames(web))) {
    c_names <- colnames(web)
  } else {
    c_names <- 1:nc
  }

  # Reverse the web in horizontal mode,
  # so that the first species are plotted at the top
  if (horizontal) {
    r_names <- rev(r_names)
    c_names <- rev(c_names)
    web <- web[r_names, c_names]
  }

  # Extract the higher and lower label texts
  # if they are defined. Otherwise take the
  # column and row names of the web.
  if (is.null(higher_labels)) {
    higher_labels <- c_names
  } else if (length(higher_labels) > 1) {
    stopifnot(length(higher_labels) == nc)
    # If the higher_labels is a named list
    # reorder it according to the web column order
    if (!is.null(names(higher_labels))) {
      higher_labels <- higher_labels[c_names]
      higher_labels <- as.expression(higher_labels)
    }
  } else if (higher_labels == FALSE) {
    higher_labels <- rep_len("", nc)
  } else {
    stop("Unrecognized value for argument higher_labels: ",
         higher_labels)
  }

  # Lower label texts
  if (is.null(lower_labels)) {
    lower_labels <- r_names
  } else if (length(lower_labels) > 1) {
    stopifnot(length(lower_labels) == nr)
    # If the lower_labels is a named list
    # reorder it according to the web row order
    if (!is.null(names(lower_labels))) {
      lower_labels <- lower_labels[r_names]
      lower_labels <- as.expression(lower_labels)
    }
  } else if (lower_labels == FALSE) {
    lower_labels <- rep_len("", nr)
  } else {
    stop("Unrecognized value for argument lower_labels: ",
         lower_labels)
  }

  # Recycle the color vectors
  if (length(higher_color) < nc) {
    higher_color <- rep_len(higher_color, nc)
  }
  if (length(lower_color) < nr) {
    lower_color <- rep_len(lower_color, nr)
  }

  # Check whether the color vector is a named vector
  if (!is.null(names(higher_color))) {
    # Sort the higher color vector the same way as the web
    higher_color <- higher_color[c_names]
  }
  # Check whether the color vector is a named vector
  if (!is.null(names(lower_color))) {
    # Sort the higher color vector the same way as the web
    lower_color <- lower_color[r_names]
  }

  # lab_distance contains the distance between.
  # the abundance boxes and the labels.
  # If 2 values are given as a vector,
  # split the label distance to higher and lower.
  if (length(lab_distance) == 1) {
    u_lab_distance <- lab_distance
    l_lab_distance <- lab_distance
  } else if (length(lab_distance) == 2) {
    u_lab_distance <- lab_distance[1]
    l_lab_distance <- lab_distance[2]
  }

  # Clean up text rotation inputs
  srt <- srt %% 360

  if (horizontal) {
    theta <- (srt + 90) %% 360
  } else {
    theta <- srt
  }

  theta_pi <- theta * pi / 180

  ## TODO: Documentation + change 0.2 to actual spacing values
  if (text_size == "auto") {
    # If the bipartite plot is added to an existing plot
    # calculate the size of the plotting area in that plot in inches.
    if (add == TRUE) {
      dev_width <- diff(grconvertX(x_lim), from = "user", to = "inches")
      dev_height <- diff(grconvertY(y_lim), from = "user", to = "inches")
    } else {
      # Get the size of the plotting device in inches
      dev_size <- dev.size("in")

      c_height_1 <- strwidth(higher_labels[1], units = "inches")
      c_height_1 <- cos(theta_pi) * c_height_1
      r_height_1 <- strwidth(lower_labels[1], units = "inches")
      r_height_1 <- cos(theta_pi) * r_height_1

      c_height_n <- strwidth(higher_labels[nc], units = "inches")
      c_height_n <- cos(theta_pi) * c_height_n
      r_height_n <- strwidth(lower_labels[nr], units = "inches")
      r_height_n <- cos(theta_pi) * r_height_n

      max_height_1 <- max(c_height_1, -r_height_1)
      max_height_n <- max(-c_height_n, r_height_n)
      # Substract the margings from the total device size
      # to get the actual plotting area.
      dev_width <- dev_size[1] - mai[2] - mai[4] - max_height_1 - max_height_n
      dev_height <- dev_size[2] - mai[1] - mai[3] - max_height_1 - max_height_n
    }
    if (horizontal) {
      dev_max <- spacing * dev_height
    } else {
      dev_max <- spacing * dev_width
    }
    a <- theta_pi - pi / 2
    # TODO: Trust the math
    min_distance_at_angle <- tan(a)^2 * abs(cos(a)) + abs(cos(a))
    # if (theta < 45 || theta > 315 || (theta > 135 && theta < 225)) {
    if (min_distance_at_angle > max(nchar(c(higher_labels, lower_labels)))) {
      sum_str_h <- sum(strwidth(higher_labels, units = "inches"))
      sum_str_l <- sum(strwidth(lower_labels, units = "inches"))
    } else {
      sum_str_h <- 1.25 * min_distance_at_angle * sum(strheight(higher_labels, units = "inches"))
      sum_str_l <- 1.25 * min_distance_at_angle * sum(strheight(lower_labels, units = "inches"))
    }
    if (sum_str_h > dev_max || sum_str_l > dev_max) {
      text_size <- min(dev_max / sum_str_h, dev_max / sum_str_l)
    } else {
      text_size <- 1
    }
  }

  if (add == FALSE) {
    # Get the maximal width of the higher and lower labels
    # in inches to set the margin accordingly
    if (theta != 0 && theta != 180) {
      c_m_t_width <- max(strwidth(higher_labels,
                                  cex = text_size,
                                  units = "inches"))
      r_m_t_width <- max(strwidth(lower_labels,
                                  cex = text_size,
                                  units = "inches"))

      c_t_width <- abs(c_m_t_width * sin(theta_pi))
      r_t_width <- abs(r_m_t_width * sin(theta_pi))
    } else { # If the labels are horizontal use their height.
      c_t_width <- max(strheight(higher_labels,
                                 cex = text_size,
                                 units = "inches"))
      r_t_width <- max(strheight(lower_labels,
                                 cex = text_size,
                                 units = "inches"))
    }

    c_height_1 <- strwidth(higher_labels[1], cex = text_size, units = "inches")
    c_height_1 <- cos(theta_pi) * c_height_1
    r_height_1 <- strwidth(lower_labels[1], cex = text_size, units = "inches")
    r_height_1 <- cos(theta_pi) * r_height_1

    c_height_n <- strwidth(higher_labels[nc], cex = text_size, units = "inches")
    c_height_n <- cos(theta_pi) * c_height_n
    r_height_n <- strwidth(lower_labels[nr], cex = text_size, units = "inches")
    r_height_n <- cos(theta_pi) * r_height_n

    max_height_1 <- max(c_height_1, -r_height_1)
    max_height_n <- max(-c_height_n, r_height_n)

    # Add the label widths of the labels so they always fit on the plot
    if (horizontal) {
      mai[1] <- mai[1] + max_height_1
      mai[2] <- mai[2] + c_t_width + u_lab_distance
      mai[3] <- mai[3] + max_height_n
      mai[4] <- mai[4] + r_t_width + l_lab_distance
    } else {
      mai[1] <- mai[1] + r_t_width + l_lab_distance
      mai[2] <- mai[2] + max_height_1
      mai[3] <- mai[3] + c_t_width + u_lab_distance
      mai[4] <- mai[4] + max_height_n
    }

    # Set the new margins in inches
    par(mai = mai)

    # Call plot() to create the empty plot object with the correct size
    plot(0, type = "n",
         ylim = y_lim,
         xlim = x_lim,
         axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }

  # Here the space of the plot along the x- or y-axis is calculated.
  if (horizontal) {
    space_size <- y_lim[2] - y_lim[1]
    space_start <- y_lim[1]
    # Convert the distance between boxes and labels from inches to user unit
    # for latex use in plotting with text()
    x_zero_inch <-  grconvertX(0, from = "inches")
    u_lab_distance <- grconvertX(u_lab_distance, from = "inches") - x_zero_inch
    l_lab_distance <- grconvertX(l_lab_distance, from = "inches") - x_zero_inch
  } else {
    space_size <- x_lim[2] - x_lim[1]
    space_start <- x_lim[1]
    # Convert the distance between boxes and labels from inches to user unit
    # for latex use in plotting with text()
    y_zero_inch <- grconvertY(0, from = "inches")
    u_lab_distance <- grconvertY(u_lab_distance, from = "inches") - y_zero_inch
    l_lab_distance <- grconvertY(l_lab_distance, from = "inches") - y_zero_inch
  }

  # If no independent abundances are given for the higher species
  # calculate the sum of all columns as abundances
  if (is.null(higher_abundances)) {
    higher_abundances <- colSums(web)
  } else {
    higher_abundances <- higher_abundances[colnames(web)]
  }
  if (!is.null(add_higher_abundances)) {
    # Sort the additional higher abundances 
    # according to the column-order in the web
    add_higher_abundances <- add_higher_abundances[colnames(web)]
    higher_abundances <- c(rbind(higher_abundances, add_higher_abundances))
    # Include the custom colors into the color vector
    if (length(higher_add_color) == 1) {
      if (higher_add_color == "same") {
        higher_add_color <- higher_color
      } else {
        higher_add_color <- rep_len(higher_add_color, nc)
      }
    } else if (!is.null(names(higher_add_color))) {
      stopifnot(length(higher_add_color) == nc)
      if (!setequal(names(higher_add_color), c_names)) {
        stop("Names of higher_add_color does not match higher species names.")
      }
      higher_add_color <- higher_add_color[colnames(web)]
    } else if (length(higher_add_color) < nc) {
      higher_add_color <- rep_len(higher_add_color, nc)
    }
    higher_color <- c(rbind(higher_color, higher_add_color))
  }

  # If no independent abundances are given for the lower species
  # calculate the sum of all columns as abundances
  if (is.null(lower_abundances)) {
    lower_abundances <- rowSums(web)
  } else {
    lower_abundances <- lower_abundances[rownames(web)]
  }
  if (!is.null(add_lower_abundances)) {
    # Sort the additional lower abundances according to the row-order in the web
    add_lower_abundances <- add_lower_abundances[rownames(web)]
    # Merge the abundances and additional abundances vector
    # by alternating indices
    lower_abundances <- c(rbind(lower_abundances, add_lower_abundances))
    # Include the custom colors into the color vector
    if (length(lower_add_color) == 1) {
      if (lower_add_color == "same") {
        lower_add_color <- lower_color
      } else {
        lower_add_color <- rep_len(lower_add_color, nr)
      }
    } else if (!is.null(names(lower_add_color))) {
      stopifnot(length(lower_add_color) == nr)
      if (!setequal(names(lower_add_color), r_names)) {
        stop("Names of lower_add_color does not match lower species names.")
      }
      lower_add_color <- lower_add_color[rownames(web)]
    } else if (length(lower_add_color) < nr) {
      lower_add_color <- rep_len(lower_add_color, nr)
    }
    lower_color <- c(rbind(lower_color, lower_add_color))
  }

  # if (horizontal) {
  #   theta <- srt * pi / 180
  # } else {
  #   theta <- (srt + 90) * pi / 180
  # }

  # Get the not rotated height and width of the higher labels
  H <- strheight(higher_labels, cex = text_size)
  W <- strwidth(higher_labels, cex = text_size)

  # Calculate the effective rotated height and width
  # str_h_c <- abs(H * cos(theta_pi)) + abs(W * sin(theta_pi))
  str_h_c <- abs(H * sin(theta_pi)) + abs(W * cos(theta_pi))

  # Get the not rotated height and width of the lower labels
  H <- strheight(lower_labels, cex = text_size)
  W <- strwidth(lower_labels, cex = text_size)

  # Calculate the effective rotated height and width
  # str_h_r <- abs(H * cos(theta_pi)) + abs(W * sin(theta_pi))
  str_h_r <- abs(H * sin(theta_pi)) + abs(W * cos (theta_pi))


  # TODO: Document scaling
  if (scaling == "relative") {
    u_scaling_factor <- sum(higher_abundances)
    l_scaling_factor <- sum(lower_abundances)
  } else if (scaling == "absolute") {
    max_abundances <- max(sum(higher_abundances), sum(lower_abundances))
    u_scaling_factor <- max_abundances
    l_scaling_factor <- max_abundances
  }

  # TODO: Document
  if (length(spacing) == 2) {
    c_space <- space_size * spacing[1] / (nc - 1)
    r_space <- space_size * spacing[2] / (nr - 1)
  } else if (length(spacing) == 1) {
    if (is.numeric(spacing)) {
      spacing <- c(spacing, spacing)
      c_space <- space_size * spacing[1] / (nc - 1)
      r_space <- space_size * spacing[2] / (nr - 1)
    } else if (spacing == "auto") {
      # Calculate the maximum overlap of a label with its corresponding box
      # and set the space to that value times 1.05
      space_c <- nc * max(str_h_c - (space_size * higher_abundances / u_scaling_factor))
      space_r <- nr * max(str_h_r - (space_size * lower_abundances / l_scaling_factor))
      space <- 1.1 * max(space_c, space_r, 0.05)
      if (space > space_size) {
        # text_size <- text_size / space * space_size * 0.9
        # space <- 0.9
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

  # Scale the abundances according to the scaling factor
  # (see relative or absolute)
  c_prop_sizes <- (space_size - spacing[1]) * higher_abundances / u_scaling_factor
  r_prop_sizes <- (space_size - spacing[2]) * lower_abundances / l_scaling_factor

  # TODO: Document scaling
  if (scaling == "relative") {
    # c_prop_sizes <- (space_size - spacing[1]) * higher_abundances / sum(higher_abundances)
    # r_prop_sizes <- (space_size - spacing[2]) * lower_abundances / sum(lower_abundances)
  } else if (scaling == "absolute") {
    # max_abundances <- max(sum(higher_abundances), sum(lower_abundances))
    # c_prop_sizes <- (space_size - spacing[1]) * higher_abundances / max_abundances
    # r_prop_sizes <- (space_size - spacing[2]) * lower_abundances / max_abundances
    c_space <- (space_size - sum(c_prop_sizes)) / (nc - 1)
    r_space <- (space_size - sum(r_prop_sizes)) / (nr - 1)
  }

  # TODO: Document additional abundances
  if (!is.null(add_higher_abundances)) {
    nc <- nc * 2
    c_space <- rep(c(0, c_space), length.out = (nc - 1))
  }
  if (!is.null(add_lower_abundances)) {
    nr <- nr * 2
    r_space <- rep(c(0, r_space), length.out = (nr - 1))
  }

  # Calculate the edge positions of the higher boxes
  c_xl <- c(space_start, space_start + cumsum(c_prop_sizes[-nc] + c_space))
  c_xr <- c(space_start + c_prop_sizes[1],
            space_start + c_prop_sizes[1] + cumsum(c_prop_sizes[-1] + c_space))

  # Calculate the edge positions of the lower boxes
  r_xl <- c(space_start, space_start + cumsum(r_prop_sizes[-nr] + r_space))
  r_xr <- c(space_start + r_prop_sizes[1],
            space_start + r_prop_sizes[1] + cumsum(r_prop_sizes[-1] + r_space))

  # If additional abundances are given 
  # the label has the be aligned in the middle between two connected boxes
  if (!is.null(add_higher_abundances)) {
    c_tx <- (c_xl[seq(1, nc, 2)] + c_xr[seq(2, nc, 2)]) / 2
  } else {
    c_tx <- (c_xl + c_xr) / 2
  }
  if (!is.null(add_lower_abundances)) {
    r_tx <- (r_xl[seq(1, nr, 2)] + r_xr[seq(2, nr, 2)]) / 2
  } else {
    r_tx <- (r_xl + r_xr) / 2
  }

  if (higher_border == "same") {
    higher_border <- higher_color
  }
  if (lower_border == "same") {
    lower_border <- lower_color
  }

  c_m_t_width <- max(strwidth(higher_labels, srt = srt, cex = text_size))
  r_m_t_width <- max(strwidth(higher_labels, srt = srt, cex = text_size))

  # Apply italics to all higher labels
  if (higher_italic) {
    higher_labels <- lapply(higher_labels, function(x) bquote(italic(.(x))))
    higher_labels <- as.expression(higher_labels)
  }
  # Apply italics to all lower labels
  if (lower_italic) {
    lower_labels <- lapply(lower_labels, function(x) bquote(italic(.(x))))
    lower_labels <- as.expression(lower_labels)
  }

  if (length(box_size) == 1) {
    higher_box_size <- box_size
    lower_box_size <- box_size
  } else if (length(box_size) == 2) {
    higher_box_size <- box_size[1]
    lower_box_size <- box_size[2]
  }

  # TODO: Document this whole text label alignment mess!
  # Draw the boxes and labels either horizontal or vertical to each other.
  if (horizontal) {
    x_start <- x_lim[1]
    x_end <- x_lim[2]
    rect(x_start, c_xl, x_start + higher_box_size, c_xr,
         col = higher_color, border = higher_border)
    rect(x_end - lower_box_size, r_xl, x_end, r_xr,
         col = lower_color, border = lower_border)
    # If the text is rotated more than 45 degrees
    # center align the labels 
    if (srt == 90 | srt == 270) {
      u_adj <- c(0.5, 0.5)
      l_adj <- c(0.5, 0.5)
    } else if (srt > 90 && srt < 270) {
      u_adj <- c(0, 0.5)
      l_adj <- c(1, 0.5)
    } else {
      u_adj <- c(1, 0.5)
      l_adj <- c(0, 0.5)
    }
    text(x_start - u_lab_distance, c_tx, higher_labels, adj = u_adj, srt = srt,
         cex = text_size, xpd = TRUE, font = font, family = family,
         col = higher_text_color)
    text(x_end + l_lab_distance, r_tx, lower_labels, adj = l_adj, srt = srt,
         cex = text_size, xpd = TRUE, font = font, family = family,
         col = lower_text_color)

  } else {
    y_start <- y_lim[1]
    y_end <- y_lim[2]
    rect(r_xl, y_start, r_xr, y_start + lower_box_size, col = lower_color, border = lower_border)
    rect(c_xl, y_end - higher_box_size, c_xr, y_end, col = higher_color, border = higher_border)

    if (srt == 0 | srt == 180) {
      text(r_tx, y_start - l_lab_distance, lower_labels, pos = 1, cex = text_size,
           xpd = TRUE, srt = srt, font = font, family = family, col = lower_text_color)
      text(c_tx, y_end + u_lab_distance, higher_labels,
           pos = 3, cex = text_size, xpd = TRUE, srt = srt,
           font = font, family = family, col = higher_text_color)
    } else {
      if (srt > 180) {
        adj_l <- c(0, 0.5)
        adj_h <- c(1, 0.5)
      } else {
        adj_l <- c(1, 0.5)
        adj_h <- c(0, 0.5)
      }
      text(r_tx, y_start - l_lab_distance, lower_labels,
           adj = adj_l, cex = text_size, xpd = TRUE, srt = srt,
           font = font, family = family)
      text(c_tx, y_end + u_lab_distance, higher_labels,
           adj = adj_h, cex = text_size, xpd = TRUE, srt = srt,
           font = font, family = family)
    }
  }

  # Interactions
  # TODO: short explanation what is going on
  if (!is.null(add_lower_abundances) && !is.null(add_higher_abundances)) {
    web.df <- data.frame(row = rep(seq(1, nr, 2), nc/2),
                         col = rep(seq(1, nc, 2), each = nr/2),
                         weight = c(web))
  } else if (!is.null(add_lower_abundances)) {
    web.df <- data.frame(row = rep(seq(1, nr, 2), nc),
                         col = rep(1:nc, each = nr/2),
                         weight = c(web))
  } else if (!is.null(add_higher_abundances)) {
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
      x1 <- x_start + higher_box_size
      x2 <- x_end - lower_box_size
    } else {
      x1 <- y_end - lower_box_size#0.9
      x2 <- y_start + higher_box_size#0.1
    }
    y1 <- link$xcoord.tl
    y2 <- link$xcoord.bl
    y3 <- link$xcoord.tr
    y4 <- link$xcoord.br
    if (link_color == "lower") {
      l_col <- lower_color[link$row]
    } else if (link_color == "higher") {
      l_col <- higher_color[link$col]
    } else {
      l_col <- link_color
    }
    draw_link(x1, x2, y1, y2, y3, y4, l_col,
              horizontal = horizontal,
              alpha = link_alpha,
              curved = curved_links,
              border = link_border)
  }

  ## TODO: Implement legend plotting
}


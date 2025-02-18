plot2webs <- function(web1,
                         web2,
                         middle = "higher",
                         sort_by = 1,
                         sorting = "normal",
                         higher_abundances = NULL,
                         lower_abundances = NULL,
                         add_higher_abundances = NULL,
                         add_lower_abundances = NULL,
                         scaling = "absolute",
                         font = NULL,
                         family = NULL,
                         srt = 0,
                         higher_italic = FALSE,
                         lower_italic = FALSE,
                         text_size = 1,
                         x_lim = c(0, 1),
                         y_lim = c(0, 1),
                         lower_color = "black",
                         lower_border = "same",
                         add_lower_color = "same",
                         higher_color = "black",
                         higher_border = "same",
                         add_higher_color = "same",
                         higher_text_color = "black",
                         lower_text_color = "black",
                         horizontal = FALSE,
                         abbr_names = FALSE,
                         link_color = "lower",
                         link_border = "same",
                         link_alpha = 0.6,
                         curved_links = FALSE,
                         arrow = "no",
                         spacing = "auto",
                         space_lower = 0.2,
                         space_higher = 0.2,
                         mar = c(1, 1, 1, 1),
                         mai = NULL,
                         plot_axes = FALSE) {


  # By default the the higher species (columns)
  # are plotted in the middle. So to enable plotting the lower
  # species (rows) in the middle a trick is used. That is,
  # both webs are transposed and all corresponding variables are swapped.
  if (middle == "lower") {
    web1 <- t(web1)
    web2 <- t(web2)
    # swap independent abundances
    tmp <- higher_abundances 
    higher_abundances <- lower_abundances
    lower_abundances <- tmp
    # swap additional abundances
    tmp <- add_higher_abundances
    add_higher_abundances <- add_lower_abundances
    add_lower_abundances <- tmp
    # swap italics
    tmp <- higher_italic
    higher_italic <- lower_italic
    lower_italic <- tmp
    # swap color
    tmp <- higher_color
    higher_color <- lower_color
    lower_color <- tmp
    # swap border color
    tmp <- higher_border
    higher_border <- lower_border
    lower_border <- tmp
    # swap additional color
    tmp <- add_higher_color
    add_higher_color <- add_lower_color
    add_lower_color <- tmp
    # swap text color
    tmp <- higher_text_color
    higher_text_color <- lower_text_color
    lower_text_color <- tmp

    if (link_color == "higher") {
      link_color <- "lower"
    } else if (link_color == "lower") {
      link_color <- "higher"
    }

  }

  if (sort_by == 1) {
    web1 <- sortweb2(web1, sort.order = sorting)
  } else if (sort_by == 2) {
    web2 <- sortweb2(web2, sort.order = sorting)
  }

  r_names_1 <- rownames(web1)
  r_names_2 <- rownames(web2)
  c_names_1 <- colnames(web1)
  c_names_2 <- colnames(web2)

  if (sort_by == 1) {
    c_names <- union(c_names_1, c_names_2)
    r_names <- union(r_names_1, r_names_2)
  } else if (sort_by == 2) {
    c_names <- union(c_names_2, c_names_1)
    r_names <- union(r_names_2, r_names_1)
  }

  nr_1 <- nrow(web1)
  nr_2 <- nrow(web2)

  # Calculate the set difference between the column names
  c_diff_1 <- setdiff(c_names_2, c_names_1)
  c_fill_1 <- matrix(0, nrow = nr_1,
                     ncol = length(c_diff_1), dimnames = list(c(), c_diff_1))
  web1 <- cbind(web1, c_fill_1)
  c_diff_2 <- setdiff(c_names_1, c_names_2)
  c_fill_2 <- matrix(0, nrow = nr_2,
                     ncol = length(c_diff_2), dimnames = list(c(), c_diff_2))
  web2 <- cbind(web2, c_fill_2)

  if (sort_by == 1) {
    #print("SORT")
    web2 <- web2[intersect(r_names, r_names_2), colnames(web1)]
    r_names_2 <- rownames(web2)
  } else if (sort_by == 2) {
    web1 <- web1[intersect(r_names, r_names_1), colnames(web2)]
    r_names_1 <- rownames(web1)
  }

  nc <- length(c_names)
  nr <- length(union(r_names_1, r_names_2))
  nc_1 <- ncol(web1)
  nc_2 <- ncol(web2)

  # Recycle the color vectors
  if (!is.null(names(higher_color))) {
    stopifnot(length(higher_color) == nc)
    if (!setequal(names(higher_color), c_names)) {
      # TODO: fix bug when middle = "lower"
      stop("Names of higher_color does not match names of higher species.")
    }
    #higher_color <- higher_color[colnames(web1)]
    higher_color <- higher_color[c_names]
  } else {
    if (length(higher_color) < nc) {
      higher_color <- rep_len(higher_color, nc)
    }
  }

  if (!is.null(names(lower_color))) {
    stopifnot(length(lower_color) == nr)
    if (!setequal(names(lower_color), r_names)) {
      # TODO: fix bug when middle = "lower"
      stop("Names of lower_color does not match names of lower species.")
    }
    #higher_color <- higher_color[colnames(web1)]
    lower_color_1 <- lower_color[r_names_1]
    lower_color_2 <- lower_color[r_names_2]
  } else if (length(lower_color) < nr) {
    lower_color_1 <- rep_len(lower_color, nr_1)
    lower_color_2 <- rep_len(lower_color, nr_2)
  }

  #c_sums_1 <- rep(0, nc)
  #c_sums_2 <- rep(0, nc)

  # for (i in seq_along(c_names)) {
  #   c_name <- c_names[i]
  #   if (c_name %in% c_names_1) {
  #     c_sums_1[i] <- sum(web1[, c_name])
  #   }
  #   if (c_name %in% c_names_2) {
  #     c_sums_2[i] <- sum(web2[, c_name])
  #   }
  # }

  c_abuns_1 <- colSums(web1)
  c_abuns_2 <- colSums(web2)

  if (!is.null(higher_abundances)) {
    c_abuns_1 <- higher_abundances[[1]]
    c_abuns_2 <- higher_abundances[[2]]
  }
  if (!is.null(add_higher_abundances)) {
    # Assert the abundances are given as a list of length 2
    stopifnot(is.list(add_higher_abundances),
              length(add_higher_abundances) == 2)
    # Extract the abundances for each web
    add_c_abundances_1 <- add_higher_abundances[[1]]
    add_c_abundances_2 <- add_higher_abundances[[2]]
    add_c_abun_names_1 <- names(add_c_abundances_1)
    add_c_abun_names_2 <- names(add_c_abundances_2)

    # Check if vector 1 contains all species
    if (!is.null(add_c_abun_names_1)) {
      if (length(add_c_abun_names_1) == nc) {
        if (setequal(add_c_abun_names_1, c_names)) {
          # Sort the additional higher abundances
          # according to the column-order in the web
          add_c_abundances_1 <- add_c_abundances_1[colnames(web1)]
        }
      } else { # The length does not fit
        warning(paste("The length of the higher additional abundances of web1",
                      "does not match with the number of higher species."))
        stop()
      }
    } else { # unnamed list
      warning(
        paste("The higher additional abundances of web1 are an unnamed list.",
              "Assuming the they are given in the same order as the web.",
              "This easily leads to unwanted behavior."))
      if (length(add_c_abundances_1) == nc) {
        # TODO sort the vector correctly
        add_c_abundances_1 <- add_c_abundances_1
      } else {
        stop(paste("The length of the higher additional abundances of web1",
                   "does not match with the number of higher species."))
      }
    }

    # Check if vector 2 contains all species
    if (!is.null(add_c_abun_names_2)) {
      if (length(add_c_abun_names_2) == nc) {
        if (setequal(add_c_abun_names_2, c_names)) {
          # Sort the additional higher abundances
          # according to the column-order in the web
          add_c_abundances_2 <- add_c_abundances_2[colnames(web1)]
        }
      } else {
        warning(paste("The length of the higher additional abundances of web2",
                      "does not match with the number of higher species."))
        stop()
      }
    } else { # unnamed list
      warning(
        paste("The higher additional abundances of web2 are an unnamed list.",
              "Assuming the they are given in the same order as the web.",
              "This easily leads to unwanted behavior."))
      if (length(add_c_abundances_2) == nc) {
        # TODO sort the vector correctly
        add_c_abundances_2 <- add_c_abundances_2
      } else {
        stop(paste("The length of the higher additional abundances of web2",
                   "does not match with the number of higher species."))
      }
    }
    add_c_abundances_1[is.na(add_c_abundances_1)] <- 0
    add_c_abundances_2[is.na(add_c_abundances_2)] <- 0
    c_abuns_1 <- c(rbind(c_abuns_1, add_c_abundances_1))
    c_abuns_2 <- c(rbind(c_abuns_2, add_c_abundances_2))

    if (add_higher_color == "same") {
      higher_color <- rep(higher_color, each = 2)
    } else if (!is.null(names(add_higher_color))) {
      stopifnot(length(add_higher_color) == nc)
      if (!setequal(names(add_higher_color), c_names)) {
        stop("Names of add_higher_color does not match higher species names.")
      }
      add_higher_color <- add_higher_color[colnames(web1)]
      higher_color <- c(rbind(higher_color, add_higher_color))
    } else if (length(add_higher_color) < nc) {
      add_higher_color <- rep_len(higher_color, nc)
      higher_color <- c(rbind(higher_color, add_higher_color))
    }
  }

  r_abuns_1 <- rowSums(web1)
  r_abuns_2 <- rowSums(web2)
  if (!is.null(lower_abundances)) {
    r_abuns_1 <- lower_abundances[[1]]
    r_abuns_2 <- lower_abundances[[2]]
  }
  if (!is.null(add_lower_abundances)) {
    stopifnot(is.list(add_lower_abundances),
              length(add_lower_abundances) == 2)
    # Sort the additional lower abundances according to the row-order in the web
    add_r_abundances_1 <- add_lower_abundances[[1]]#[rownames(web1)]
    add_r_abundances_2 <- add_lower_abundances[[2]]#[rownames(web2)]
    add_r_abundances_1[is.na(add_r_abundances_1)] <- 0
    add_r_abundances_2[is.na(add_r_abundances_2)] <- 0
    # Merge the abundances and additional abundances vector
    # by alternating indices
    #lower_abundances <- c(rbind(lower_abundances, add_lower_abundances))
    r_abuns_1 <- c(rbind(r_abuns_1, add_r_abundances_1))
    r_abuns_2 <- c(rbind(r_abuns_2, add_r_abundances_2))
    # TODO: fix bug with new lower_color_1 and lower_color_2
    if (add_lower_color == "same") {
      lower_color <- rep(lower_color, each = 2)
    } else if (!is.null(names(add_lower_color))) {
      stopifnot(length(add_lower_color) == nr)
      lower_color <- c(rbind(lower_color, add_lower_color))
    } else {
      lower_color <- c(rbind(lower_color, add_lower_color))
    }
  }

  if (scaling == "relative") {
    c_prop_sizes_1 <- c_abuns_1 / sum(c_abuns_1)
    c_prop_sizes_2 <- c_abuns_2 / sum(c_abuns_2)
  } else if (scaling == "absolute") {
    total_abuns <- sum(c(c_abuns_1, c_abuns_2, r_abuns_1, r_abuns_2))
    c_prop_sizes_1 <- c_abuns_1 / total_abuns
    c_prop_sizes_2 <- c_abuns_2 / total_abuns
  }

  # c_prop_sizes_1 <- colSums(web1) / sum(colSums(web1))
  # c_prop_sizes_2 <- colSums(web2) / sum(colSums(web2))

  c_prop_sizes_max <- pmax(c_prop_sizes_1, c_prop_sizes_2)

  c_prop_sizes_1 <- (1 - space_higher) * c_prop_sizes_1
  c_prop_sizes_2 <- (1 - space_higher) * c_prop_sizes_2

  # Adjust the lower space so that it matches the necessary higher space
  if (scaling == "relative") {
    space_lower <- 1 - ((1 - space_higher) / sum(c_prop_sizes_max))
    r_prop_sizes_1 <- (1 - space_lower) * r_abuns_1 / sum(r_abuns_1)
    r_prop_sizes_2 <- (1 - space_lower) * r_abuns_2 / sum(r_abuns_2)
    r_space_1 <- space_lower / (nr_1 - 1)
    r_space_2 <- space_lower / (nr_2 - 1)
  } else if (scaling == "absolute") {
    r_prop_sizes_1 <- r_abuns_1 / total_abuns# / sum(c_prop_sizes_max)
    r_prop_sizes_2 <- r_abuns_2 / total_abuns# / sum(c_prop_sizes_max)
  }

  max_prop_size <- max(sum(c_prop_sizes_max), sum(r_prop_sizes_1), sum(r_prop_sizes_2))
  c_prop_sizes <- (1 - space_higher) * c_prop_sizes_max / max_prop_size

  c_space <- space_higher / (nc - 1)

  if (scaling == "absolute") {
    r_prop_sizes_1 <- (1 - space_lower) * r_prop_sizes_1 / max_prop_size
    r_prop_sizes_2 <- (1 - space_lower) * r_prop_sizes_2 / max_prop_size
    r_space_1 <- (1 - sum(r_prop_sizes_1)) / (nr_1 - 1)
    r_space_2 <- (1 - sum(r_prop_sizes_2)) / (nr_2 - 1)
    # c_prop_sizes_1 <- (1 - space_higher) * c_prop_sizes_1 / max_prop_size
    # c_prop_sizes_2 <- (1 - space_higher) * c_prop_sizes_2 / max_prop_size
    c_prop_sizes_1 <- c_prop_sizes_1 / max_prop_size
    c_prop_sizes_2 <- c_prop_sizes_2 / max_prop_size
    c_space <- (1 - sum(c_prop_sizes)) / (nc - 1)
  } else if (scaling == "relative") {
    c_prop_sizes_1 <- c_prop_sizes_1 / sum(c_prop_sizes_max)
    c_prop_sizes_2 <- c_prop_sizes_2 / sum(c_prop_sizes_max)
    #print(sum(r_prop_sizes_1))
    #print(sum(r_prop_sizes_2))
    #print(sum(c_prop_sizes_1))
    #print(sum(c_prop_sizes_2))
  }

  # if (scaling == "relative") {
  # } else if (scaling == "absolute") {
  #   r_prop_sizes_1 <- (1 - space_lower) * rowSums(web1) / sum(c(rowSums(web1), rowSums(web2)))
  #   r_prop_sizes_2 <- (1 - space_lower) * rowSums(web2) / sum(c(rowSums(web1), rowSums(web2)))
  # }

  if (!is.null(add_higher_abundances)) {
    c_space <- c(0, c_space)
    nc <- nc * 2
  }
  if (!is.null(add_lower_abundances)) {
    r_space_1 <- c(0, r_space_1)
    r_space_2 <- c(0, r_space_2)
    nr_1 <- nr_1 * 2
    nr_2 <- nr_2 * 2
  }

  # r_space_1 <- space_lower / (nr_1 - 1)
  # r_space_2 <- space_lower / (nr_2 - 1)

  c_xl <- c(0, cumsum(c_prop_sizes[-nc] + c_space))
  c_xr <- c(c_prop_sizes[1],
            c_prop_sizes[1] + cumsum(c_prop_sizes[-1] + c_space))

  r_xl_1 <- c(0, cumsum(r_prop_sizes_1[-nr_1] + r_space_1))
  r_xr_1 <- c(r_prop_sizes_1[1],
              r_prop_sizes_1[1] + cumsum(r_prop_sizes_1[-1] + r_space_1))

  r_xl_2 <- c(0, cumsum(r_prop_sizes_2[-nr_2] + r_space_2))
  r_xr_2 <- c(r_prop_sizes_2[1],
              r_prop_sizes_2[1] + cumsum(r_prop_sizes_2[-1] + r_space_2))


  # c_m_t_width <- 2 * max(strwidth(c_names))
  # r_m_t_width_1 <- 2 * max(strwidth(r_names_1))
  # r_m_t_width_2 <- 2 * max(strwidth(r_names_2))

  c_m_t_width <- max(strwidth(c_names, units = "inches", cex = text_size))
  r_m_t_width_1 <- max(strwidth(r_names_1, units = "inches", cex = text_size))
  r_m_t_width_2 <- max(strwidth(r_names_2, units = "inches", cex = text_size))
  #print(c_m_t_width)
  #print(r_m_t_width_1)
  #print(r_m_t_width_2)

  if (horizontal) {
    par(fig=c(0,0.5,0,1), mai = c(0.5, r_m_t_width_1 + 0.1, 0.5, c_m_t_width/2 + 0.1))
  } else {
    par(fig=c(0,1,0,0.5), mai = c(r_m_t_width_1 + 0.1, 0.5, c_m_t_width/2 + 0.1, 0.5))
  }
  plot(0, type = "n", ylim = c(0, 1), xlim = c(0, 1),
       axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")


  if (!is.null(add_higher_abundances)) {
    c_xl_1 <- rep(0, length(c_prop_sizes))
    c_xr_1 <- rep(0, length(c_prop_sizes))
    c_xl_1[c(TRUE, FALSE)] <- c_xl[c(TRUE, FALSE)] + 0.5 * (c_prop_sizes[c(TRUE, FALSE)] - c_prop_sizes_1[c(TRUE, FALSE)])
    c_xr_1[c(TRUE, FALSE)] <- c_xr[c(TRUE, FALSE)] - 0.5 * (c_prop_sizes[c(TRUE, FALSE)] - c_prop_sizes_1[c(TRUE, FALSE)])
    c_xl_1[c(FALSE, TRUE)] <- c_xr_1[c(TRUE, FALSE)]
    c_xr_1[c(FALSE, TRUE)] <- c_xr_1[c(TRUE, FALSE)] + c_prop_sizes_1[c(FALSE, TRUE)]

    c_xl_2 <- rep(0, length(c_prop_sizes))
    c_xr_2 <- rep(0, length(c_prop_sizes))
    c_xl_2[c(TRUE, FALSE)] <- c_xl[c(TRUE, FALSE)] + 0.5 * (c_prop_sizes[c(TRUE, FALSE)] - c_prop_sizes_2[c(TRUE, FALSE)])
    c_xr_2[c(TRUE, FALSE)] <- c_xr[c(TRUE, FALSE)] - 0.5 * (c_prop_sizes[c(TRUE, FALSE)] - c_prop_sizes_2[c(TRUE, FALSE)])
    c_xl_2[c(FALSE, TRUE)] <- c_xr_2[c(TRUE, FALSE)]
    c_xr_2[c(FALSE, TRUE)] <- c_xr_2[c(TRUE, FALSE)] + c_prop_sizes_2[c(FALSE, TRUE)]
  } else {
    c_xl_1 <- c_xl + 0.5 * (c_prop_sizes - c_prop_sizes_1)
    c_xr_1 <- c_xr - 0.5 * (c_prop_sizes - c_prop_sizes_1)
    c_xl_2 <- c_xl + 0.5 * (c_prop_sizes - c_prop_sizes_2)
    c_xr_2 <- c_xr - 0.5 * (c_prop_sizes - c_prop_sizes_2)
  }
  # rect(0.2, c_xl_1, 0.3, c_xr_1, col = higher_color, border = higher_color)

  # rect(0.4, c_xl_2, 0.5, c_xr_2, col = higher_color, border = higher_color)

  if (higher_border == "same") {
    higher_border <- higher_color
  }
  if (lower_border == "same") {
    lower_border_1 <- lower_color_1
    lower_border_2 <- lower_color_2
  }


  if (!is.null(add_higher_abundances)) {
    c_tx <- (c_xl[seq(1, nc, 2)] + c_xr[seq(2, nc, 2)]) / 2
  } else{
    c_tx <- (c_xl + c_xr) / 2
  }
  if (!is.null(add_lower_abundances)) {
    r_tx_1 <- (r_xl_1[seq(1, nr_1, 2)] + r_xr_1[seq(2, nr_1, 2)]) / 2
  } else {
    r_tx_1 <- (r_xl_1 + r_xr_1) / 2
  }

  # r_tx_1 <- (r_xr_1 + r_xl_1) / 2
  # c_tx <- (c_xr + c_xl) / 2

  if (higher_italic) {
    c_names <- lapply(c_names, function(x) bquote(italic(.(x))))
    c_names <- as.expression(c_names)
  }
  if (lower_italic) {
    r_names_1 <- lapply(r_names_1, function(x) bquote(italic(.(x))))
    r_names_1 <- as.expression(r_names_1)
  }

  if (horizontal) {
    rect(0, r_xl_1, 0.1, r_xr_1, col = lower_color_1, border = lower_border_1)
    rect(0.9, c_xl_1, 1, c_xr_1, col = higher_color, border = higher_border)
    text(-0.01, r_tx_1, r_names_1, adj = c(1, 0.5), xpd=T,
         cex = text_size, srt = srt, col = lower_text_color)
    text(1 + 0.5 * grconvertX(c_m_t_width + r_m_t_width_1 + 0.3, from="inches"), 
         c_tx, c_names, adj = c(0.5, 0.5), xpd=NA, cex = text_size,
         srt = srt, col = higher_text_color)
  } else {
    rect(r_xl_1, 0, r_xr_1, 0.1, col = lower_color_1, border = lower_border_1)
    rect(c_xl_1, 0.9, c_xr_1, 1, col = higher_color, border = higher_border)
    text(r_tx_1, -0.01, r_names_1, adj = c(1, 0.5), xpd=T,
         cex = text_size, srt = srt + 90, col = lower_text_color)
    text(c_tx, 1 + 0.5 * grconvertY(c_m_t_width + r_m_t_width_1 + 0.3, from="inches"),
         c_names, adj = c(0.5, 0.5), xpd=NA, cex = text_size,
         srt = srt + 90, col = higher_text_color)
  }
  # text(1, c_tx, c_names, adj = c(0, 0.5), xpd=NA)

  # # Interactions
  # web.df_1 <- data.frame(row = rep(1:nr_1, nc), col = rep(1:nc, each = nr_1), weight = c(web1))
  # web.df_1 <- data.frame(row = rep(1:nr_1, nc_1), col = rep(which(c_names %in% c_names_1), each = nr_1), weight = c(web1))
  web.df_1 <- data.frame(row = rep(1:nr_1, nc), col = rep(1:nc, each = nr_1), weight = c(web1))
  if (!is.null(add_higher_abundances) && !is.null(add_lower_abundances)) {
    web.df_1 <- data.frame(row = rep(seq(1, nr_1, 2), nc / 2),
                           col = rep(seq(1, nc, 2), each = nr_1 / 2),
                           weight = c(web1))
  } else if (!is.null(add_lower_abundances)) {
    web.df_1 <- data.frame(row = rep(seq(1, nr_1, 2), nc),
                           col = rep(1:nc, each = nr_1 / 2),
                           weight = c(web1))
  } else if (!is.null(add_higher_abundances)) {
    web.df_1 <- data.frame(row = rep(1:nr_1, nc / 2),
                           col = rep(seq(1, nc, 2), each = nr_1),
                           weight = c(web1))
  }
  web.df_1 <- web.df_1[web.df_1$weight > 0, ]
  # x-coordinates of interactions: tl=topleft, etc
  web.df_1[, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")] <- NA 

  # # low coordinates for interactions (in order of the web.df)
  for (i in unique(web.df_1$row)) { # for i in lower species
    # i <- 3
    links.i <- web.df_1[web.df_1$row == i, ]
    relpos <- cumsum(links.i$weight) / sum(links.i$weight)
    coords.int.low <- (r_xl_1[i] + relpos * (r_xr_1[i] - r_xl_1[i]))
    web.df_1[web.df_1$row == i, "xcoord.bl"] <- c(r_xl_1[i], coords.int.low[-nrow(links.i)])
    web.df_1[web.df_1$row == i, "xcoord.br"] <- c(coords.int.low)
  }

  # # high coordinates for interactions (in order of the web.df)
  for (j in unique(web.df_1$col)) { # for j in higher species
    # j <- 3
    links.j <- web.df_1[web.df_1$col == j, ]
    relpos <- cumsum(links.j$weight) / sum(links.j$weight)
    coords.int.high <- (c_xl_1[j] + relpos * (c_xr_1[j] - c_xl_1[j]))
    web.df_1[web.df_1$col == j, "xcoord.tl"] <- c(c_xl_1[j], coords.int.high[-nrow(links.j)])
    web.df_1[web.df_1$col == j, "xcoord.tr"] <- c(coords.int.high)
  }
  for (linki in order(-web.df_1$weight)) {
    link <- web.df_1[linki, ]
    x1 <- 0.1
    x2 <- 0.9
    y1 <- link$xcoord.tr
    y2 <- link$xcoord.br
    y3 <- link$xcoord.tl
    y4 <- link$xcoord.bl
    if (link_color == "lower") {
      l_col <- lower_color_1[link$row]
    } else if (link_color == "higher") {
      l_col <- higher_color[link$col]
    } else {
      l_col <- link_color
    }
    draw_link(x2, x1, y1, y2, y3, y4,
              l_col, horizontal = horizontal,
              curved = curved_links)
  }

  # par(fig=c(0.5, 0.5, 0, 1), mai = c(0, c_m_t_width/2, 0, c_m_t_width/2), new = TRUE)
  # plot(0, type = "n", ylim = c(0, 1), xlim = c(0, 1),
  #      axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  if(horizontal) {
    par(fig=c(0.5, 1, 0, 1), mai = c(0.5, c_m_t_width / 2 + 0.1, 0.5, r_m_t_width_2 + 0.1), new = TRUE)
  } else {
    par(fig=c(0, 1, 0.5, 1), mai = c(c_m_t_width / 2 + 0.1, 0.5, r_m_t_width_2 + 0.1, 0.5), new = TRUE)
  }
  plot(0, type = "n", ylim = c(0, 1), xlim = c(0, 1),
       axes = plot_axes, xlab = "", ylab = "", xaxs = "i", yaxs = "i")


  # r_tx_2 <- (r_xr_2 + r_xl_2) / 2
  if (!is.null(add_lower_abundances)) {
    r_tx_2 <- (r_xl_2[seq(1, nr_2, 2)] + r_xr_2[seq(2, nr_2, 2)]) / 2
  } else {
    r_tx_2 <- (r_xl_2 + r_xr_2) / 2
  }


  if (lower_italic) {
    r_names_2 <- lapply(r_names_2, function(x) bquote(italic(.(x))))
    r_names_2 <- as.expression(r_names_2)
  }

  if (horizontal) {
    rect(0, c_xl_2, 0.1, c_xr_2, col = higher_color, border = higher_border)
    rect(0.9, r_xl_2, 1, r_xr_2, col = lower_color_2, border = lower_border_2)
    text(1.01, r_tx_2, r_names_2, adj = c(0, 0.5), cex = text_size,
         xpd=T, srt = srt, col = lower_text_color)
  } else {
    rect(c_xl_2, 0, c_xr_2, 0.1, col = higher_color, border = higher_border)
    rect(r_xl_2, 0.9, r_xr_2, 1, col = lower_color_2, border = lower_border_2)
    text(r_tx_2, 1.01, r_names_2, adj = c(0, 0.5), cex = text_size,
         xpd=T, srt = srt + 90, col = lower_text_color)
  }
  #text(-0.5 * c_m_t_width, c_tx, c_names, adj = c(0.5, 0.5), xpd=T)

  # Interactions
  # web.df_2 <- data.frame(row = rep(1:nr_2, nc_2), col = rep(which(c_names %in% c_names_2), each = nr_2), weight = c(web2))
  if (!is.null(add_higher_abundances) && !is.null(add_lower_abundances)) {
    web.df_2 <- data.frame(row = rep(seq(1, nr_2, 2), nc / 2),
                           col = rep(seq(1, nc, 2), each = nr_2 / 2),
                           weight = c(web2))
  } else if (!is.null(add_lower_abundances)) {
    web.df_2 <- data.frame(row = rep(seq(1, nr_2, 2), nc),
                           col = rep(1:nc, each = nr_2 / 2),
                           weight = c(web2))
  } else if (!is.null(add_higher_abundances)) {
    web.df_2 <- data.frame(row = rep(1:nr_2, nc / 2),
                           col = rep(seq(1, nc, 2), each = nr_2),
                           weight = c(web2))
  } else {
    web.df_2 <- data.frame(row = rep(1:nr_2, nc), col = rep(1:nc, each = nr_2), weight = c(web2))
  }
  web.df_2 <- web.df_2[web.df_2$weight > 0, ]
  web.df_2[, c("xcoord.tl", "xcoord.tr", "xcoord.br", "xcoord.bl")] <- NA # x-coordinates of interactions: tl=topleft, etc

  # low coordinates for interactions (in order of the web.df)
  for (i in unique(web.df_2$row)) { # for i in lower species
    # i <- 3
    links.i <- web.df_2[web.df_2$row == i, ]
    relpos <- cumsum(links.i$weight) / sum(links.i$weight)
    coords.int.low <- (r_xl_2[i] + relpos * (r_xr_2[i] - r_xl_2[i]))
    web.df_2[web.df_2$row == i, "xcoord.bl"] <- c(r_xl_2[i], coords.int.low[-nrow(links.i)])
    web.df_2[web.df_2$row == i, "xcoord.br"] <- c(coords.int.low)
  }

  # high coordinates for interactions (in order of the web.df)
  for (j in unique(web.df_2$col)) { # for j in higher species
    # j <- 3
    links.j <- web.df_2[web.df_2$col == j, ]
    relpos <- cumsum(links.j$weight) / sum(links.j$weight)
    coords.int.high <- (c_xl_2[j] + relpos * (c_xr_2[j] - c_xl_2[j]))
    web.df_2[web.df_2$col == j, "xcoord.tl"] <- c(c_xl_2[j], coords.int.high[-nrow(links.j)])
    web.df_2[web.df_2$col == j, "xcoord.tr"] <- c(coords.int.high)
  }
  for (linki in order(-web.df_2$weight)) {
    link <- web.df_2[linki, ]
    x1 <- 0.1
    x2 <- 0.9
    y1 <- link$xcoord.tr
    y2 <- link$xcoord.br
    y3 <- link$xcoord.tl
    y4 <- link$xcoord.bl
    if (link_color == "lower") {
      l_col <- lower_color_2[link$row]
    } else if (link_color == "higher") {
      l_col <- higher_color[link$col]
    } else {
      l_col <- link_color
    }
    draw_link(x1, x2, y1, y2, y3, y4,
              l_col, horizontal = horizontal,
              curved = curved_links)
  }
}
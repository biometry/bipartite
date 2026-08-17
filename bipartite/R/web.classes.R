###
# @author: Tobias Bauer

# Generics ------------------------------------------------
# Access and replacement functions are S3 generics
higher_attributes <- function(x, web_id = NULL) {
  UseMethod("higher_attributes")
}
`higher_attributes<-` <- function(x, web_id = NULL, value) {
  UseMethod("higher_attributes<-")
}


lower_attributes <- function(x, web_id = NULL) {
  UseMethod("lower_attributes")
}
`lower_attributes<-` <- function(x, web_id = NULL, value) {
  UseMethod("lower_attributes<-")
}

meta_attributes <- function(x, web_id = NULL) {
  UseMethod("meta_attributes")
}
`meta_attributes<-` <- function(x, web_id = NULL, value) {
  UseMethod("meta_attributes<-")
}

# Base class bipartite_web --------------------------------
## Basic constructor --------------------------------------
bipartite_web <- function(web,
                          higher_attributes = NULL,
                          lower_attributes = NULL,
                          meta_attributes = list()) {
  # Validate the underlying object
  if (!is.matrix(web)) {
    web <- as.matrix(web)
  }

  # Default metadata
  if (is.null(lower_attributes)) {
    lower_attributes <- data.frame(
      row.names = rownames(web)
    )
  }

  if (is.null(higher_attributes)) {
    higher_attributes <- data.frame(
      row.names = colnames(web)
    )
  }

  # Validate metadata dimensions
  if (!is.data.frame(lower_attributes) || nrow(lower_attributes) != nrow(web)) {
    stop("`lower_attributes` must be a data frame with one row per matrix row.")
  }

  if (!is.data.frame(higher_attributes) || nrow(higher_attributes) != ncol(web)) {
    stop("`higher_attributes` must be a data frame with one row per matrix column.")
  }

  if (!is.list(meta_attributes)) {
    stop("`meta_attributes` must be a list.")
  }

  structure(
    web,
    lower_attributes = lower_attributes,
    higher_attributes = higher_attributes,
    meta_attributes = meta_attributes,
    class = c("bipartite_web", "matrix", "array")
  )
}

## Subsetting ---------------------------------------------
`[.bipartite_web` <- function(x, i, j, drop = TRUE) {

  # Subset the actual matrix
  y <- NextMethod("[", drop = drop)

  # If dimensions have been dropped, we no longer have
  # a my_matrix in the intended sense.
  if (drop && !is.matrix(y)) {
    return(y)
  }

  lower_attributes <- attr(x, "lower_attributes")
  higher_attributes <- attr(x, "higher_attributes")
  meta_attributes <- attr(x, "meta_attributes")

  if (!missing(i)) {
    lower_attributes <- lower_attributes[i, , drop = FALSE]
  }

  if (!missing(j)) {
    higher_attributes <- higher_attributes[j, , drop = FALSE]
  }

  bipartite_web(
    y,
    lower_attributes = lower_attributes,
    higher_attributes = higher_attributes,
    meta_attributes = meta_attributes
  )
}

## Access functions ---------------------------------------

higher_attributes.bipartite_web <- function(x, web_id) {
  attr(x, "higher_attributes")
}

lower_attributes.bipartite_web <- function(x, web_id) {
  attr(x, "lower_attributes")
}

meta_attributes.bipartite_web <- function(x, web_id) {
  attr(x, "meta_attributes")
}

## Replacement functions ----------------------------------
`higher_attributes<-.bipartite_web` <- function(x, web_id, value) {
  if (!is.data.frame(value)) {
    stop("`higher_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`higher_attributes` must have one row per matrix column")
  }

  attr(x, "higher_attributes") <- value
  x
}

`lower_attributes<-.bipartite_web` <- function(x, web_id, value) {
  if (!is.data.frame(value)) {
    stop("`lower_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`lower_attributes` must have one row per matrix column")
  }

  attr(x, "lower_attributes") <- value
  x
}

`meta_attributes<-.bipartite_web` <- function(x, web_id, value) {
  if (!is.list(value)) {
    stop("`meta_attributes` must be a list.")
  }

  attr(x, "meta_attributes") <- value
  x
}

## Print function -----------------------------------------
print.bipartite_web <- function(x, ...) {
  cat("<bipartite_web object>\n")
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
  cat(strrep("-", 50), "\n")
  cat("Higher species attributes:\n")
  print(higher_attributes(x))
  cat("Lower species attributes:\n")
  print(lower_attributes(x))
  cat("Meta attributes:\n")
  print(meta_attributes(x))
  cat(strrep("-", 50), "\n")
  cat("Matrix:\n")
  print(matrix(x,
               nrow = nrow(x),
               ncol = ncol(x),
               dimnames = list(rownames(x), colnames(x))))
}

## Plot function prototype --------------------------------
## Wrapper for plotweb function.
## TODO: Decide whether to keep / extent the auto color functionality
plot.bipartite_web <- function(x,
                               higher_color = "black",
                               higher_color_attr = NULL,
                               lower_color_attr = NULL,
                               ...) {
  if (!is.null(higher_color_attr)) {
    stopifnot(higher_color_attr %in% colnames(higher_attributes(x)))
    higher_color_attr_vec <- higher_attributes(x)[[higher_color_attr]]
    higher_names <- rownames(higher_attributes(x))
    print(typeof(higher_color_attr_vec))
    if (is.numeric(higher_color_attr_vec)) {
      ramp <- colorRamp(c("white", "black"))
      max_v <- max(higher_color_attr_vec)
      min_v <- min(higher_color_attr_vec)
      normalized_vector <- (higher_color_attr_vec - min_v) / (max_v - min_v)
      print(normalized_vector)
      higher_color <- rgb(ramp(normalized_vector), maxColorValue = 255)
      names(higher_color) <- higher_names
      print(higher_color)
    }
  }
  plotweb(x, higher_color = higher_color, ...)
}


# bipartite_webarray class --------------------------------
## Wide approach ------------------------------------------
## Separate data.frame for each web and each level.
## So w * h + w * l local data.frames in total.
## Additonally one global data.frame for higher and lower each.
## TODO: Add constructor for list of bipartite_web objects
## TODO: Maybe make all constructors generics?
bipartite_webarray <- function(webs, ...,
                               global_higher_attributes = NULL,
                               global_lower_attributes = NULL,
                               local_higher_attributes = NULL,
                               local_lower_attributes = NULL,
                               global_meta_attributes = NULL,
                               local_meta_attributes = NULL) {
  # Combining the matrix inputs into a 3D-array
  if (is.array(webs) && length(dim(webs)) == 3) {
    warray <- webs
  } else {
    # Reused code from web2array by Carsten and Jochen
    if (is(webs, "list")) {
      web_names <- names(webs)
      if (is.null(web_names)) {
        web_names <- paste0("web", seq_along(webs))
      }
    } else {
      # combine all ... and x into a single list:
      webinput <- substitute(list(webs, ...))
      web_names <- as.character(webinput)[-1L]
      webs <- eval(webinput)
    }

    if (length(webs) < 2) stop("Please provide more than one web as input!")
    all.names.HL <- sort(unique(unlist(lapply(webs, colnames)))) # this was sapply before, which caused duplicate names (unique didn't work)
    all.names.LL <- sort(unique(unlist(lapply(webs, rownames)))) # was sapply before
    warray <- array(0,
                    dim = c(length(all.names.LL),
                            length(all.names.HL),
                            length(web_names)),
                    dimnames = list(all.names.LL, all.names.HL, web_names))
    for (i in seq_along(web_names)){
      ri <- rownames(webs[[i]])
      ci <- colnames(webs[[i]])
      warray[ri, ci, i] <- webs[[i]]
    }
  }
  # web_names <- dimnames(warray)[[3]]
  # Default attributes
  if (is.null(global_higher_attributes)) {
    global_higher_attributes <- data.frame(
      row.names = colnames(warray)
    )
  }
  if (is.null(global_lower_attributes)) {
    global_lower_attributes <- data.frame(
      row.names = rownames(warray)
    )
  }
  if (is.null(local_higher_attributes)) {
    local_higher_attributes <- list()
    for (name in web_names) {
      local_higher_attributes[[name]] <- data.frame(
        row.names = colnames(warray)
      )
    }
  }
  if (is.null(local_lower_attributes)) {
    local_lower_attributes <- list()
    for (name in web_names) {
      local_lower_attributes[[name]] <- data.frame(
        row.names = rownames(warray)
      )
    }
  }
  if (is.null(global_meta_attributes)) {
    global_meta_attributes <- list()
  }
  if (is.null(local_meta_attributes)) {
    local_meta_attributes <- list()
    for (name in web_names) {
      local_meta_attributes[[name]] <- list()
    }
  }
  structure(
    warray,
    global_higher_attributes = global_higher_attributes,
    global_lower_attributes  = global_lower_attributes,
    local_higher_attributes  = local_higher_attributes,
    local_lower_attributes   = local_lower_attributes,
    global_meta_attributes   = global_meta_attributes,
    local_meta_attributes    = local_meta_attributes,
    class                    = c("bipartite_webarray", "array")
  )
}

### Access functions --------------------------------------
higher_attributes.bipartite_webarray <- function(x, web_id = NULL) {
  if (is.null(web_id)) {
    return(attr(x, "global_higher_attributes"))
  }  else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    local_attrs <- attr(x, "local_higher_attributes")
    return(local_attrs[[web_id]])
  }
}

lower_attributes.bipartite_webarray <- function(x, web_id = NULL) {
  if (is.null(web_id)) {
    return(attr(x, "global_lower_attributes"))
  }  else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    local_attrs <- attr(x, "local_lower_attributes")
    return(local_attrs[[web_id]])
  }
}

meta_attributes.bipartite_webarray <- function(x, web_id = NULL) {
  if (is.null(web_id)) {
    return(attr(x, "global_meta_attributes"))
  }  else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    local_attrs <- attr(x, "local_meta_attributes")
    return(local_attrs[[web_id]])
  }
}

### Replacement functions ---------------------------------
`higher_attributes<-.bipartite_webarray` <- function(x, web_id = NULL, value) {
  if (!is.data.frame(value)) {
    stop("`higher_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`higher_attributes` must have one row per matrix column")
  }

  if (is.null(web_id)) {
    attr(x, "global_higher_attributes") <- value
  } else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    attr(x, "local_higher_attributes")[[web_id]] <- value
  }
  x
}

`lower_attributes<-.bipartite_webarray` <- function(x, web_id = NULL, value) {
  if (!is.data.frame(value)) {
    stop("`lower_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`lower_attributes` must have one row per matrix column")
  }

  if (is.null(web_id)) {
    attr(x, "global_lower_attributes") <- value
  } else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    attr(x, "local_lower_attributes")[[web_id]] <- value
  }
  x
}

`meta_attributes<-.bipartite_webarray` <- function(x, web_id = NULL, value) {
  if (!is.list(value)) {
    stop("`meta_attributes` must be a list.")
  }

  if (is.null(web_id)) {
    attr(x, "global_meta_attributes") <- value
  } else if (!web_id %in% dimnames(x)[[3]] && !web_id %in% seq_len(dim(x)[[3]])) {
    stop("Index ", web_id, " not in ", as.character(substitute(x)), ".")
  } else {
    attr(x, "local_meta_attributes")[[web_id]] <- value
  }
  x
}

### Subsetting ---------------------------------------------
`[.bipartite_webarray` <- function(x, i, j, k, drop = TRUE) {

  # Subset the actual array
  y <- NextMethod("[", drop = drop)

  global_higher_attributes <- attr(x, "global_higher_attributes")
  local_higher_attributes  <- attr(x, "local_higher_attributes")
  global_lower_attributes  <- attr(x, "global_lower_attributes")
  local_lower_attributes   <- attr(x, "local_lower_attributes")
  global_meta_attributes   <- attr(x, "global_meta_attributes")
  local_meta_attributes    <- attr(x, "local_meta_attributes")

  if (!missing(k)) {
    local_higher_attributes <- local_higher_attributes[k]
    local_lower_attributes  <- local_lower_attributes[k]
    local_meta_attributes   <- local_meta_attributes[k]
    web_indices <- k
  } else {
    web_indices <- dimnames(y)[[3]]
  }

  if (!missing(i)) {
    global_lower_attributes <- global_lower_attributes[i, , drop = FALSE]
    for (web in web_indices) {
      local_lower_attributes[[web]]  <- local_lower_attributes[[web]][i, , drop = FALSE]
    }
  }

  if (!missing(j)) {
    global_higher_attributes <- global_higher_attributes[j, , drop = FALSE]
    for (web in web_indices) {
      local_higher_attributes[[web]]  <- local_higher_attributes[[web]][j, , drop = FALSE]
    }
  }

  ## If k has 1 element only 1 web remains.
  ## Return that web as a bipartite_web object instead of an array.
  if (!missing(k) && length(k) == 1) {
    higher_attributes <- cbind(global_higher_attributes, local_higher_attributes)
    lower_attributes  <- cbind(global_lower_attributes, local_lower_attributes)
    meta_attributes   <- c(global_meta_attributes, local_meta_attributes)
    z <- bipartite_web(
      y,
      higher_attributes = higher_attributes,
      lower_attributes  = lower_attributes,
      meta_attributes   = meta_attributes
    )
    return(z)
  }

  bipartite_webarray(
    y,
    global_higher_attributes = global_higher_attributes,
    global_lower_attributes  = global_lower_attributes,
    local_higher_attributes  = local_higher_attributes,
    local_lower_attributes   = local_lower_attributes,
    global_meta_attributes   = global_meta_attributes,
    local_meta_attributes    = local_meta_attributes,
  )
}

### Print -------------------------------------------------
print.bipartite_webarray <- function(x, ...) {
  cat("<bipartite_webarray object>\n")
  cat("Dimensions:", dim(x)[1L], "x", dim(x)[2L], "x", dim(x) [3L], "\n")
  web_names <- dimnames(x)[[3]]
  cat("Web names:", web_names, "\n")
  cat(strrep("-", 50), "\n")
  cat("Global higher attributes:\n")
  print(higher_attributes(x))
  cat("Global lower attributes:\n")
  print(lower_attributes(x))
  cat("Global meta attributes:\n")
  print(meta_attributes(x))
  nr <- nrow(x)
  nc <- ncol(x)
  for (web in web_names) {
    cat(strrep("-", 50), "\n")
    cat(web, "\n")
    cat("Higher attributes:\n")
    print(higher_attributes(x, web))
    cat("Lower attributes:\n")
    print(lower_attributes(x, web))
    cat("Meta attributes:\n")
    print(meta_attributes(x, web))
    cat("Matrix:\n")
    web_data <- x[,,web]
    web_matrix <- matrix(web_data,
                         nrow = nr,
                         ncol = nc,
                         dimnames = dimnames(web_data))
    print(web_matrix)
  }
}

# warray <- bipartite_webarray(Safariland, vazarr)
# higher_attributes(warray)$rand <- runif(41)
# higher_attributes(warray, "Safariland")$rand <- runif(41)
# warray[1:3, 3:5, "Safariland"]
# print(warray)
# web_names <- dimnames(warray)[[3]]
# web_names

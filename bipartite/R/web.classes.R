###
# @author: Tobias Bauer

# Basic constructor ------------------------------------
bipartite_web <- function(x,
                          higher_attributes = NULL,
                          lower_attributes = NULL,
                          web_attributes = list()) {
  # Validate the underlying object
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Default metadata
  if (is.null(lower_attributes)) {
    lower_attributes <- data.frame(
      row.names = rownames(x)
    )
  }

  if (is.null(higher_attributes)) {
    higher_attributes <- data.frame(
      row.names = colnames(x)
    )
  }

  # Validate metadata dimensions
  if (!is.data.frame(lower_attributes) || nrow(lower_attributes) != nrow(x)) {
    stop("`lower_attributes` must be a data frame with one row per matrix row.")
  }

  if (!is.data.frame(higher_attributes) || nrow(higher_attributes) != ncol(x)) {
    stop("`higher_attributes` must be a data frame with one row per matrix column.")
  }

  if (!is.list(web_attributes)) {
    stop("`info` must be a list.")
  }

  structure(
    x,
    lower_attributes = lower_attributes,
    higher_attributes = higher_attributes,
    web_attributes = web_attributes,
    class = c("bipartite_web", "matrix", "array")
  )
}

# Access functions ----------------------------------
higher_attributes <- function(web) {
  attr(web, "higher_attributes")
}

lower_attributes <- function(web) {
  attr(web, "lower_attributes")
}

# Replacement functions -----------------------------
`higher_attributes<-` <- function(x, value) {
  if (!is.data.frame(value)) {
    stop("`higher_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`higher_attributes` must have one row per matrix column")
  }

  attr(x, "higher_attributes") <- value
  x
}

`lower_attributes<-` <- function(x, value) {
  if (!is.data.frame(value)) {
    stop("`lower_attributes` must be a data frame.")
  }

  if (nrow(value) != ncol(x)) {
    stop("`lower_attributes` must have one row per matrix column")
  }

  attr(x, "lower_attributes") <- value
  x
}

# Print function ------------------------------------
print.bipartite_web <- function(web) {
  cat("<bipartite_web object>\n")
  cat("Dimensions:", nrow(web), "x", ncol(web), "\n")
  cat(strrep("-", 50), "\n")
  cat("Higher species attributes:\n")
  print(higher_attributes(web))
  cat("Lower species attributes:\n")
  print(lower_attributes(web))
  cat(strrep("-", 50), "\n")
  cat("Matrix:\n")
  print(web[, ])
}
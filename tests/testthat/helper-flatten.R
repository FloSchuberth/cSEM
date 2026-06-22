# tests/testthat/helper-flatten.R

# Walk an arbitrarily nested structure and return a named list:
#   names  = dotted/bracketed path to each leaf
#   values = the leaf (atomic vector, matrix, etc.)
flatten_object <- function(x, path = "") {
  # A "leaf" = anything we want to compare directly rather than recurse into.
  is_leaf <- function(v) {
    is.null(v) ||
      is.atomic(v) ||                          # numeric/char/logical vectors & matrices
      is.factor(v) ||
      inherits(v, c("Date", "POSIXct", "formula"))
  }
  
  if (is_leaf(x)) {
    out <- list(x)
    names(out) <- if (path == "") "<root>" else path
    return(out)
  }
  
  # Recurse into lists / S3 / S4 objects
  if (isS4(x)) {
    nms <- slotNames(x)
    children <- setNames(lapply(nms, function(s) slot(x, s)), nms)
  } else {
    children <- as.list(x)
    nms <- names(children)
    if (is.null(nms)) nms <- as.character(seq_along(children))
    nms[nms == ""] <- as.character(which(nms == ""))
  }
  
  out <- list()
  for (i in seq_along(children)) {
    child_path <- if (path == "") nms[i] else paste0(path, "$", nms[i])
    out <- c(out, flatten_object(children[[i]], child_path))
  }
  out
}
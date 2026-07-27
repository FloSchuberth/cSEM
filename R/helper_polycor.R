#' Internal: Fast polychoric correlation
#'
#' Estimate the polychoric correlation between two ordinal (categorical)
#' variables `x` and `y`, i.e., the correlation between the two continuous,
#' bivariate normal latent variables assumed to underlie `x` and `y`
#' \insertCite{Drasgow1988}{cSEM}. Implemented by Kjell S. Slupphaug.
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal
polychor <- function(x, y,
                     control = list(),
                     maxrho =.999,
                     start = rawcor(x, y),
                     thresholds = FALSE) {
  if (!is.integer(x))
    x <- as.integer(as.ordered(x))
  if (!is.integer(y))
    y <- as.integer(as.ordered(y))

  freq <- fastIntTab(x, y)
  zerorows <- rowSums(freq) == 0
  zerocols <- colSums(freq) == 0

  zr <- sum(zerorows)
  zc <- sum(zerocols)

  if (zr > 0) {
    warning2(paste0(
      zr, " row", suffix <- if(zr == 1) "" else "s",
      " with zero marginal", suffix," removed"
    ))
  }

  if (zc > 0) {
    warning2(paste0(
      zc, " column", suffix <- if(zc == 1) "" else "s",
      " with zero marginal", suffix, " removed"
    ))
  }

  freq <- freq[!zerorows, ,drop=FALSE]  
  freq <- freq[, !zerocols, drop=FALSE] 

  r <- nrow(freq)
  c <- ncol(freq)

  if (r < 2) {
    warning2("the cross table has fewer than 2 rows")
    return(NA)
  }

  if (c < 2) {
    warning2("the cross table has fewer than 2 columns")
    return(NA)
  }

  n <- sum(freq)
  rc <- qnorm(cumsum(rowSums(freq))/n)[-r]
  cc <- qnorm(cumsum(colSums(freq))/n)[-c]
  kx <- length(rc)
  ky <- length(cc)

  # We can ignore computing a corner probability if none of the four
  # adjacent cells is nonzero.
  nzero <- freq > 0 & !is.na(freq)
  nz00  <- rbind(cbind(nzero, FALSE), FALSE)   # top-left corner of each nonzero cell
  nz10  <- rbind(cbind(FALSE, nzero), FALSE)   # top-right corner of each nonzero cell
  nz01  <- rbind(FALSE, cbind(nzero, FALSE))   # bottom-left corner of each nonzero cell
  nz11  <- rbind(FALSE, cbind(FALSE, nzero))   # bottom-right corner of each nonzero cell
  keep  <- nz11 | nz10 | nz01 | nz00

  # Keep only interior corners between finite thresholds.
  innerCorners <- keep[seq_len(kx) + 1, seq_len(ky) + 1, drop = FALSE]
  outerCorners <- unname(rbind(
    FALSE, cbind(FALSE, innerCorners, FALSE), FALSE
  ))

  keep.inner.idx <- which(innerCorners)
  keep.outer.idx <- which(outerCorners)
  keep.freq.idx  <- which(nzero)

  pcorners <- unname(rbind(
    0,
    cbind(0, matrix(NA, nrow = kx, ncol = ky), pnorm(rc)),
    c(0, pnorm(cc), 1)
  ))

  gcorners <- unname(rbind(
    0, cbind(0, matrix(NA, nrow = kx, ncol = ky), 0), 0)
  )
    
  cache.rho  <- NA_real_ # for now
  P          <- NULL
  G          <- NULL
  t          <- freq[keep.freq.idx]
  upper.x    <- rep(rc, times = length(cc))[keep.inner.idx]
  upper.y    <- rep(cc, each = length(rc))[keep.inner.idx]
  nr         <- nrow(pcorners)

  freq.row <- ((keep.freq.idx - 1L) %% r) + 1L
  freq.col <- ((keep.freq.idx - 1L) %/% r) + 1L

  # Precompute corners for each nonzero frequency cell.
  idx11 <- freq.row + 1L + freq.col * nr
  idx10 <- freq.row + 1L + (freq.col - 1L) * nr
  idx01 <- freq.row + freq.col * nr
  idx00 <- freq.row + (freq.col - 1L) * nr
  
  updateCache <- function(rho) {
    if (!is.na(cache.rho) && identical(rho, cache.rho))
      return(list(P = P, G = G))
  
    cache.rho <<- rho
  
    # Get probabilities for corners
    pcorners[keep.outer.idx] <- pbivnorm::pbivnorm(
      x   = upper.x,
      y   = upper.y,
      rho = rho
    )

    P <<- pcorners[idx11] - pcorners[idx10] -
      pcorners[idx01] + pcorners[idx00]

    # Get densities for corners
    gcorners[keep.outer.idx] <- dbinorm(
      u   = upper.x,
      v   = upper.y,
      rho = rho,
      force.zero = TRUE # numerical stability
    )
    
    # Get (truncated) densites from corners (for gradient)
    G <<- gcorners[idx11] - gcorners[idx10] -
      gcorners[idx01] + gcorners[idx00]

    list(G = G, P = P)
  }

  polycorObjective <- function(rho) {
    cache <- updateCache(rho = rho)
    -sum(t * log(cache$P), na.rm = TRUE)
  }

  polycorGradient <- function(rho) {
    cache <- updateCache(rho = rho)
    -sum(t * cache$G / cache$P, na.rm = TRUE)
  }
  
  # try 1
  optim <- suppressWarnings(nlminb(
    objective = polycorObjective,
    gradient = polycorGradient,
    start = start, control = control,
    lower = -abs(maxrho), upper = abs(maxrho)
  ))

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- suppressWarnings(nlminb(
      objective = polycorObjective,
      gradient = polycorGradient,
      start = 0.0, control = control,
      lower = -abs(maxrho), upper = abs(maxrho)
    ))
  }

  # check convergence
  if (optim$convergence != 0L) {
    warning2("estimation of polychoric correlation did not converge!")
  }

  rho <- optim$par

  if (thresholds) {
    attr(rho, "thr.x") <- rc
    attr(rho, "thr.y") <- cc
  }

  rho
}


#' Internal: Fast polyserial correlation
#'
#' Estimate the polyserial correlation between a continuous variable `x`
#' and an ordinal (categorical) variable `y`, i.e., the correlation between
#' `x` and the continuous latent variable assumed to underlie `y`
#' \insertCite{Drasgow1988}{cSEM}. Implemented by Kjell S. Slupphaug.
#'
#' @references
#'   \insertAllCited{}
#'    
#' @keywords internal 
polyserial <- function(x, y,
                       control = list(),
                       maxrho =.999,
                       start = rawcor(x, y),
                       thresholds = FALSE) {
  if (!is.integer(y))
    y <- as.integer(as.ordered(y))
  if (!is.numeric(x))
    x <- as.numeric(x)

  freq <- fastIntTab(y) # y is ordinal

  if (length(freq) <= 1)
    stop2("Ordinal variable must have at least two (observed) categories!")

  stopifnot(length(x) == length(y))
  
  n <- length(y)
  thr.inner <- qnorm(cumsum(freq)/n)[-length(freq)]
  thr <- c(-Inf, thr.inner, Inf)
  tau0 <- thr[y]
  tau1 <- thr[y+1]
  finite.lower <- is.finite(tau0)
  finite.upper <- is.finite(tau1)

  cache.rho  <- NA_real_ # for now
  z          <- (x - mean(x)) / sd(x)
  logy       <- NULL
  lower      <- NULL
  upper      <- NULL
  plower     <- NULL
  pupper     <- NULL
  prob       <- NULL
  zrlower    <- NULL
  zrupper    <- NULL
  sigma      <- NULL

  updateCache <- function(rho) {
    if (!is.na(cache.rho) && identical(rho, cache.rho))
      return(invisible(NULL))
  
    cache.rho <<- rho
    ey    <- z * rho
    var   <- 1 - rho * rho
    sigma <<- sqrt(var)
    lower <<- tau0 - ey
    upper <<- tau1 - ey

    plower <<- pnorm(lower, sd = sigma)
    pupper <<- pnorm(upper, sd = sigma)
    prob   <<- pupper - plower

    logy <<- log(prob)

    rho.var <- rho / var
    zrlower <<- z - rho.var * lower
    zrupper <<- z - rho.var * upper

    invisible(NULL)
  }

  plsPolyserialObjective <- function(rho) {
    updateCache(rho)
    -sum(logy)
  }

  plsPolyserialGradient <- function(rho) {
    updateCache(rho)

    lowerTerm <- dnorm(lower, 0, sigma) * zrlower
    upperTerm <- dnorm(upper, 0, sigma) * zrupper

    lowerTerm[!finite.lower] <- 0
    upperTerm[!finite.upper] <- 0

    -sum((lowerTerm - upperTerm) / prob)
  }
  
  # try 1
  optim <- suppressWarnings(nlminb(
    objective = plsPolyserialObjective,
    gradient = plsPolyserialGradient,
    start = start, control = control,
    lower = -abs(maxrho), upper = abs(maxrho)
  ))

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- suppressWarnings(nlminb(
      objective = plsPolyserialObjective,
      gradient = plsPolyserialGradient,
      start = 0.0, control = control,
      lower = -abs(maxrho), upper = abs(maxrho)
    ))
  }

  # check convergence
  if (optim$convergence != 0L)
    warning2("estimation of polyserial correlation did not converge!")

  rho <- optim$par

  if (thresholds)
    attr(rho, "thr.y") <- thr.inner

  rho
}

#' Internal: Bivariate standard normal density
#'
#' Evaluate the density of the bivariate standard normal distribution with
#' correlation `rho` at `(u, v)`. Used internally by [polychor()] to
#' compute the gradient of the polychoric log-likelihood. Implemented by Kjell S. Slupphaug.
#'
#' @references
#'   \insertAllCited{}
#'    
#' @keywords internal 
dbinorm <- function(u, v, rho, force.zero = FALSE, rho.lim = 0.9999) {
  # dirty hack to handle extreme large values for rho
  # note that u, v, and rho are vectorized!
  abs.rho <- abs(rho)
  idx <- which(abs.rho > rho.lim)
  if (length(idx) > 0L)
    rho[idx] <- sign(rho[idx]) * rho.lim

  r <- 1 - rho * rho
  out <- 1 / (2 * pi * sqrt(r)) *
       exp(-0.5 * (u * u - 2 * rho * u * v + v * v) / r)

  # if abs(u) or abs(v) are very large (say, >10), set result equal
  # to exactly zero
  idx <- which(abs(u) > 10 | abs(v) > 10)
  if (length(idx) > 0L && force.zero)
    out[idx] <- 0

  out
}

#' Internal: Raw starting value for [polychor()]/[polyserial()]
#'
#' Compute the ordinary Bravais-Pearson correlation between `x` and `y`
#' after coercing both to numeric. Used as the default starting value for
#' the `rho` optimization in [polychor()] and [polyserial()]; for
#' [polychor()], `x` and `y` are integer category codes rather than the
#' underlying latent variables, so this is only a rough approximation of
#' the polychoric correlation. Implemented by Kjell S. Slupphaug.
#'
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal 
rawcor <- function(x, y) {
  if (!is.numeric(x)) x <- as.numeric(x)
  if (!is.numeric(y)) y <- as.numeric(y)
  cor(x, y)
}

#' Internal: Fast (cross-)tabulation of integer-coded categorical variables
#'
#' A faster, more specialized alternative to [base::table()] for
#' (cross-)tabulating one or two variables that are (or can be coerced to)
#' small positive integer codes, as used by [polychor()] and
#' [polyserial()]. Rows with a missing value in either `x` or `y` are
#' dropped before tabulating.Implemented by Kjell S. Slupphaug.
#'
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal 
fastIntTab <- function(x, y = NULL) {
  if (is.null(y)) {
    ok <- !is.na(x)
    x <- as.integer(x[ok])

    nr <- max(x)
    return(tabulate(x, nbins = nr))
  }

  ok <- !is.na(x) & !is.na(y)
  x <- as.integer(x[ok])
  y <- as.integer(y[ok])

  nr <- max(x)
  nc <- max(y)

  matrix(
    tabulate(x + (y - 1L) * nr, nbins = nr * nc),
    nrow = nr,
    ncol = nc
  )
}

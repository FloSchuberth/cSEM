### R code from vignette source 'libcoin.Rnw'

###################################################
### code chunk number 1: 1dex-1
###################################################
isequal <-
function(a, b) {
    attributes(a) <- NULL
    attributes(b) <- NULL
    if (!isTRUE(all.equal(a, b))) {
        print(a, digits = 10)
        print(b, digits = 10)
        FALSE
    } else
        TRUE
}
library("libcoin")
set.seed(290875)
x <- gl(5, 20)
y <- round(runif(length(x)), 1)
ls1 <- LinStatExpCov(X = model.matrix(~ x - 1), Y = matrix(y, ncol = 1))
ls1$LinearStatistic
tapply(y, x, sum)


###################################################
### code chunk number 2: 1dex-2
###################################################
ls2 <- LinStatExpCov(X = x, Y = matrix(y, ncol = 1))
all.equal(ls1[-grep("Xfactor", names(ls1))],
          ls2[-grep("Xfactor", names(ls2))])


###################################################
### code chunk number 3: 2dex-1
###################################################
X <- rbind(0, diag(nlevels(x)))
ix <- unclass(x)
ylev <- sort(unique(y))
Y <- rbind(0, matrix(ylev, ncol = 1))
iy <- .bincode(y, breaks = c(-Inf, ylev, Inf))
ls3 <- LinStatExpCov(X = X, ix = ix, Y = Y, iy = iy)
all.equal(ls1[-grep("Table", names(ls1))],
          ls3[-grep("Table", names(ls3))])
### works also with factors
ls3 <- LinStatExpCov(X = X, ix = factor(ix), Y = Y, iy = factor(iy))
all.equal(ls1[-grep("Table", names(ls1))],
          ls3[-grep("Table", names(ls3))])


###################################################
### code chunk number 4: 2dex-2
###################################################
ls4 <- LinStatExpCov(ix = ix, Y = Y, iy = iy)
all.equal(ls3[-grep("Xfactor", names(ls3))],
          ls4[-grep("Xfactor", names(ls4))])


###################################################
### code chunk number 5: 2dex-3
###################################################
ls3$Table
xtabs(~ ix + iy)


###################################################
### code chunk number 6: vcov-1
###################################################
ls1$Covariance
vcov(ls1)


###################################################
### code chunk number 7: doTest-1
###################################################
### c_max test statistic
### no p-value
doTest(ls1, teststat = "maximum", pvalue = FALSE)
### p-value
doTest(ls1, teststat = "maximum")
### log(p)-value
doTest(ls1, teststat = "maximum", log = TRUE)
### (1-p)-value
doTest(ls1, teststat = "maximum", lower = TRUE)
### log(1 - p)-value
doTest(ls1, teststat = "maximum", lower = TRUE, log = TRUE)
### quadratic
doTest(ls1, teststat = "quadratic")


###################################################
### code chunk number 8: Contrasts-1
###################################################
set.seed(29)
ls1d <- LinStatExpCov(X = model.matrix(~ x - 1), Y = matrix(y, ncol = 1),
                      nresample = 10, standardise = TRUE)
set.seed(29)
ls1s <- LinStatExpCov(X = as.double(1:5)[x], Y = matrix(y, ncol = 1),
                      nresample = 10, standardise = TRUE)
ls1c <- lmult(1:5, ls1d)
stopifnot(isequal(ls1c, ls1s))
set.seed(29)
ls1d <- LinStatExpCov(X = model.matrix(~ x - 1), Y = matrix(c(y, y), ncol = 2),
                      nresample = 10, standardise = TRUE)
set.seed(29)
ls1s <- LinStatExpCov(X = as.double(1:5)[x], Y = matrix(c(y, y), ncol = 2),
                      nresample = 10, standardise = TRUE)
ls1c <- lmult(1:5, ls1d)
stopifnot(isequal(ls1c, ls1s))


###################################################
### code chunk number 9: ctabsex-1
###################################################
t1 <- ctabs(ix = ix, iy = iy)
t2 <- xtabs(~ ix + iy)
max(abs(t1[-1, -1] - t2))


###################################################
### code chunk number 10: ex-setup
###################################################
N <- 20L
P <- 3L
Lx <- 10L
Ly <- 5L
Q <- 4L
B <- 2L
iX2d <- rbind(0, matrix(runif(Lx * P), nrow = Lx))
ix <- sample(1:Lx, size = N, replace = TRUE)
levels(ix) <- 1:Lx
ixf <- factor(ix, levels = 1:Lx, labels = 1:Lx)
x <- iX2d[ix + 1,]
Xfactor <- diag(Lx)[ix,]
iY2d <- rbind(0, matrix(runif(Ly * Q), nrow = Ly))
iy <- sample(1:Ly, size = N, replace = TRUE)
levels(iy) <- 1:Ly
iyf <- factor(iy, levels = 1:Ly, labels = 1:Ly)
y <- iY2d[iy + 1,]
weights <- sample(0:5, size = N, replace = TRUE)
block <- sample(gl(B, ceiling(N / B))[1:N])
subset <- sort(sample(1:N, floor(N * 1.5), replace = TRUE))
subsety <- sample(1:N, floor(N * 1.5), replace = TRUE)
r1 <- rep(1:ncol(x), ncol(y))
r1Xfactor <- rep(1:ncol(Xfactor), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))
r2Xfactor <- rep(1:ncol(y), each = ncol(Xfactor))


###################################################
### code chunk number 11: Rlibcoin
###################################################
LSEC <-
function(X, Y, weights = integer(0), subset = integer(0), block = integer(0))
{
    if (length(weights) == 0) weights <- rep.int(1, NROW(X))
    if (length(subset) == 0) subset <- seq_len(NROW(X))

    X <- X[subset,, drop = FALSE]
    Y <- Y[subset,, drop = FALSE]
    weights <- weights[subset]

    if (length(block) == 0) {
        w. <- sum(weights)
        wX <- weights * X
        wY <- weights * Y
        ExpX <- colSums(wX)
        ExpY <- colSums(wY) / w.
        CovX <- crossprod(X, wX)
        Yc <- t(t(Y) - ExpY)
        CovY <- crossprod(Yc, weights * Yc) / w.
        T <- crossprod(X, wY)
        Exp <- kronecker(ExpY, ExpX)
        Cov <- w. / (w. - 1) * kronecker(CovY, CovX) -
                1 / (w. - 1) * kronecker(CovY, tcrossprod(ExpX))

        list(LinearStatistic = as.vector(T), Expectation = as.vector(Exp),
             Covariance = Cov, Variance = diag(Cov))
    } else {
        block <- block[subset]
        ret <- list(LinearStatistic = 0, Expectation = 0,
                    Covariance = 0, Variance = 0)
        for (b in levels(block)) {
            tmp <- LSEC(X = X, Y = Y, weights = weights, subset = which(block == b))
            for (l in names(ret)) ret[[l]] <- ret[[l]] + tmp[[l]]
        }
        ret
    }
}


###################################################
### code chunk number 12: cmpr
###################################################
cmpr <-
function(ret1, ret2)
{
    if (inherits(ret1, "LinStatExpCov")) {
        if (!ret1$varonly)
            ret1$Covariance <- vcov(ret1)
    }
    ret1 <- ret1[!sapply(ret1, is.null)]
    ret2 <- ret2[!sapply(ret2, is.null)]
    nm1 <- names(ret1)
    nm2 <- names(ret2)
    nm <- c(nm1, nm2)
    nm <- names(table(nm))[table(nm) == 2]
    isequal(ret1[nm], ret2[nm])
}


###################################################
### code chunk number 13: benchmarks
###################################################
LECVxyws <- LinStatExpCov(x, y, weights = weights, subset = subset)
LEVxyws <- LinStatExpCov(x, y, weights = weights, subset = subset, varonly = TRUE)


###################################################
### code chunk number 14: tests
###################################################
### with X given
testit <-
function(...)
{
    a <- LinStatExpCov(x, y, ...)
    b <- LSEC(x, y, ...)
    d <- LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy, ...)
    cmpr(a, b) && cmpr(d, b)
}

stopifnot(
    testit() && testit(weights = weights) &&
    testit(subset = subset) && testit(weights = weights, subset = subset) &&
    testit(block = block) && testit(weights = weights, block = block) &&
    testit(subset = subset, block = block) &&
    testit(weights = weights, subset = subset, block = block)
)

### without dummy matrix X
testit <-
function(...)
{
    a <- LinStatExpCov(X = ix, y, ...)
    b <- LSEC(Xfactor, y, ...)
    d <- LinStatExpCov(X = integer(0), ix = ix, Y = iY2d, iy = iy, ...)
    cmpr(a, b) && cmpr(d, b)
}

stopifnot(
    testit() && testit(weights = weights) &&
    testit(subset = subset) && testit(weights = weights, subset = subset) &&
    testit(block = block) && testit(weights = weights, block = block) &&
    testit(subset = subset, block = block) &&
    testit(weights = weights, subset = subset, block = block)
)


###################################################
### code chunk number 15: permutations-2d
###################################################
LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy,
              weights = weights, subset = subset, nresample = 10)$PermutedLinearStatistic


###################################################
### code chunk number 16: quadform
###################################################
MPinverse <-
function(x, tol = sqrt(.Machine$double.eps))
{
    SVD <- svd(x)
    pos <- SVD$d > max(tol * SVD$d[1L], 0)
    inv <- SVD$v[, pos, drop = FALSE] %*%
             ((1/SVD$d[pos]) * t(SVD$u[, pos, drop = FALSE]))
    list(MPinv = inv, rank = sum(pos))
}

quadform <-
function(linstat, expect, MPinv)
{
    censtat <- linstat - expect
    censtat %*% MPinv %*% censtat
}

linstat <- ls1$LinearStatistic
expect <- ls1$Expectation
MPinv <- MPinverse(vcov(ls1))$MPinv
MPinv_sym <- MPinv[lower.tri(MPinv, diag = TRUE)]
qf1 <- quadform(linstat, expect, MPinv)
qf2 <- .Call(libcoin:::R_quadform, linstat, expect, MPinv_sym)

stopifnot(isequal(qf1, qf2))


###################################################
### code chunk number 17: ExpectationInfluence
###################################################
sumweights <- sum(weights[subset])
expecty <- colSums(y[subset, ] * weights[subset]) / sumweights

a0 <- expecty
a1 <- .Call(libcoin:::R_ExpectationInfluence, y, weights, subset)
a2 <- .Call(libcoin:::R_ExpectationInfluence, y, as.double(weights), as.double(subset))
a3 <- .Call(libcoin:::R_ExpectationInfluence, y, weights, as.double(subset))
a4 <- .Call(libcoin:::R_ExpectationInfluence, y, as.double(weights), subset)
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset)$ExpectationInfluence

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4) &&
          isequal(a0, a5))


###################################################
### code chunk number 18: CovarianceInfluence
###################################################
sumweights <- sum(weights[subset])
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
a0 <- colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sumweights
a0 <- matrix(a0, ncol = ncol(y))
vary <- diag(a0)

a0 <- a0[lower.tri(a0, diag = TRUE)]
a1 <- .Call(libcoin:::R_CovarianceInfluence, y, weights, subset, 0L)
a2 <- .Call(libcoin:::R_CovarianceInfluence, y, as.double(weights), as.double(subset), 0L)
a3 <- .Call(libcoin:::R_CovarianceInfluence, y, weights, as.double(subset), 0L)
a4 <- .Call(libcoin:::R_CovarianceInfluence, y, as.double(weights), subset, 0L)
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset)$CovarianceInfluence

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4) &&
          isequal(a0, a5))

a0 <- vary
a1 <- .Call(libcoin:::R_CovarianceInfluence, y, weights, subset, 1L)
a2 <- .Call(libcoin:::R_CovarianceInfluence, y, as.double(weights), as.double(subset), 1L)
a3 <- .Call(libcoin:::R_CovarianceInfluence, y, weights, as.double(subset), 1L)
a4 <- .Call(libcoin:::R_CovarianceInfluence, y, as.double(weights), subset, 1L)
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset, varonly = TRUE)$VarianceInfluence

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4) &&
          isequal(a0, a5))


###################################################
### code chunk number 19: ExpectationCovarianceX
###################################################
a0 <- colSums(x[subset, ] * weights[subset])
a1 <- .Call(libcoin:::R_ExpectationX, x, P, weights, subset);
a2 <- .Call(libcoin:::R_ExpectationX, x, P, as.double(weights), as.double(subset))
a3 <- .Call(libcoin:::R_ExpectationX, x, P, weights, as.double(subset))
a4 <- .Call(libcoin:::R_ExpectationX, x, P, as.double(weights), subset)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4) &&
          isequal(a0, LECVxyws$ExpectationX))

a0 <- colSums(x[subset, ]^2 * weights[subset])
a1 <- .Call(libcoin:::R_CovarianceX, x, P, weights, subset, 1L)
a2 <- .Call(libcoin:::R_CovarianceX, x, P, as.double(weights), as.double(subset), 1L)
a3 <- .Call(libcoin:::R_CovarianceX, x, P, weights, as.double(subset), 1L)
a4 <- .Call(libcoin:::R_CovarianceX, x, P, as.double(weights), subset, 1L)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))

a0 <- as.vector(colSums(Xfactor[subset, ] * weights[subset]))
a1 <- .Call(libcoin:::R_ExpectationX, ix, Lx, weights, subset)
a2 <- .Call(libcoin:::R_ExpectationX, ix, Lx, as.double(weights), as.double(subset))
a3 <- .Call(libcoin:::R_ExpectationX, ix, Lx, weights, as.double(subset))
a4 <- .Call(libcoin:::R_ExpectationX, ix, Lx, as.double(weights), subset)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))

a1 <- .Call(libcoin:::R_CovarianceX, ix, Lx, weights, subset, 1L)
a2 <- .Call(libcoin:::R_CovarianceX, ix, Lx, as.double(weights), as.double(subset), 1L)
a3 <- .Call(libcoin:::R_CovarianceX, ix, Lx, weights, as.double(subset), 1L)
a4 <- .Call(libcoin:::R_CovarianceX, ix, Lx, as.double(weights), subset, 1L)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))

r1x <- rep(1:ncol(Xfactor), ncol(Xfactor))
r2x <- rep(1:ncol(Xfactor), each = ncol(Xfactor))
a0 <- colSums(Xfactor[subset, r1x] * Xfactor[subset, r2x] * weights[subset])
a0 <- matrix(a0, ncol = ncol(Xfactor))
vary <- diag(a0)

a0 <- a0[lower.tri(a0, diag = TRUE)]
a1 <- .Call(libcoin:::R_CovarianceX, ix, Lx, weights, subset, 0L)
a2 <- .Call(libcoin:::R_CovarianceX, ix, Lx, as.double(weights), as.double(subset), 0L)
a3 <- .Call(libcoin:::R_CovarianceX, ix, Lx, weights, as.double(subset), 0L)
a4 <- .Call(libcoin:::R_CovarianceX, ix, Lx, as.double(weights), subset, 0L)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 20: SimpleSums
###################################################
a0 <- sum(weights[subset])
a1 <- .Call(libcoin:::R_Sums, N, weights, subset)
a2 <- .Call(libcoin:::R_Sums, N, as.double(weights), as.double(subset))
a3 <- .Call(libcoin:::R_Sums, N, weights, as.double(subset))
a4 <- .Call(libcoin:::R_Sums, N, as.double(weights), subset)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 21: KronSums
###################################################
a0 <- colSums(x[subset, r1] * y[subset, r2] * weights[subset])
a1 <- .Call(libcoin:::R_KronSums, x, P, y, weights, subset, 0L)
a2 <- .Call(libcoin:::R_KronSums, x, P, y, as.double(weights), as.double(subset), 0L)
a3 <- .Call(libcoin:::R_KronSums, x, P, y, weights, as.double(subset), 0L)
a4 <- .Call(libcoin:::R_KronSums, x, P, y, as.double(weights), subset, 0L)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))

a0 <- as.vector(colSums(Xfactor[subset, r1Xfactor] *
                        y[subset, r2Xfactor] * weights[subset]))
a1 <- .Call(libcoin:::R_KronSums, ix, Lx, y, weights, subset, 0L)
a2 <- .Call(libcoin:::R_KronSums, ix, Lx, y, as.double(weights), as.double(subset), 0L)
a3 <- .Call(libcoin:::R_KronSums, ix, Lx, y, weights, as.double(subset), 0L)
a4 <- .Call(libcoin:::R_KronSums, ix, Lx, y, as.double(weights), subset, 0L)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 22: KronSums-Permutation
###################################################
a0 <- colSums(x[subset, r1] * y[subsety, r2])
a1 <- .Call(libcoin:::R_KronSums_Permutation, x, P, y, subset, subsety)
a2 <- .Call(libcoin:::R_KronSums_Permutation, x, P, y, as.double(subset), as.double(subsety))

stopifnot(isequal(a0, a1) && isequal(a0, a2))

a0 <- as.vector(colSums(Xfactor[subset, r1Xfactor] * y[subsety, r2Xfactor]))
a1 <- .Call(libcoin:::R_KronSums_Permutation, ix, Lx, y, subset, subsety)
a2 <- .Call(libcoin:::R_KronSums_Permutation, ix, Lx, y, as.double(subset), as.double(subsety))

stopifnot(isequal(a0, a1) && isequal(a0, a2))


###################################################
### code chunk number 23: colSums
###################################################
a0 <- colSums(x[subset, ] * weights[subset])
a1 <- .Call(libcoin:::R_colSums, x, weights, subset)
a2 <- .Call(libcoin:::R_colSums, x, as.double(weights), as.double(subset))
a3 <- .Call(libcoin:::R_colSums, x, weights, as.double(subset))
a4 <- .Call(libcoin:::R_colSums, x, as.double(weights), subset)

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 24: OneTableSum
###################################################
a0 <- as.vector(xtabs(weights ~ ixf, subset = subset))
a1 <- ctabs(ix, weights = weights, subset = subset)[-1]
a2 <- ctabs(ix, weights = as.double(weights), subset = as.double(subset))[-1]
a3 <- ctabs(ix, weights = weights, subset = as.double(subset))[-1]
a4 <- ctabs(ix, weights = as.double(weights), subset = subset)[-1]

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 25: TwoTableSum
###################################################
a0 <- xtabs(weights ~ ixf + iyf, subset = subset)
class(a0) <- "matrix"
dimnames(a0) <- NULL
attributes(a0)$call <- NULL
a1 <- ctabs(ix, iy, weights = weights, subset = subset)[-1, -1]
a2 <- ctabs(ix, iy, weights = as.double(weights),
            subset = as.double(subset))[-1, -1]
a3 <- ctabs(ix, iy, weights = weights, subset = as.double(subset))[-1, -1]
a4 <- ctabs(ix, iy, weights = as.double(weights), subset = subset)[-1, -1]

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 26: ThreeTableSum
###################################################
a0 <- xtabs(weights ~ ixf + iyf + block, subset = subset)
class(a0) <- "array"
dimnames(a0) <- NULL
attributes(a0)$call <- NULL
a1 <- ctabs(ix, iy, block, weights, subset)[-1, -1,]
a2 <- ctabs(ix, iy, block, as.double(weights), as.double(subset))[-1,-1,]
a3 <- ctabs(ix, iy, block, weights, as.double(subset))[-1,-1,]
a4 <- ctabs(ix, iy, block, as.double(weights), subset)[-1,-1,]

stopifnot(isequal(a0, a1) && isequal(a0, a2) &&
          isequal(a0, a3) && isequal(a0, a4))


###################################################
### code chunk number 27: blocks
###################################################
sb <- sample(block)

ns1 <- do.call(c, tapply(subset, sb[subset], function(i) i))
ns2 <- .Call(libcoin:::R_order_subset_wrt_block, y, integer(0), subset, sb)

stopifnot(isequal(ns1, ns2))


###################################################
### code chunk number 28: kronecker
###################################################
A <- matrix(runif(12), ncol = 3)
B <- matrix(runif(10), ncol = 2)

K1 <- kronecker(A, B)
K2 <- .Call(libcoin:::R_kronecker, A, B)

stopifnot(isequal(K1, K2))


###################################################
### code chunk number 29: MPinv
###################################################
covar <- vcov(ls1)
covar_sym <- ls1$Covariance
n <- (sqrt(1 + 8 * length(covar_sym)) - 1) / 2

tol <- sqrt(.Machine$double.eps)
MP1 <- MPinverse(covar, tol)
MP2 <- .Call(libcoin:::R_MPinv_sym, covar_sym, as.integer(n), tol)

lt <- lower.tri(covar, diag = TRUE)
stopifnot(isequal(MP1$MPinv[lt], MP2$MPinv) &&
          isequal(MP1$rank, MP2$rank))


###################################################
### code chunk number 30: unpack
###################################################
m <- matrix(c(3, 2, 1,
              2, 4, 2,
              1, 2, 5),
            ncol = 3)

s <- m[lower.tri(m, diag = TRUE)]
u1 <- .Call(libcoin:::R_unpack_sym, s, NULL, 0L)
u2 <- .Call(libcoin:::R_unpack_sym, s, NULL, 1L)

stopifnot(isequal(m, u1) && isequal(diag(m), u2))


###################################################
### code chunk number 31: pack
###################################################
m <- matrix(c(4, 3, 2, 1,
              3, 5, 4, 2,
              2, 4, 6, 5,
              1, 2, 5, 7),
            ncol = 4)

s <- m[lower.tri(m, diag = TRUE)]
p <- .Call(libcoin:::R_pack_sym, m)

stopifnot(isequal(s, p))


###################################################
### code chunk number 32: bib
###################################################
thisdir <- getwd()
bibfile <- system.file("REFERENCES.bib", package = "libcoin")
### bibfile may contain spaces LaTeX is unable to deal with on MacOS it seems
if (file.copy(bibfile, to = thisdir, overwrite = TRUE)) {
    bibfile <- "REFERENCES.bib"
} else {
    ### hope for the best
    bibfile <- file.path("..", "inst", "REFERENCES.bib")
}

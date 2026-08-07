# Anatomy of the libcoin calls underneath doTrees(.influence = "mat") ---------
#
# What this script is for: the COIN family of igsca trees does not implement a
# test of its own. doTrees() hands partykit::ctree() a `ytrafo`, ctree hands the
# influence matrix that trafo returns to libcoin, and libcoin does the actual
# testing in C. Two libcoin entry points do all of it:
#
#   LinStatExpCov()  builds the linear statistic T = sum_i X_i (x) Y_i together
#                    with its permutation-null expectation and covariance
#                    (Strasser & Weber 1999). It does NOT test anything.
#   doTest()         takes that object and turns it into a scalar statistic and
#                    a p-value -- either a chi-squared quadratic form over all
#                    columns of Y (variable selection), or a maximally selected
#                    statistic scanning every cutpoint (split selection).
#
# The call chain from the test file:
#
#   grow_tree(influence = "mat", splitter = "native", igsca_tree_control())
#     -> doTrees()                        R/postestimate_doTrees.R:71-224
#        -> partykit::ctree()             R/postestimate_doTrees.R:194
#           -> partykit:::extree_fit()
#              -> partykit:::.ctree_select()      per node
#                 -> partykit:::.ctree_test()     per covariate, SPLITONLY=FALSE
#                    -> .ctree_test_1d() -> .ctree_test_internal()
#                       -> libcoin::LinStatExpCov()  ->  libcoin::doTest()
#              -> partykit:::.ctree_split()       once per node, on the winner
#                 -> .ctree_test(SPLITONLY = TRUE)
#                    -> ... same two libcoin calls, different arguments
#
# How the script is organised:
#
#   Part 0  rebuild doTrees()'s inputs for exactly the fixture used in
#           tests/testthat/test-postestimate_doTrees.R:101
#   Part 1  SPY on a real doTrees() run and record every (LinStatExpCov,
#           doTest) argument list partykit actually produced
#   Part 2  replay the captured root-node selection call by hand, and rebuild
#           every field of the LinStatExpCov output from first principles
#   Part 3  the same for the split-selection call (maxselect)
#   Part 4  the unordered-factor covariate, which takes a third argument shape
#   Part 5  the Monte-Carlo variant (coin_distribution = "approximate")
#   Part 6  argument glossary
#
# Run: Rscript dev/igsca/tree_template/libcoin_call_anatomy.R
# Everything prints; nothing is written to disk.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({
  library(partykit)
  library(libcoin)
})
loaded <- FALSE
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  loaded <- tryCatch(
    {
      suppressMessages(pkgload::load_all(".", quiet = TRUE))
      TRUE
    },
    error = function(e) FALSE
  )
}
if (!loaded) suppressPackageStartupMessages(library(cSEM))

options(width = 120)

hr <- function(title) cat("\n", strrep("=", 78), "\n== ", title, "\n", sep = "")
sub <- function(title) cat("\n-- ", title, " ", strrep("-", max(0, 70 - nchar(title))), "\n", sep = "")


# Part 0: the inputs doTrees() is given ---------------------------------------
#
# Copied verbatim from tests/testthat/test-postestimate_doTrees.R so the numbers
# below are the numbers that test produces.

hr("Part 0 -- rebuilding doTrees()'s inputs")

load(testthat::test_path("data/igscaTrees.Rdata")) # -> dat
covs <- c("z_true", grep("^noise_", names(dat), value = TRUE))

model <- "# Latent variable model
 eta1 =~ x11 + x12 + x13
 eta2 =~ x21 + x22 + x23

 # Composite model
 eta3 <~ x31 + x32 + x33
 eta4 <~ x41 + x42 + x43 + x44

 # Structural model
 eta4 ~ eta3 + eta1
 eta2 ~ eta1 + eta4 + eta3
 "

res <- csem(
  .data = dat,
  .model = model,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = "sum_diff_absolute",
  .iter_max = 100,
  .GSCA_modes = "CCMP",
  .tolerance = 0.0001
)

## doTrees() derives all of this from `.object` alone (postestimate_doTrees.R:56-66).
control    <- igsca_tree_control()
args       <- cSEM:::csem_tree_args(res)
tree_data  <- as.data.frame(args$.data)
indicators <- cSEM:::parseModel(args$.model)$indicators
fml        <- stats::as.formula(paste(
  paste(indicators, collapse = " + "), "~", paste(covs, collapse = " + ")
))

cat("n                 :", nrow(tree_data), "\n")
cat("indicators (J=13) :", paste(indicators, collapse = " "), "\n")
cat("covariates        :", paste(covs, collapse = " "), "\n")
cat("covariate classes :", paste(vapply(tree_data[covs], function(z) class(z)[1], ""),
                                 collapse = " "), "\n")
cat("formula           :", deparse(fml), "\n")

## The ctree_control() doTrees() assembles at postestimate_doTrees.R:145-164.
## Only the fields libcoin ever sees are worth watching:
##   teststat / splitstat -> doTest(teststat=)
##   testtype             -> whether LinStatExpCov(nresample=) is > 0
##   nresample            -> LinStatExpCov(nresample=)
##   minbucket            -> doTest(minbucket=)   (maxselect only)
##   pargs                -> doTest(pargs=)       (maximum-type tests only)
##   nmax = Inf           -> no binning of the covariates
testtype <- c(
  if (isTRUE(control$bonferroni)) "Bonferroni",
  if (control$coin_distribution == "approximate") "MonteCarlo"
)
if (is.null(testtype)) testtype <- "Univariate"

cc <- partykit::ctree_control(
  teststat = "quadratic", splitstat = "quadratic", testtype = testtype,
  nresample = control$R_test, alpha = control$alpha,
  minsplit = control$minsplit, minbucket = control$minbucket,
  maxdepth = control$maxdepth, maxsurrogate = 0L,
  nmax = c(yx = Inf, z = Inf), saveinfo = TRUE
)
cat("\nctree_control() fields libcoin sees:\n")
str(cc[c("teststat", "splitstat", "testtype", "nresample", "minbucket",
         "splittest", "MIA", "tol", "intersplit", "criterion", "logmincriterion")])
cat("\nNote testtype collapsed to \"", cc$testtype, "\" with bonferroni = ", cc$bonferroni,
    ";\nwith coin_distribution = \"asymptotic\" that means nresample is never used.\n", sep = "")

## ctree() turns (fml, data) into an extree_data object before any testing
## happens, and that object is what .ctree_test() indexes with `j`. Rebuilding
## it here gives us the exact X inputs libcoin sees (Parts 3 and 4 need them).
d <- partykit::extree_data(
  fml, data = tree_data, yx = "none", nmax = c(yx = Inf, z = Inf)
)
mf <- model.frame(d)
whichvar <- which(d$variables$z > 0)

sub("what extree_data() prepares, and how `j` indexes it")
cat("`j` indexes the WHOLE model frame, responses first:\n")
cat("  ", paste0(seq_along(names(mf)), "=", names(mf), collapse = "  "), "\n", sep = "")
cat("partitioning variables (extree's `whichvar`): ", paste(whichvar, collapse = " "),
    "  ->  ", paste(names(mf)[whichvar], collapse = " "), "\n", sep = "")
cat("\nPer covariate, the two X shapes .ctree_test_1d() can hand to libcoin:\n")
for (j in whichvar) {
  z <- d$zindex[[j]]
  cat("  j=", j, " ", formatC(names(mf)[j], width = -8),
      " class(d[[j]])=", formatC(class(d[[j]])[1], width = -8),
      "  zindex: ", class(z)[1], " with ", length(attr(z, "levels")), " levels\n", sep = "")
}
cat("\nzindex is present for every covariate even though nmax = Inf -- Inf means\n")
cat("\"do not BIN\", i.e. one level per distinct observed value, not \"do not index\".\n")
cat("d$missings[[j]] is integer(0) throughout (no NAs), d[[j, type = \"scores\"]] is\n")
cat("NULL (no ordinal scores), d$yxindex is NULL (yx = \"none\") and d[[\"(cluster)\"]]\n")
cat("is NULL -- which is why ix/iy/block all arrive at libcoin empty.\n")


# Part 1: spy on the real run --------------------------------------------------
#
# partykit imports LinStatExpCov/doTest into its own imports environment, so
# trace()-ing libcoin alone would not fire. Rebinding in BOTH environments does.
# This is the only way to see the arguments partykit really builds -- everything
# after this part is a replay of what is captured here.

hr("Part 1 -- capturing the real calls made by grow_tree(\"mat\", \"native\")")

CAP <- new.env(parent = emptyenv())
CAP$lin <- list()
CAP$test <- list()
CAP$n_lin <- 0L
CAP$n_test <- 0L
CAP$max <- 32L # the whole tree fits well inside this at these settings

lin_orig  <- libcoin::LinStatExpCov
test_orig <- libcoin::doTest

lin_spy <- function(X, Y, ix = NULL, iy = NULL, weights = integer(0),
                    subset = integer(0), block = integer(0), checkNAs = TRUE,
                    varonly = FALSE, nresample = 0, standardise = FALSE,
                    tol = sqrt(.Machine$double.eps)) {
  out <- lin_orig(X = X, Y = Y, ix = ix, iy = iy, weights = weights,
                  subset = subset, block = block, checkNAs = checkNAs,
                  varonly = varonly, nresample = nresample,
                  standardise = standardise, tol = tol)
  CAP$n_lin <- CAP$n_lin + 1L
  if (CAP$n_lin <= CAP$max) {
    CAP$lin[[CAP$n_lin]] <- list(
      X = X, Y = Y, ix = ix, iy = iy, weights = weights, subset = subset,
      block = block, checkNAs = checkNAs, varonly = varonly,
      nresample = nresample, standardise = standardise, tol = tol, out = out
    )
  }
  out
}

test_spy <- function(object, teststat = c("maximum", "quadratic", "scalar"),
                     alternative = c("two.sided", "less", "greater"),
                     pvalue = TRUE, lower = FALSE, log = FALSE,
                     PermutedStatistics = FALSE, minbucket = 10L,
                     ordered = TRUE, maxselect = object$Xfactor,
                     pargs = libcoin::GenzBretz()) {
  force(maxselect)
  out <- test_orig(object = object, teststat = teststat, alternative = alternative,
                   pvalue = pvalue, lower = lower, log = log,
                   PermutedStatistics = PermutedStatistics, minbucket = minbucket,
                   ordered = ordered, maxselect = maxselect, pargs = pargs)
  CAP$n_test <- CAP$n_test + 1L
  if (CAP$n_test <= CAP$max) {
    CAP$test[[CAP$n_test]] <- list(
      lin_ref = CAP$n_lin, teststat = teststat, alternative = alternative,
      pvalue = pvalue, lower = lower, log = log,
      PermutedStatistics = PermutedStatistics, minbucket = minbucket,
      ordered = ordered, maxselect = maxselect, pargs = pargs, out = out
    )
  }
  out
}

rebind <- function(name, value) {
  envs <- list(asNamespace("libcoin"), parent.env(asNamespace("partykit")))
  for (e in envs) {
    if (!exists(name, envir = e, inherits = FALSE)) next
    if (bindingIsLocked(name, e)) unlockBinding(name, e)
    assign(name, value, envir = e)
  }
}

rebind("LinStatExpCov", lin_spy)
rebind("doTest", test_spy)

set.seed(12353)
tree <- doTrees(
  .object = res, .covariates = covs,
  .influence = "mat", .splitter = "native", .control = control
)

rebind("LinStatExpCov", lin_orig)
rebind("doTest", test_orig)

cat("tree width           :", partykit::width(partykit::node_party(tree)), "\n")
cat("LinStatExpCov calls  :", CAP$n_lin, " (kept the first ", length(CAP$lin), ")\n", sep = "")
cat("doTest calls         :", CAP$n_test, " (kept the first ", length(CAP$test), ")\n", sep = "")

## The two counts differ, so the k-th doTest is NOT the k-th LinStatExpCov.
## Each captured doTest records `lin_ref`, the LinStatExpCov call it consumed.
pair <- vapply(CAP$test, function(x) x$lin_ref, 0L)

## Which covariate a captured call belongs to. Recoverable whenever X carries
## the data values or the factor codes (every selection call); a numeric split
## call passes ranks instead, so it comes back "?" and Part 3 reads the name off
## the fitted split instead.
name_of <- function(L) {
  hit <- vapply(whichvar, function(j) {
    xv <- if (is.numeric(mf[[j]])) as.double(mf[[j]]) else as.double(as.integer(mf[[j]]))
    length(xv) == length(L$X) &&
      isTRUE(all.equal(xv, as.double(L$X), check.attributes = FALSE))
  }, NA)
  if (any(hit)) names(mf)[whichvar][which(hit)[1]] else "?"
}

## At the root, .ctree_select() loops over the three covariates (SPLITONLY =
## FALSE, one call each), then .ctree_split() re-tests the winner in split mode
## (SPLITONLY = TRUE, maxselect = TRUE). Then the same at each child.
sub("argument signature of each captured call")
sig <- do.call(rbind, lapply(seq_along(CAP$lin), function(i) {
  L <- CAP$lin[[i]]
  k <- match(i, pair)
  T <- if (is.na(k)) NULL else CAP$test[[k]]
  data.frame(
    LECV = i,
    doTest = if (is.na(k)) NA_integer_ else k,
    covariate = name_of(L),
    X.class = class(L$X)[1],
    X.nlev = length(attr(L$X, "levels")),
    dim.PxQ = paste(L$out$dimension, collapse = " x "),
    Xfactor = L$out$Xfactor,
    n.subset = length(L$subset),
    varonly = L$varonly,
    teststat = if (is.null(T)) NA_character_ else T$teststat,
    pvalue = if (is.null(T)) NA else T$pvalue,
    maxselect = if (is.null(T)) NA else T$maxselect,
    stringsAsFactors = FALSE
  )
}))
print(sig, row.names = FALSE)

## LinStatExpCov is called more often than doTest, and the gap is not noise:
## .ctree_test_internal() returns list(statistic = NA, p.value = NA) WITHOUT
## calling doTest when every entry of the variance is below ctrl$tol -- i.e.
## when the covariate is constant on that node. That is what happens to z_true
## inside its own children after the root splits on it.
orphans <- setdiff(seq_along(CAP$lin), pair)
if (length(orphans)) {
  sub("LinStatExpCov calls that were never handed to doTest")
  for (i in orphans) {
    L <- CAP$lin[[i]]
    cat("  LECV #", i, ": n = ", length(L$subset),
        ", distinct X values on the subset = ", length(unique(L$X[L$subset])),
        ", all(Variance < ctrl$tol) = ", all(L$out$Variance < cc$tol), "\n", sep = "")
  }
  cat("  -> constant covariate, so the test short-circuits before doTest().\n")
}

cat("\nRoot criteria stored by extree (saveinfo = TRUE):\n")
print(attr(tree, "igsca_info")$root_criteria)
cat("\nRoot split chosen:\n")
print(partykit::node_party(tree)$split)

## Pick the calls the rest of the script dissects, by shape rather than position.
sel_num_k <- which(vapply(seq_along(CAP$test), function(k) {
  !isTRUE(CAP$test[[k]]$maxselect) && !CAP$lin[[pair[k]]]$out$Xfactor
}, NA))[1]
sel_fac_k <- which(vapply(seq_along(CAP$test), function(k) {
  !isTRUE(CAP$test[[k]]$maxselect) && CAP$lin[[pair[k]]]$out$Xfactor
}, NA))[1]
spl_k <- which(vapply(seq_along(CAP$test),
                      function(k) isTRUE(CAP$test[[k]]$maxselect), NA))[1]


# Part 2: the variable-selection call, dissected -------------------------------
#
# This is the call partykit makes once per covariate per node to decide WHICH
# covariate to split on. It answers "is the node's influence Y independent of
# covariate z?", not "where do I cut z?".

hr("Part 2 -- variable selection (SPLITONLY = FALSE, maxselect = FALSE)")

L <- CAP$lin[[pair[sel_num_k]]]
T <- CAP$test[[sel_num_k]]
cat("using LinStatExpCov #", pair[sel_num_k], " -> doTest #", sel_num_k,
    " -- covariate `", name_of(L), "` (numeric)\n", sep = "")

sub("what partykit passed to LinStatExpCov()")
cat("X          : ", class(L$X)[1], " length ", length(L$X),
    "  -- the RAW covariate as doubles.\n", sep = "")
cat("             .ctree_test_1d() uses as.double(x) here, NOT the binned\n")
cat("             index: with maxselect = FALSE libcoin only needs one column.\n")
cat("Y          : ", nrow(L$Y), " x ", ncol(L$Y),
    "  -- model$estfun, i.e. the influence matrix.\n", sep = "")
cat("             This is exactly what doTrees()'s ytrafo returned:\n")
cat("               E  <- calculateGSCAErrors(fit at this node)\n")
cat("               h  <- influence_mat(E) = E^2      (n_node x (J + P))\n")
cat("               ef <- matrix(0, n_full, ncol(h)); ef[subset, ] <- h\n")
cat("             The zero padding is mandatory: .ctree_test() asserts\n")
cat("             NROW(Y) == length(data$zindex[[j]]) == n_full, so Y is\n")
cat("             always full-n and `subset` selects the node's rows.\n")
cat("ix, iy     : ", if (is.null(L$ix)) "NULL" else class(L$ix)[1], " / ",
    if (is.null(L$iy)) "NULL" else class(L$iy)[1],
    " -- non-NULL only on the binned\n", sep = "")
cat("             two-dimensional path (nmax finite); doTrees() sets nmax = Inf.\n")
cat("subset     : integer length ", length(L$subset), " -- the node's row indices,\n", sep = "")
cat("             1-based here (libcoin also accepts 0-based; partykit uses 1-based).\n")
cat("weights    : length ", length(L$weights),
    " -- integer(0) means \"all case weights are 1\".\n", sep = "")
cat("block      : length ", length(L$block),
    " -- ctrl$cluster; integer(0) = no blocking, so\n", sep = "")
cat("             permutations are unrestricted.\n")
cat("nresample  : ", L$nresample,
    " -- 0 because testtype is not \"MonteCarlo\". Nothing is\n", sep = "")
cat("             permuted; the null distribution is the asymptotic one.\n")
cat("varonly    : ", L$varonly,
    " -- the FULL covariance of T is computed. A quadratic-form\n", sep = "")
cat("             test needs it; varonly = TRUE would only give the diagonal.\n")
cat("checkNAs   : ", L$checkNAs,
    " -- partykit removed missings from `subset` already; it\n", sep = "")
cat("             retries with checkNAs = TRUE only if the statistic comes back NA.\n")
cat("standardise: ", L$standardise, " / tol: ", format(L$tol, digits = 4), "\n", sep = "")

sub("what LinStatExpCov() returned")
o <- L$out
cat("dimension            : ", paste(o$dimension, collapse = " x "),
    "  (P = columns of the X design, Q = columns of Y)\n", sep = "")
cat("Xfactor              : ", o$Xfactor, "\n", sep = "")
cat("Sumweights           : ", o$Sumweights, "\n", sep = "")
cat("LinearStatistic      : length ", length(o$LinearStatistic), "\n", sep = "")
cat("Expectation          : length ", length(o$Expectation), "\n", sep = "")
cat("Covariance           : length ", length(o$Covariance),
    "  (lower triangle incl. diagonal, packed)\n", sep = "")
cat("Variance             : length ", length(o$Variance), "\n", sep = "")
cat("ExpectationX         : ", format(o$ExpectationX, digits = 8), "\n", sep = "")
cat("ExpectationInfluence : length ", length(o$ExpectationInfluence), "\n", sep = "")
cat("CovarianceInfluence  : length ", length(o$CovarianceInfluence), "\n", sep = "")
cat("PermutedLinearStatistic: ", paste(dim(o$PermutedLinearStatistic), collapse = " x "),
    "  (empty unless nresample > 0)\n", sep = "")

sub("rebuilding every one of those fields in plain R")
X <- L$X[L$subset]
Yn <- L$Y[L$subset, , drop = FALSE]
n <- length(L$subset)
Ybar <- colMeans(Yn)
Yc <- sweep(Yn, 2, Ybar)
Vy <- crossprod(Yc) / n # ML, not n-1

pack <- function(M) M[lower.tri(M, diag = TRUE)]
mx <- function(a, b) max(abs(a - b))

## Strasser & Weber (1999), for a single-column X (P = 1):
##   T     = sum_i X_i Y_i                            = crossprod(X, Y)
##   mu    = (sum_i X_i) * E[Y]                       with E[Y] = colMeans(Y)
##   Sigma = sw/(sw-1) * V(Y) * sum X^2  -  1/(sw-1) * V(Y) * (sum X)^2
sw <- n
Sigma <- (sw / (sw - 1)) * Vy * sum(X^2) - (1 / (sw - 1)) * Vy * sum(X)^2

cat("LinearStatistic      vs crossprod(X, Y)          : ",
    format(mx(o$LinearStatistic, as.vector(crossprod(X, Yn))), digits = 3), "\n", sep = "")
cat("ExpectationX         vs sum(X)                   : ",
    format(mx(o$ExpectationX, sum(X)), digits = 3), "\n", sep = "")
cat("ExpectationInfluence vs colMeans(Y)              : ",
    format(mx(o$ExpectationInfluence, Ybar), digits = 3), "\n", sep = "")
cat("Expectation          vs sum(X) * colMeans(Y)     : ",
    format(mx(o$Expectation, sum(X) * Ybar), digits = 3), "\n", sep = "")
cat("CovarianceInfluence  vs pack(crossprod(Yc)/n)    : ",
    format(mx(o$CovarianceInfluence, pack(Vy)), digits = 3), "\n", sep = "")
cat("Covariance           vs pack(Strasser-Weber)     : ",
    format(mx(o$Covariance, pack(Sigma)), digits = 3), "\n", sep = "")
cat("Variance             vs diag(Strasser-Weber)     : ",
    format(mx(o$Variance, diag(Sigma)), digits = 3), "\n", sep = "")
cat("\n(ExpectationX is the weighted SUM of X over the subset, not its mean --\n")
cat(" the name is misleading. Same for the Expectation/Covariance pair, which\n")
cat(" are the moments of T itself under the permutation null.)\n")

sub("what partykit passed to doTest(), and what it means")
cat("teststat  : \"", T$teststat, "\"  <- ctrl$teststat. Quadratic form over all ",
    ncol(L$Y), " columns\n", sep = "")
cat("            of Y at once; \"maximum\" would take the max standardised entry.\n")
cat("pvalue    : ", T$pvalue, "   <- !(\"Teststatistic\" %in% ctrl$testtype)\n", sep = "")
cat("lower     : ", T$lower, "   <- return P[T <= t], i.e. 1 - p, not p\n", sep = "")
cat("log       : ", T$log, "   <- ... on the log scale, so the field named\n", sep = "")
cat("            `p.value` in the result actually holds log(1 - p).\n")
cat("ordered   : ", T$ordered, "   <- is.ordered(x) || is.numeric(x); ignored when\n", sep = "")
cat("            maxselect = FALSE.\n")
cat("maxselect : ", T$maxselect, "  <- one test for the whole covariate, no cutpoint scan.\n", sep = "")
cat("minbucket : ", T$minbucket, "  <- ignored when maxselect = FALSE.\n", sep = "")
cat("pargs     : GenzBretz() numerical-integration settings; used only by the\n")
cat("            \"maximum\" test statistic, so inert here.\n")
cat("alternative: \"", T$alternative[1], "\"  <- forced two.sided for quadratic tests.\n", sep = "")

sub("replaying the call verbatim")
lev <- libcoin::LinStatExpCov(
  X = L$X, Y = L$Y, ix = L$ix, iy = L$iy, weights = L$weights,
  subset = L$subset, block = L$block, checkNAs = L$checkNAs,
  varonly = L$varonly, nresample = L$nresample, standardise = L$standardise,
  tol = L$tol
)
tst <- libcoin::doTest(
  lev, teststat = T$teststat, alternative = T$alternative, pvalue = T$pvalue,
  lower = T$lower, log = T$log, PermutedStatistics = T$PermutedStatistics,
  minbucket = T$minbucket, ordered = T$ordered, maxselect = T$maxselect,
  pargs = T$pargs
)
cat("replayed : TestStatistic = ", format(tst$TestStatistic, digits = 10),
    "   p.value(field) = ", format(tst$p.value, digits = 6), "\n", sep = "")
cat("captured : TestStatistic = ", format(T$out$TestStatistic, digits = 10),
    "   p.value(field) = ", format(T$out$p.value, digits = 6), "\n", sep = "")
stopifnot(all.equal(tst, T$out))

sub("and what .ctree_test_internal() does with the result")
cat("It returns list(statistic = log(max(TestStatistic, eps)), p.value = <the field>),\n")
cat("so BOTH numbers reach extree on the log scale:\n")
cat("  statistic = log(", format(T$out$TestStatistic, digits = 8), ") = ",
    format(log(pmax(T$out$TestStatistic, .Machine$double.eps)), digits = 8), "\n", sep = "")
cat("  p.value   = log(1 - p) = ", format(T$out$p.value, digits = 8),
    "   =>  p = ", format(-expm1(T$out$p.value), digits = 6), "\n", sep = "")
cat(".extree_node() then multiplies that log(1 - p) by the number of covariates\n")
cat("tested (the Sidak step, ctrl$bonferroni) and compares against\n")
cat("logmincriterion = log(1 - alpha) = ", format(cc$logmincriterion, digits = 6), ".\n", sep = "")

rc <- attr(tree, "igsca_info")$root_criteria
nm <- name_of(L)
if (!is.null(rc) && nm %in% colnames(rc)) {
  k_tested <- sum(!is.na(rc["criterion", ]))
  cat("\nClosing the loop on what root_criteria(tree) reports for `", nm, "`:\n", sep = "")
  cat("  criterion stored      : ", format(rc["criterion", nm], digits = 10), "\n", sep = "")
  cat("  k * doTest p.value    : ", format(k_tested * T$out$p.value, digits = 10),
      "   (k = ", k_tested, " covariates)\n", sep = "")
  cat("  p.value stored (adj.) : ", format(rc["p.value", nm], digits = 10), "\n", sep = "")
  cat("  1 - exp(criterion)    : ", format(-expm1(rc["criterion", nm]), digits = 10), "\n", sep = "")
  cat("  statistic stored      : ", format(rc["statistic", nm], digits = 10), "\n", sep = "")
  cat("  doTest TestStatistic  : ", format(T$out$TestStatistic, digits = 10), "\n", sep = "")
  cat("So root_criteria() reports the statistic and the ADJUSTED p on their raw\n")
  cat("scales; only the `criterion` row is still log(1 - p), and that is the row\n")
  cat("the stopping rule uses.\n")
}

sub("reconstructing the statistic and the p-value by hand")
unpack <- function(v, q) {
  M <- matrix(0, q, q)
  M[lower.tri(M, diag = TRUE)] <- v
  M + t(M) - diag(diag(M))
}
S <- unpack(o$Covariance, ncol(L$Y))
sv <- svd(S)
## libcoin's Moore-Penrose inverse keeps singular values above tol * max(d);
## its count is the chi-squared degrees of freedom.
keep <- sv$d > sqrt(.Machine$double.eps) * max(sv$d)
df <- sum(keep)
MP <- sv$v[, keep, drop = FALSE] %*% diag(1 / sv$d[keep], df) %*% t(sv$u[, keep, drop = FALSE])
dvec <- o$LinearStatistic - o$Expectation
stat_manual <- as.numeric(t(dvec) %*% MP %*% dvec)
p_manual <- stats::pchisq(stat_manual, df = df, lower.tail = TRUE, log.p = TRUE)
cat("rank(Sigma) = df      : ", df, "   (< ncol(Y) = ", ncol(L$Y),
    " because the P construct\n", sep = "")
cat("                        residual columns are linear combinations of the J\n")
cat("                        indicator residual columns)\n")
cat("manual quadratic form : ", format(stat_manual, digits = 10), "\n", sep = "")
cat("libcoin TestStatistic : ", format(T$out$TestStatistic, digits = 10), "\n", sep = "")
cat("manual log P[X2 <= t] : ", format(p_manual, digits = 10), "\n", sep = "")
cat("libcoin p.value field : ", format(T$out$p.value, digits = 10), "\n", sep = "")


# Part 3: the split-selection call ---------------------------------------------
#
# Same two functions, three arguments different, and a completely different
# question: given the covariate, WHERE is the best cut?

hr("Part 3 -- split selection (SPLITONLY = TRUE, maxselect = TRUE)")

L2 <- CAP$lin[[pair[spl_k]]]
T2 <- CAP$test[[spl_k]]
j2 <- partykit::varid_split(partykit::node_party(tree)$split)
cat("LinStatExpCov #", pair[spl_k], " -> doTest #", spl_k,
    " -- the winning covariate `", names(mf)[j2], "` (", class(mf[[j2]])[1], ")\n", sep = "")

sub("what changed relative to Part 2")
cat("X        : ", class(L2$X)[1], " with ", length(attr(L2$X, "levels")),
    " levels -- NOT raw doubles.\n", sep = "")
cat("           This is data$zindex[[j]]: an integer code into the sorted unique\n")
cat("           values of the covariate, carrying them as attr(X, \"levels\").\n")
cat("           doTest(maxselect = TRUE) does stopifnot(object$Xfactor), so the\n")
cat("           cutpoint scan is reachable ONLY through this encoding -- that is\n")
cat("           the whole reason the split path re-encodes a numeric covariate\n")
cat("           that the selection path passed as plain doubles.\n")
cat("varonly  : ", L2$varonly, " -- SPLITONLY forces this. Cov(T) is P*Q x P*Q, and P is\n", sep = "")
cat("           now the number of levels (up to n), so materialising it is out of\n")
cat("           the question. libcoin does not need it: for each candidate cut it\n")
cat("           rebuilds the 2-sample covariance from CovarianceInfluence, which\n")
cat("           is Q x Q and is stored in full either way (length ",
    length(L2$out$CovarianceInfluence), ").\n", sep = "")
cat("nresample: ", L2$nresample, " -- SPLITONLY never resamples.\n", sep = "")
cat("teststat : \"", T2$teststat, "\" <- ctrl$SPLITstat now, not ctrl$teststat.\n", sep = "")
cat("pvalue   : ", T2$pvalue, " -- a cutpoint does not need a p-value, only an argmax.\n", sep = "")
cat("maxselect: ", T2$maxselect, " -- scan all admissible cuts.\n", sep = "")
cat("minbucket: ", T2$minbucket, " -- NOW it bites: cuts leaving fewer than this on\n", sep = "")
cat("           either side are excluded from the scan.\n")
cat("ordered  : ", T2$ordered, " -- is.ordered(x) || is.numeric(x). TRUE = cut the\n", sep = "")
cat("           levels in order; FALSE = search subsets of an unordered factor.\n")

sub("what came back")
cat("dimension        : ", paste(L2$out$dimension, collapse = " x "),
    "  (P = number of levels now)\n", sep = "")
cat("LinearStatistic  : length ", length(L2$out$LinearStatistic),
    " = P * Q, the per-level column sums of Y,\n", sep = "")
cat("                   which libcoin cumulates to get every candidate split.\n")
cat("doTest returns   :\n")
str(T2$out)

ux2 <- attr(L2$X, "levels")
cat("\n$index is read differently depending on `ordered`:\n")
if (isTRUE(T2$ordered)) {
  sp <- T2$out$index
  brk <- if (isTRUE(cc$intersplit) && sp < length(ux2)) (ux2[sp] + ux2[sp + 1]) / 2 else ux2[sp]
  cat("  ordered = TRUE: a scalar POSITION in the level table.\n")
  cat("    index = ", sp, "  ->  levels[", sp, "] = ", format(brk, digits = 8),
      "   (intersplit = ", cc$intersplit, ")\n", sep = "")
  cat("    tree's breakpoint: ",
      format(partykit::breaks_split(partykit::node_party(tree)$split), digits = 8), "\n", sep = "")
} else {
  cat("  ordered = FALSE: a per-LEVEL 0/1 membership vector, one entry per level.\n")
  cat("    index = (", paste(T2$out$index, collapse = ", "), ") over levels (",
      paste(ux2, collapse = ", "), ")\n", sep = "")
  cat("    partykit builds partysplit(index = as.integer(index) + 1L) with\n")
  cat("    breaks = NULL, i.e. level -> child mapping rather than a cutpoint:\n")
  cat("    tree's split index: (",
      paste(partykit::index_split(partykit::node_party(tree)$split), collapse = ", "),
      "), and breaks_split() is ",
      if (is.null(partykit::breaks_split(partykit::node_party(tree)$split))) "NULL" else "set",
      "\n", sep = "")
  cat("    rows going left   : ",
      sum(as.integer(mf[[j2]]) %in% which(T2$out$index == 0L)), " of ",
      length(L2$subset), "\n", sep = "")
}

sub("replaying it verbatim")
lev2 <- libcoin::LinStatExpCov(
  X = L2$X, Y = L2$Y, ix = L2$ix, iy = L2$iy, weights = L2$weights,
  subset = L2$subset, block = L2$block, checkNAs = L2$checkNAs,
  varonly = L2$varonly, nresample = L2$nresample, standardise = L2$standardise,
  tol = L2$tol
)
tst2 <- libcoin::doTest(
  lev2, teststat = T2$teststat, alternative = T2$alternative, pvalue = T2$pvalue,
  lower = T2$lower, log = T2$log, PermutedStatistics = T2$PermutedStatistics,
  minbucket = T2$minbucket, ordered = T2$ordered, maxselect = T2$maxselect,
  pargs = T2$pargs
)
cat("replayed index (", paste(tst2$index, collapse = ", "), ") / statistic ",
    format(tst2$TestStatistic, digits = 10), "\n", sep = "")
stopifnot(identical(tst2$index, T2$out$index))

## The winning covariate here happens to be an unordered factor, so the ordered
## cutpoint scan -- the common case -- never ran above. Build it by hand for a
## numeric covariate, with the SAME Y and subset the root node used.
sub("the same call for an ORDERED (numeric) covariate, built by hand")
jn <- whichvar[vapply(whichvar, function(j) is.numeric(mf[[j]]), NA)][1]
Yroot <- CAP$lin[[1]]$Y
sub_root <- CAP$lin[[1]]$subset
Xn <- d$zindex[[jn]] # <- exactly what .ctree_test_1d() assigns when SPLITONLY
cat("covariate `", names(mf)[jn], "`  X = zindex with ", length(attr(Xn, "levels")),
    " levels (all values distinct)\n", sep = "")

lev3 <- libcoin::LinStatExpCov(
  X = Xn, Y = Yroot, ix = NULL, iy = NULL, weights = integer(0),
  subset = sub_root, block = integer(0), checkNAs = FALSE,
  varonly = TRUE, nresample = 0L
)
tst3 <- libcoin::doTest(
  lev3, teststat = cc$splitstat, pvalue = FALSE, lower = TRUE, log = TRUE,
  ordered = TRUE, maxselect = TRUE, minbucket = cc$minbucket, pargs = cc$pargs
)
ux3 <- attr(Xn, "levels")
sp3 <- tst3$index
brk3 <- if (isTRUE(cc$intersplit) && sp3 < length(ux3)) (ux3[sp3] + ux3[sp3 + 1]) / 2 else ux3[sp3]
cat("dimension ", paste(lev3$dimension, collapse = " x "),
    "  LinearStatistic length ", length(lev3$LinearStatistic), "\n", sep = "")
cat("doTest -> index ", sp3, "  statistic ", format(tst3$TestStatistic, digits = 8),
    "  cut at levels[", sp3, "] = ", format(brk3, digits = 8), "\n", sep = "")
cat("rows going left: ", sum(mf[[jn]] <= brk3), " of ", length(sub_root), "\n", sep = "")

## Hard check: partykit's own SPLITONLY path, which is these two calls plus the
## index -> breakpoint mapping, must land on the same number.
pk_split <- partykit:::.ctree_test(
  model = list(estfun = Yroot), trafo = NULL, data = d, subset = sub_root,
  weights = integer(0), j = jn, SPLITONLY = TRUE, ctrl = cc
)
cat("partykit:::.ctree_test(SPLITONLY = TRUE) breakpoint: ",
    format(partykit::breaks_split(pk_split), digits = 8), "\n", sep = "")
stopifnot(isTRUE(all.equal(brk3, partykit::breaks_split(pk_split))))

sub("the effect of minbucket, re-running that LECV with other values")
for (mb in c(1L, 30L, 200L, 450L)) {
  t_mb <- libcoin::doTest(lev3, teststat = "quadratic", pvalue = FALSE, lower = TRUE,
                          log = TRUE, ordered = TRUE, maxselect = TRUE, minbucket = mb)
  cat("  minbucket = ", formatC(mb, width = 3), " -> index ", formatC(t_mb$index, width = 4),
      "  cut at ", formatC(ux3[t_mb$index], width = 10, format = "f", digits = 5),
      "  statistic ", format(t_mb$TestStatistic, digits = 7), "\n", sep = "")
}


# Part 4: the unordered-factor covariate ---------------------------------------
#
# z_true in this fixture is a two-level factor, so it takes a third argument
# shape: X is the zindex even during VARIABLE SELECTION, and ordered = FALSE.

hr("Part 4 -- a factor covariate takes X = zindex even for selection")

if (!is.na(sel_fac_k)) {
  Lf <- CAP$lin[[pair[sel_fac_k]]]
  Tf <- CAP$test[[sel_fac_k]]
  fname <- name_of(Lf)
  cat("LinStatExpCov #", pair[sel_fac_k], " -> doTest #", sel_fac_k, " -- `", fname,
      "` is a ", class(mf[[fname]])[1], " with levels ",
      paste(levels(mf[[fname]]), collapse = "/"), "\n", sep = "")
  cat("X          : ", class(Lf$X)[1], ", ", length(attr(Lf$X, "levels")),
      " levels -> Xfactor = ", Lf$out$Xfactor, "\n", sep = "")
  cat("dimension  : ", paste(Lf$out$dimension, collapse = " x "),
      "  (P = 2 dummy columns, one per level)\n", sep = "")
  cat("ordered    : ", Tf$ordered, "  <- is.ordered(x) || is.numeric(x) is FALSE here\n", sep = "")
  cat("maxselect  : ", Tf$maxselect, "  <- still a plain selection test\n", sep = "")
  cat("statistic  : ", format(Tf$out$TestStatistic, digits = 10),
      "   log(1 - p) = ", format(Tf$out$p.value, digits = 6), "\n", sep = "")
  cat("\nSo there are three argument shapes in total, not two:\n")
  cat("  numeric covariate, selection : X = as.double(x),    Xfactor FALSE, ordered TRUE\n")
  cat("  numeric covariate, split     : X = zindex (n lev),  Xfactor TRUE,  ordered TRUE\n")
  cat("  factor covariate,  either    : X = zindex (k lev),  Xfactor TRUE,  ordered FALSE\n")
  cat("Only `maxselect` distinguishes selection from split for a factor -- which is\n")
  cat("why the two z_true rows in the Part 1 table share a statistic (",
      format(Tf$out$TestStatistic, digits = 7), ").\n", sep = "")
} else {
  cat("(no factor covariate among the captured calls)\n")
}


# Part 5: the Monte-Carlo variant ----------------------------------------------
#
# igsca_tree_control(coin_distribution = "approximate") puts "MonteCarlo" in
# testtype, which is the ONLY thing that makes nresample non-zero.

hr("Part 5 -- nresample > 0 (coin_distribution = \"approximate\")")

lev_mc <- libcoin::LinStatExpCov(
  X = L$X, Y = L$Y, subset = L$subset, weights = L$weights, block = L$block,
  nresample = 500L, varonly = FALSE, checkNAs = FALSE
)
cat("PermutedLinearStatistic: ", paste(dim(lev_mc$PermutedLinearStatistic), collapse = " x "),
    "  (P*Q rows x nresample columns)\n", sep = "")
cat("Each column is the linear statistic under one permutation of the subset\n")
cat("row order, drawn inside libcoin's C code -- this is where R_test goes on\n")
cat("the COIN path, and it is why R_test is inert under \"asymptotic\".\n")
tst_mc <- libcoin::doTest(lev_mc, teststat = "quadratic", pvalue = TRUE, lower = TRUE,
                          log = TRUE, ordered = TRUE, maxselect = FALSE,
                          minbucket = cc$minbucket, pargs = cc$pargs)
cat("\nasymptotic : statistic ", format(T$out$TestStatistic, digits = 8),
    "   p = ", format(-expm1(T$out$p.value), digits = 6), "\n", sep = "")
cat("MonteCarlo : statistic ", format(tst_mc$TestStatistic, digits = 8),
    "   p = ", format(-expm1(tst_mc$p.value), digits = 6),
    "   (same statistic, permutation p)\n", sep = "")
cat("\nblock = would restrict those permutations to within-block; partykit passes\n")
cat("ctrl$cluster there, which doTrees() never sets, hence integer(0).\n")


# Part 6: argument glossary ----------------------------------------------------

hr("Part 6 -- glossary")

cat("
LinStatExpCov(X, Y, ix, iy, weights, subset, block, checkNAs, varonly,
              nresample, standardise, tol)

  X          n x P design for the covariate. Doubles for a numeric covariate
             under variable selection; an `enum integer` (zindex, levels =
             sorted unique values) for factors and for every split scan.
  Y          n x Q influence matrix. Here: casewise squared GSCA residuals,
             influence_mat(E) = E^2, zero-padded to full n by doTrees()'s ytrafo.
             Q = J indicators + P constructs = 17; only 13 dimensions are
             linearly independent, which is the chi-squared df.
  ix, iy     integer index vectors for the binned 2-d path. NULL under
             nmax = Inf, which doTrees() always sets.
  weights    integer case weights; integer(0) == all ones.
  subset     1-based row indices defining the node. THE node boundary lives
             here, not in X or Y -- both are always full-n.
  block      blocking/cluster factor restricting permutations; integer(0) = none.
  checkNAs   scan for NAs. FALSE on the first try because partykit already
             stripped them from `subset`.
  varonly    TRUE = diagonal of Cov(T) only. Set on every split scan, where the
             full P x Q covariance would be enormous and is not needed.
  nresample  number of permutations to draw. > 0 only under testtype MonteCarlo.
  standardise standardise the permuted statistics; partykit leaves it FALSE.
  tol        numerical zero for variances.

  Returned: LinearStatistic / Expectation / Covariance / Variance are T and its
  permutation-null moments (Covariance and CovarianceInfluence are packed lower
  triangles). ExpectationX is the weighted SUM of X columns. Expectation-
  Influence / CovarianceInfluence are the moments of Y alone. Xfactor drives
  whether doTest(maxselect = TRUE) is even legal.

doTest(object, teststat, alternative, pvalue, lower, log, PermutedStatistics,
       minbucket, ordered, maxselect, pargs)

  object     the LinStatExpCov() result, unmodified.
  teststat   \"quadratic\" (chi-squared form over all Q columns, what doTrees()
             uses), \"maximum\" (max standardised entry), \"scalar\" (P*Q == 1).
  alternative forced \"two.sided\" for quadratic and for maxselect.
  pvalue     compute a p-value at all. FALSE on the split scan.
  lower      TRUE -> report P[T <= t] instead of the upper-tail p-value.
  log        TRUE -> on the log scale. lower + log together mean the returned
             `p.value` field is log(1 - p). partykit keeps it that way: the
             Sidak adjustment is a multiplication by k on this scale, and the
             stop rule is a comparison against log(1 - alpha).
  PermutedStatistics return the per-permutation statistics too.
  minbucket  minimum node size on either side of a candidate cut. Only read
             when maxselect = TRUE.
  ordered    TRUE -> cutpoints respect level order; FALSE -> search level
             subsets of an unordered factor.
  maxselect  FALSE -> one test for the covariate (variable selection).
             TRUE  -> scan every admissible cutpoint and return $index, the
             argmax. Requires object$Xfactor.
  pargs      GenzBretz() integration settings for the multivariate normal
             probabilities behind teststat = \"maximum\"; unused for quadratic.
")

hr("done")

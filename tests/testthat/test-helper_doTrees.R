# Tests for R/helper_doTrees.R. 
model_pls <- "
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33

eta2 ~ eta1
eta3 ~ eta1 + eta2
"

res_pls <- csem(.data = threecommonfactors, .model = model_pls)

test_that("bdiagFit() repeats a single-group implied VCV block-diagonally", {
  for(type in c("indicator", "construct")) {
    Sigma <- fit(res_pls, .saturated = FALSE, .type_vcv = type)
    B3    <- bdiagFit(res_pls, .n_blocks = 3, .type_vcv = type)

    # One comparison covers dimensions, block content, and zero off-blocks
    expect_equal(unname(B3), kronecker(diag(3), unname(Sigma)))
    # Names: block-suffixed and unique
    expect_identical(rownames(B3)[seq_len(nrow(Sigma))],
                     paste0(rownames(Sigma), "_1"))
    expect_identical(anyDuplicated(rownames(B3)), 0L)
  }
})

test_that("bdiagFit() block-diagonalizes multigroup fits in group order", {
  # Model syntax and csem() arguments as in tests/testthat/test-csem_fit.R
  model_IGSCA <- "
  OrgPres =~ cei1 + cei2 + cei3
  OrgIden <~ ma1 + ma2 + ma3
  AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
  AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

  OrgIden ~ OrgPres
  AffLove ~ OrgIden
  AffJoy  ~ OrgIden"

  res_mg <- csem(
    .data = BergamiBagozzi2000,
    .model = model_IGSCA,
    .approach_weights = "GSCA",
    .id = "gender",
    .tolerance = 1e-4,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "NCMP",
    .iter_max = 600
  )

  Sigma_list <- fit(res_mg, .saturated = FALSE, .type_vcv = "indicator")
  B          <- bdiagFit(res_mg, .type_vcv = "indicator")

  expect_identical(dim(B), rep(sum(sapply(Sigma_list, nrow)), 2))
  pos <- 0
  for(g in names(Sigma_list)) {
    idx <- pos + seq_len(nrow(Sigma_list[[g]]))
    expect_equal(B[idx, idx], Sigma_list[[g]], ignore_attr = TRUE)
    pos <- pos + nrow(Sigma_list[[g]])
  }

  # .n_blocks is ignored (with a warning) for multigroup objects
  expect_warning(bdiagFit(res_mg, .n_blocks = 2, .type_vcv = "indicator"))
})

test_that("bdiagFit() validates .n_blocks", {
  expect_error(bdiagFit(res_pls, .n_blocks = 0))
  expect_error(bdiagFit(res_pls, .n_blocks = 1.5))
})

# doTrees node-level helpers ---------------------------------------------
# These call the partition machinery directly, one node at a time. At the tree
# level the root gives a kernel no freedom -- z_true is a 2-level factor, so
# candidate_partitions() returns exactly one candidate -- so this is the only
# layer where cutpoint choice and the scan cache can be observed at all.
load(testthat::test_path("data/igscaTrees.Rdata")) # Creates dat

dat$noise_1 <- round(dat$noise_1)
dat$noise_2 <- cut(
  dat$noise_2,
  breaks = stats::quantile(dat$noise_2, probs = seq(0, 1, length.out = 21L)),
  include.lowest = TRUE,
  ordered_result = TRUE
)

model_trees <- "# Latent variable model
 eta1 =~ x11 + x12 + x13
 eta2 =~ x21 + x22 + x23

 # Composite model
 eta3 <~ x31 + x32 + x33
 eta4 <~ x41 + x42 + x43 + x44

 # Structural model
 eta4 ~ eta3 + eta1
 eta2 ~ eta1 + eta4 + eta3
 "

# The argument list doTrees() replays at every node, as csem() stores it.
args_trees <- csem(
  .data = dat,
  .model = model_trees,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = "sum_diff_absolute",
  .iter_max = 100,
  .GSCA_modes = "CCMP",
  .tolerance = 0.0001
)$Information$Arguments

# One pooled node fit, packaged the way the trafo hands it to a kernel.
node_fixture <- function(minbucket = 100L) {
  ind <- parseModel(args_trees$.model)$indicators
  ft <- try_fit(dat[, ind, drop = FALSE], args_trees)
  list(
    model = list(object = ft$fit),
    mf = dat,
    subset = seq_len(nrow(dat)),
    ctrl = list(
      collector = new_collector(),
      args = args_trees,
      indicators = ind,
      minbucket = minbucket
    )
  )
}

test_that("fit_csem() replays the arguments of the fit it was given", {
  ind <- parseModel(args_trees$.model)$indicators
  node <- fit_csem(dat[1:500, ind, drop = FALSE], args_trees)
  used <- node$Information$Arguments
  expect_identical(used$.approach_weights, "GSCA")
  expect_identical(used$.GSCA_modes, args_trees$.GSCA_modes)
  expect_identical(used$.tolerance, args_trees$.tolerance)
  # .id is dropped, not passed as NULL, so csem()'s own default applies
  expect_s3_class(node, "cSEMResults_default")

  # The MGA fit the split kernels compare is the same call with .id set
  g <- cbind(
    group = factor(rep(1:2, each = 500L)),
    dat[, ind, drop = FALSE]
  )
  expect_s3_class(fit_csem(g, args_trees, .id = "group"), "cSEMResults_multi")
})

# A kernel that throws is caught by argmax_split() and turned into NA, which
# partykit reads as "no admissible split" -- so at the tree level a broken kernel
# and a genuinely unsplittable node look identical. Called directly, a kernel
# either returns a finite statistic or it does not, and the failure names it.
test_that("every split kernel returns a finite statistic", {
  fx <- node_fixture()
  # A median split on a noise covariate: no group difference is expected here,
  # but every kernel must still put a number on it.
  goes_left <- fx$mf$noise_1 <= stats::median(fx$mf$noise_1)
  for (k in c("FITdiff", "DLi", "DGi")) {
    expect_true(
      is.finite(
        partition_stat(stat_kind = k, 
          model = fx$model, 
          mf = fx$mf,
          subset = fx$subset,
          goes_left = goes_left,
          ctrl = fx$ctrl)
      ),
      info = k
    )
  }
})

test_that("partition_stat() records a failed auxiliary fit instead of throwing", {
  local_mocked_bindings(
    try_fit = function(.data, .args, .id = NULL) list(fit = NULL, ok = FALSE)
  )
  coll <- new_collector()
  ctrl <- list(collector = coll, args = args_trees, indicators = "x11")
  goes_left <- rep(c(TRUE, FALSE), length.out = nrow(dat))

  # An unfittable candidate partition must come back as NA, not an exception:
  # anything that escapes here would make SimDesign redraw the whole rep.
  expect_identical(
    partition_stat(
      stat_kind = "DLi",
      model = list(object = NULL),
      mf = dat,
      subset = seq_len(nrow(dat)),
      goes_left = goes_left,
      ctrl = ctrl
    ),
    NA_real_
  )
  expect_identical(coll$n_fail_candidate, 1L)
})

# ndt_dists() selectivity -------------------------------------------------
# The point of `dists` is negative -- what is NOT computed -- so it is mocked
# rather than fitted: a call that never happens is only observable from inside.
test_that("ndt_dists() computes only the distance it is asked for", {
  seen <- character(0)
  local_mocked_bindings(
    bdiagFit = function(.object = NULL, .n_blocks = 1L,
                        .type_vcv = "indicator", .saturated = FALSE) {
      seen <<- c(seen, paste0("bdiag_", .type_vcv))
      diag(2)
    },
    calculateDG = function(...) { seen <<- c(seen, "DG"); 1 },
    calculateDL = function(...) { seen <<- c(seen, "DL"); 2 }
  )
  d <- ndt_dists(diag(2), diag(2), mga_fit = NULL, dists = "DLi")

  # Still the full named double vector: partition_stat() subscripts it by name,
  # and NA_real_ (not a bare logical NA) is what argmax_split()'s
  # vapply(..., numeric(1)) accepts.
  expect_identical(names(d), c("DGc", "DGi", "DLc", "DLi"))
  expect_identical(typeof(d), "double")
  expect_identical(unname(d[["DLi"]]), 2)
  expect_true(all(is.na(d[c("DGc", "DGi", "DLc")])))
  # calculateDG() was never reached and the construct VCV was never built.
  expect_identical(seen, c("bdiag_indicator", "DL"))
})

test_that("a complex geodesic distance stumps the scan instead of killing the tree", {
  # calculateDG() eigen-decomposes a non-symmetric product, so it can return a
  # complex scalar with no error at all. Unguarded it coerced the whole distance
  # vector to complex, and vapply(..., numeric(1)) in argmax_split() rejects
  # that AFTER the per-candidate tryCatch has returned -- so the type error
  # escaped into partykit rather than degrading to "no admissible split".
  fx <- node_fixture()
  local_mocked_bindings(
    calculateDG = function(...) complex(real = 1, imaginary = 2)
  )
  expect_no_error(
    sp <- argmax_split(
      splitter = split_max_dgi,
      collector = fx$ctrl$collector,
      model = fx$model,
      trafo = NULL,
      mf = fx$mf,
      subset = fx$subset,
      whichvar = match("noise_1", names(fx$mf)),
      ctrl = fx$ctrl
    )
  )
  expect_null(sp) # nothing finite: no split, but a diagnosed one
  expect_gt(fx$ctrl$collector$n_fail_candidate, 0L)
})

test_that("the DLi path never reaches calculateDG()", {
  fx <- node_fixture()
  local_mocked_bindings(
    calculateDG = function(...) stop("calculateDG() must not be called")
  )
  goes_left <- fx$mf$noise_1 <= stats::median(fx$mf$noise_1)
  expect_true(is.finite(
    partition_stat(
      stat_kind = "DLi", model = fx$model, mf = fx$mf,
      subset = fx$subset, goes_left = goes_left, ctrl = fx$ctrl
    )
  ))
})

test_that("partition_stat() refuses an unknown stat_kind loudly", {
  # ds[stat_kind] turned a typo into a silent NA, which partykit reads as "no
  # admissible split". The check has to sit ahead of partition_stat()'s own
  # tryCatch, or it is swallowed into an NA and a bogus counter increment.
  fx <- node_fixture()
  expect_error(
    partition_stat(
      stat_kind = "DLx", model = fx$model, mf = fx$mf,
      subset = fx$subset, goes_left = rep(TRUE, length(fx$subset)),
      ctrl = fx$ctrl
    ),
    "stat_kind"
  )
  expect_identical(fx$ctrl$collector$n_fail_candidate, 0L)
})

# This is where cutpoint choice has to be tested: at the tree level the root
# offers no freedom, since z_true is a 2-level factor and candidate_partitions()
# returns exactly one candidate for it. noise_1 is numeric, and rounding leaves
# it seven distinct values -- two admissible cuts at minbucket = 100 -- so a
# kernel that runs but picks the wrong cutpoint is visible here and nowhere
# else. Two is thin: it is enough for the kernels to disagree on this fixture,
# but a coarser noise_1 would leave them nothing to disagree about.
test_that("argmax_split() returns an admissible cutpoint the kernels disagree on", {
  fx <- node_fixture()
  j <- match("noise_1", names(fx$mf)) # column number of noise_1 in the data
  z <- fx$mf[[j]]
  # Function for getting the cutpoint in mf$j to split the data on
  cutpoint <- function(kern) {
    sp <- argmax_split(
      splitter = kern,
      collector = fx$ctrl$collector,
      model = fx$model,
      mf = fx$mf,
      subset = fx$subset,
      whichvar = j,
      ctrl = fx$ctrl
    )
    # NULL is what argmax_split() returns when no candidate produced a finite
    # statistic -- the dead kernel this layer exists to localise.
    expect_true(inherits(sp, "partysplit"))
    partykit::breaks_split(sp)
  }

  brs <- vapply(
    list(FIT = split_max_fitdiff, DLi = split_max_dli, DGi = split_max_dgi),
    cutpoint,
    numeric(1)
  )
  for (k in names(brs)) {
    # Strictly below the maximum: the top value sends every row left, which
    # partykit accepts and grows a degenerate child from.
    expect_true(is.finite(brs[[k]]) && brs[[k]] < max(z), info = k)
    # ctree's default reports the observed value below the gap, so a break
    # that is not one of the covariate's own values is the midpoint rule.
    expect_true(brs[[k]] %in% z, info = k)
  }

  # The kernels are only worth pairing with a selector if they disagree about
  # where to cut. Asserted for FIT vs DLi only: DLi and DGi are both distances
  # between the same model-implied indicator VCVs and may legitimately land on
  # the same candidate, so distinctness across all three would be an assertion
  # about this fixture rather than about the kernels.
  expect_false(isTRUE(all.equal(brs[["FIT"]], brs[["DLi"]])))
})

# candidate_partitions() conventions vs ctree ----------------------------
# These are pure-function tests: the point is that the candidates handed to a
# kernel, and the partysplit handed back to partykit, follow the same
# conventions partykit:::.ctree_test_internal() follows on the native route.

test_that("candidate_partitions() breaks at observed values, as ctree does", {
  z <- c(1, 2, 4, 8)
  brk <- function(cands) {
    vapply(cands, function(cc) partykit::breaks_split(cc$split), numeric(1))
  }
  # ctree's default is intersplit = FALSE, which reports the observed value
  # below the gap rather than the midpoint.
  obs <- candidate_partitions(1L, z, z, minbucket = 1L)
  expect_identical(brk(obs), c(1, 2, 4))
})

test_that("candidate_partitions() emits splits partykit routes as told", {
  mf <- data.frame(
    num = c(1, 2, 4, 8, 16, 32),
    ord = ordered(c("a", "b", "c", "d", "e", "f")),
    fac = factor(c("a", "b", "c", "a", "b", "c"))
  )
  # minbucket = 2L leaves every branch non-empty but does make keep_min() drop
  # candidates, so the routing convention is checked on a filtered list too --
  # at minbucket = 1L nothing is ever filtered and the guard goes unexercised.
  for (mb in c(1L, 2L)) {
    for (j in seq_along(mf)) {
      z <- mf[[j]]
      cands <- candidate_partitions(j, z, z, minbucket = mb)
      expect_gt(length(cands), 0L)
      for (cc in cands) {
        # The kernel is told `goes_left`; partykit later routes rows through
        # kidids_split(). A candidate whose two disagree scores one partition
        # and grows another.
        expect_identical(
          partykit::kidids_split(cc$split, mf) == 1L,
          cc$goes_left,
          info = paste0(names(mf)[j], ", minbucket = ", mb)
        )
      }
    }
  }
  # keep_min() bit: at minbucket = 2 the two extreme cuts are inadmissible.
  expect_length(candidate_partitions(1L, mf$num, mf$num, minbucket = 1L), 5L)
  expect_length(candidate_partitions(1L, mf$num, mf$num, minbucket = 2L), 3L)

  # ctree sets index = 1:2 on every break-valued split it emits.
  for (v in c("num", "ord")) {
    j <- match(v, names(mf))
    sp <- candidate_partitions(j, mf[[j]], mf[[j]], minbucket = 1L)[[1L]]$split
    expect_identical(partykit::index_split(sp), 1:2, info = v)
  }
})

test_that("candidate_partitions() enumerates the node's levels, not the column's", {
  # `z` and `zs` differ at every node below the root: a level is dropped by a
  # parent split, or lost with the covariate's missing rows in argmax_split().
  # Sized off nlevels(z) rather than the node's own levels, this node costs
  # 2^15 - 1 candidates -- each one a two-group IGSCA fit -- for the three
  # bipartitions it actually has, and slips past the K > 11L guard doing it,
  # since that guard counts the levels the node has.
  z  <- factor(letters[1:16], levels = letters[1:16])
  zs <- factor(rep(letters[1:3], each = 40L), levels = letters[1:16])
  cands <- candidate_partitions(1L, z, zs, minbucket = 1L)
  expect_length(cands, 3L) # 2^(3 - 1) - 1

  # Levels absent from the node stay NA, which is how partykit says "route this
  # one by `prob`". A real kid id would send a level never seen here down an arm
  # picked by a bit pattern.
  idx <- partykit::index_split(cands[[1L]]$split)
  expect_identical(idx[4:16], rep(NA_integer_, 13L))
  expect_true(all(idx[1:3] %in% 1:2))

  # The routing convention has to survive z != zs, which is the case the
  # same-vector fixtures above cannot see.
  for (cc in cands) {
    expect_identical(
      partykit::kidids_split(cc$split, data.frame(f = zs)) == 1L,
      cc$goes_left
    )
  }
})

test_that("candidate_partitions() refuses an unordered factor it cannot afford", {
  # 11 levels is 1023 bipartitions -- the last scan a kernel can be asked for.
  ok <- factor(rep(letters[1:11], each = 2L))
  expect_length(candidate_partitions(1L, ok, ok, minbucket = 1L), 1023L)

  too_many <- factor(rep(letters[1:12], each = 2L))
  expect_error(
    candidate_partitions(1L, too_many, too_many, minbucket = 1L),
    "unordered factor"
  )
})

test_that("argmax_split() scans a covariate that has missing values", {
  # partykit drops the covariate's own missing rows before testing it
  # (`subsetNArm` in partykit:::.ctree_test()) and routes them by the split's
  # `prob` afterwards. Left in, `zs <= ct` is NA, sum() is NA, and Filter()
  # silently drops every candidate -- the covariate goes unsplit with no
  # diagnostic at all.
  mf <- data.frame(z = c(rep(1:5, each = 20L), rep(NA_real_, 10L)))
  seen <- integer(0)
  kern <- function(model, mf, subset, goes_left, ctrl) {
    seen <<- c(seen, length(subset))
    as.numeric(sum(goes_left))
  }
  sp <- argmax_split(
    splitter = kern,
    collector = new_collector(),
    model = list(object = NULL),
    mf = mf,
    subset = seq_len(nrow(mf)),
    whichvar = 1L,
    ctrl = list(minbucket = 20L)
  )
  expect_s3_class(sp, "partysplit")
  # The kernel must never be handed a row whose covariate is missing.
  expect_identical(unique(seen), 100L)
})

test_that("validate_tree_input() refuses covariate types partykit cannot scan", {
  ind <- parseModel(args_trees$.model)$indicators
  d <- dat
  d$flag <- dat$z_true == levels(dat$z_true)[1L]
  d$label <- as.character(dat$z_true)

  # partykit dies on both deep inside libcoin -- "object 'X' not found" for a
  # logical, "cannot handle objects of class 'character'" for a string -- and a
  # non-native kernel silently never splits them.
  expect_error(
    validate_tree_input(d, ind, c("z_true", "flag"), "mat", "FIT", args_trees),
    "flag"
  )
  expect_error(
    validate_tree_input(d, ind, c("z_true", "label"), "mat", "FIT", args_trees),
    "label"
  )
  expect_true(validate_tree_input(d, ind, "z_true", "mat", "FIT", args_trees))
})

# partykit:::.split() branches argmax_split() has to stand in for -------
# doTrees() replaces ctrl$splitfun outright, so .split() never runs and the two
# things it does around the cutpoint search are the splitter's to do instead.
# igsca_tree_control() exposes neither setting, so these are reachable only
# through ctrl -- which is exactly how partykit hands them over.

test_that("argmax_split() splits an unordered factor multiway, without the kernel", {
  # Two levels below minbucket, which partykit lumps into one extra kid
  # rather than dropping: index becomes c(1, 2, 5, 5) before it is re-coded.
  f <- factor(rep(c("a", "b", "c", "d"), times = c(40L, 40L, 3L, 3L)))
  mf <- data.frame(f = f)
  called <- FALSE
  kern <- function(model, mf, subset, goes_left, ctrl) {
    called <<- TRUE
    1
  }
  sp <- argmax_split(
    splitter = kern,
    collector = new_collector(),
    model = list(object = NULL),
    trafo = NULL,
    mf = mf,
    subset = seq_len(nrow(mf)),
    whichvar = 1L,
    ctrl = list(minbucket = 10L, multiway = TRUE, maxsurrogate = 0L)
  )
  # multiway is structural: partykit picks the kids without a statistic, so no
  # cutpoint rule is consulted and none of the kernels can override it.
  expect_false(called)
  expect_identical(partykit::index_split(sp), c(1L, 2L, 3L, 3L))
})

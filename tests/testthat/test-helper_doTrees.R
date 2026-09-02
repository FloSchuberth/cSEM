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

# One multigroup IGSCA fit, shared by the bdiagFit() and ndt_dists() tests
# below. Model syntax and csem() arguments as in tests/testthat/test-csem_fit.R.
# The grouping variable is deliberately not coded 1/2: labelling blocks by group
# name and by block index coincide when it is, which is why the two branches
# disagreeing went unnoticed for so long.
model_IGSCA <- "
OrgPres =~ cei1 + cei2 + cei3
OrgIden <~ ma1 + ma2 + ma3
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

dat_mg <- BergamiBagozzi2000
dat_mg$grp <- factor(paste0("g", dat_mg$gender))

res_mg <- csem(
  .data = dat_mg,
  .model = model_IGSCA,
  .approach_weights = "GSCA",
  .id = "grp",
  .tolerance = 1e-4,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP",
  .iter_max = 600
)

test_that("bdiagFit() block-diagonalizes multigroup fits in group order", {
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

test_that("bdiagFit() aligns a pooled matrix with an MGA through .block_names", {
  # One group of res_mg stands in for a pooled fit: same indicators, and it
  # takes the replication branch, which is the branch a pooled object takes.
  pooled <- res_mg[[1]]
  mga    <- bdiagFit(res_mg, .type_vcv = "indicator")

  # The defaults disagree by construction -- group names against block indices.
  # That is the disagreement .block_names exists to resolve.
  expect_false(identical(
    dimnames(bdiagFit(pooled, .n_blocks = 2L, .type_vcv = "indicator")),
    dimnames(mga)
  ))

  aligned <- bdiagFit(pooled, .n_blocks = 2L, .type_vcv = "indicator",
                      .block_names = names(res_mg))
  expect_identical(dimnames(aligned), dimnames(mga))
  # Only the labels moved: the blocks are the same replication as before.
  expect_equal(
    unname(aligned),
    kronecker(diag(2), unname(fit(pooled, .saturated = FALSE,
                                  .type_vcv = "indicator")))
  )
})

test_that("bdiagFit() validates .block_names", {
  expect_error(bdiagFit(res_pls, .n_blocks = 2, .block_names = "one"))
  expect_error(bdiagFit(res_pls, .n_blocks = 2, .block_names = c("a", "a")))
  expect_error(bdiagFit(res_pls, .n_blocks = 2, .block_names = c("a", NA)))
})

# doTrees node-level helpers ---------------------------------------------
# These call the partition machinery directly, one node at a time.
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

# At the tree level a broken kernel and an unsplittable node look identical (see
# ?warn_dead_splitter). Called directly, a kernel either returns a finite
# statistic or names its failure.
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

test_that("ndt_dists() refuses a pooled matrix the MGA does not conform to", {
  # calculateDL() indexes positionally and returns a number for a wrong
  # .n_blocks, a mismatched indicator set or a flipped group order alike, so
  # matching dimnames is the only check there is.
  pooled  <- bdiagFit(res_mg[[1]], .n_blocks = 2L, .type_vcv = "indicator")
  aligned <- bdiagFit(res_mg[[1]], .n_blocks = 2L, .type_vcv = "indicator",
                      .block_names = names(res_mg))

  expect_error(
    ndt_dists(Sc_pool = NULL, Si_pool = pooled, mga_fit = res_mg, dists = "DLi"),
    "dimnames"
  )
  expect_true(is.finite(
    ndt_dists(Sc_pool = NULL, Si_pool = aligned, mga_fit = res_mg,
              dists = "DLi")[["DLi"]]
  ))
})

test_that("a complex geodesic distance stumps the scan instead of killing the tree", {
  # calculateDG() can return a complex scalar with no error (see ndt_dists()).
  # Unguarded it coerced the whole distance vector to complex, and
  # vapply(..., numeric(1)) in argmax_split() rejects that AFTER the
  # per-candidate tryCatch has returned -- so the type error escaped into
  # partykit rather than degrading to "no admissible split".
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
  j <- match("noise_1", names(fx$mf))
  z <- fx$mf[[j]]
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

brk <- function(cands) {
  vapply(cands, function(cc) partykit::breaks_split(cc$split), numeric(1))
}

test_that("candidate_partitions() breaks where ctree's intersplit says to", {
  z <- c(1, 2, 4, 8)
  # ctree's default is intersplit = FALSE, which reports the observed value
  # below the gap rather than the midpoint.
  obs <- candidate_partitions(1L, z, z, minbucket = 1L)
  mid <- candidate_partitions(1L, z, z, minbucket = 1L, intersplit = TRUE)
  expect_identical(brk(obs), c(1, 2, 4))
  expect_identical(brk(mid), c(1.5, 3, 6))
  # Only the reported break moves. Every kernel scores `goes_left`, so if
  # intersplit reached that too it would change which partition wins, and not
  # merely how the winner is written down.
  expect_identical(
    lapply(obs, `[[`, "goes_left"),
    lapply(mid, `[[`, "goes_left")
  )
})

test_that("candidate_partitions() reports an observed break when no midpoint fits", {
  # Two adjacent doubles have no double between them, so their exact midpoint
  # is a tie, and IEEE rounds a tie to the even mantissa. For the pair below
  # that is the larger of the two: a break lands on b, where right = TRUE makes
  # the interval (-Inf, b] and kidids_split() routes b left -- while goes_left,
  # keyed on the observed value, sends it right. partykit's own interpolation
  # has no guard for this.
  #
  # The parity matters, which is why this pair and not the obvious one: at
  # c(1, 1 + double.eps) the tie rounds the other way, onto the smaller value,
  # where the fallback would have put the break anyway.
  a <- 1 - .Machine$double.eps / 2 # the predecessor of 1
  skip_if_not((a + 1) / 2 == 1, "platform rounds the tie away from the larger value")
  z <- c(a, 1, 3)
  cands <- candidate_partitions(1L, z, z, minbucket = 1L, intersplit = TRUE)
  # First cut falls back to a -- unguarded it would report 1. The second gap is
  # wide enough to hold a midpoint, and gets one.
  expect_identical(brk(cands), c(a, 2))
  for (cc in cands) {
    expect_identical(
      partykit::kidids_split(cc$split, data.frame(v = z)) == 1L,
      cc$goes_left
    )
  }
})

test_that("argmax_split() takes intersplit off the control object", {
  # doTrees() sets it once, on the ctree_control() object partykit hands back to
  # splitfun as `ctrl`. Reading it anywhere else would leave the native route
  # honouring the setting and the refitting ones ignoring it.
  mf <- data.frame(num = c(1, 2, 4, 8, 16, 32))
  args <- list(
    # Maximized by the largest admissible left kid, so both runs settle on the
    # same candidate and the reported break is the only thing free to differ.
    splitter = function(model, mf, subset, goes_left, ctrl) {
      as.numeric(sum(goes_left))
    },
    collector = new_collector(),
    model = list(object = NULL),
    trafo = NULL,
    mf = mf,
    subset = seq_len(nrow(mf)),
    whichvar = 1L
  )
  sp <- function(ctrl) {
    partykit::breaks_split(do.call(argmax_split, c(args, list(ctrl = ctrl))))
  }
  # Absent, as in the bare ctrl lists the tests above build: isTRUE() reads
  # that as the default rather than erroring on NULL.
  expect_identical(sp(list(minbucket = 1L)), 16)
  expect_identical(sp(list(minbucket = 1L, intersplit = FALSE)), 16)
  expect_identical(sp(list(minbucket = 1L, intersplit = TRUE)), 24)
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
  # Both intersplit settings, because moving the break is exactly the way to
  # break the agreement this test is about.
  for (isp in c(FALSE, TRUE)) {
    for (mb in c(1L, 2L)) {
      for (j in seq_along(mf)) {
        z <- mf[[j]]
        cands <- candidate_partitions(j, z, z, minbucket = mb, intersplit = isp)
        expect_gt(length(cands), 0L)
        for (cc in cands) {
          # The kernel is told `goes_left`; partykit later routes rows through
          # kidids_split(). A candidate whose two disagree scores one partition
          # and grows another.
          expect_identical(
            partykit::kidids_split(cc$split, mf) == 1L,
            cc$goes_left,
            info = paste0(
              names(mf)[j], ", minbucket = ", mb, ", intersplit = ", isp
            )
          )
        }
      }
    }
  }
  # keep_min() bit: at minbucket = 2 the two extreme cuts are inadmissible.
  # It filters on goes_left, which intersplit leaves alone, so the admissible
  # set is the same under both.
  expect_length(candidate_partitions(1L, mf$num, mf$num, minbucket = 1L), 5L)
  expect_length(candidate_partitions(1L, mf$num, mf$num, minbucket = 2L), 3L)
  expect_length(
    candidate_partitions(1L, mf$num, mf$num, minbucket = 2L, intersplit = TRUE),
    3L
  )
  # Ordered and unordered factors are untouched by it, as in partykit, whose
  # own interpolation is guarded with !is.ordered(x).
  for (v in c("ord", "fac")) {
    j <- match(v, names(mf))
    expect_identical(
      candidate_partitions(j, mf[[j]], mf[[j]], minbucket = 1L, intersplit = TRUE),
      candidate_partitions(j, mf[[j]], mf[[j]], minbucket = 1L),
      info = v
    )
  }

  # ctree sets index = 1:2 on every break-valued split it emits.
  for (v in c("num", "ord")) {
    j <- match(v, names(mf))
    sp <- candidate_partitions(j, mf[[j]], mf[[j]], minbucket = 1L)[[1L]]$split
    expect_identical(partykit::index_split(sp), 1:2, info = v)
  }
})

test_that("candidate_partitions() scans the numeric grid exhaustively", {
  z <- as.numeric(1:200)
  expect_length(candidate_partitions(1L, z, z, minbucket = 1L), 199L)
  expect_length(candidate_partitions(1L, z, z, minbucket = 20L), 161L)
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

test_that("validate_tree_input() refuses an unaffordable factor by name, up front", {
  # candidate_partitions() also refuses it, but from inside partykit's splitfun:
  # after arbitrary tree growth, and naming a model-frame column index rather
  # than the covariate. The native scan does not refit, so it is exempt.
  ind <- parseModel(args_trees$.model)$indicators
  d <- dat
  d$many <- factor(rep(letters[1:12], length.out = nrow(dat)))
  expect_error(
    validate_tree_input(d, ind, c("z_true", "many"), "mat", "DLi", args_trees),
    "many"
  )
  expect_true(
    validate_tree_input(d, ind, c("z_true", "many"), "mat", "native", args_trees)
  )
})

test_that("validate_tree_input() refuses a covariate with missing values", {
  # With maxsurrogate = 0 partykit routes a missing row by an unrepeated
  # sample() draw, once while growing the node and again when recording
  # tree$fitted, so the two partitions of one tree need not agree. Refused up
  # front rather than routed; see the note at the check.
  ind <- parseModel(args_trees$.model)$indicators
  d <- dat
  d$z_na <- dat$z_true
  d$z_na[1L] <- NA

  expect_error(
    validate_tree_input(d, ind, c("z_true", "z_na"), "mat", "DLi", args_trees),
    "z_na"
  )
  # The native scan routes them the same way, so it is not exempt.
  expect_error(
    validate_tree_input(d, ind, "z_na", "mat", "native", args_trees),
    "z_na"
  )
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

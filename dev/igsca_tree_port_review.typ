#set document(
  title: "doTrees() Port Review",
  author: "Claude Opus 5",
)
#set page(numbering: "1", margin: 2.4cm)
#set text(font: "New Computer Modern", size: 10.5pt)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")
#show raw.where(block: true): set block(
  fill: luma(245), inset: 8pt, radius: 3pt, width: 100%,
)

#let ok = text(fill: rgb("#1a7f37"), weight: "bold")[resolved]
#let open = text(fill: rgb("#b35900"), weight: "bold")[open]
#let part = text(fill: rgb("#8250df"), weight: "bold")[changed]

#align(center)[
  #text(size: 18pt, weight: "bold")[Porting `igsca_tree` / `igsca_ctree` into `doTrees()`]\
  #text(size: 10pt, style: "italic")[
    What is missing, what is wrong, and what to change #sym.dash.em branch `igscaTrees`
  ]
]

= Scope and method

This reviews the port of `dev/igsca/tree_template/igsca_tree.R` and
`igsca_testers.R` into `R/postestimate_doTrees.R` + `R/helper_doTrees.R`, and the
new `tests/testthat/test-postestimate_doTrees.R`.

Every finding below was *executed*, not read off. All five `influence` values were
run against `tests/testthat/data/igscaTrees.Rdata` (n = 1000, covariates `z_true`
(2-level factor), `noise_1`, `noise_2`) under `devtools::load_all()`, with reduced
`R_test`/`maxdepth`/`max_cuts` where the default control was too slow to finish.
Verbatim error text is quoted.

*Revision history.* The first pass found the partition family unrunnable. Since
then:

- *Second pass* --- `idx_permutation()` and `ndt_dists()` were added
  (`R/helper_doTrees.R:928`, `:959`) and `devtools::document()` was run. Verified
  in #link(<b1>)[B1] and #link(<p1>)[P1].
- *Third pass* --- `splitter` became a `match.arg()` string and the test file was
  updated to match. Verified in #link(<b2>)[B2] across all 17 configurations the
  tests exercise, plus both guard paths.
- *Fourth pass* --- the influence function moved to its own `influence_fn` local
  (#link(<b2>)[B2]), and S1 was fixed: `argmax_split()` now counts failed kernel
  scans, `doTrees()` warns on a wholly dead kernel, and the three tests of
  #link(<negctl>)[section 5.4] are implemented and passing.
- *Fifth pass* --- S2's dropped `cc$bonferroni` line was restored and verified on
  both distributions. #link(<testtype>)[Section 3.2.1] now answers the `testtype`
  TODO properly, including why the asymptotic route is genuine COIN inference and
  why partykit's "Bonferroni" is really Šidák.
- *Sixth pass* --- S3 fixed: leaves are refit after growing and attached to the
  terminal nodes, so `coef()` and `plot()` work and the tree finally carries a
  per-leaf IGSCA fit. Inner-node objects are kept, with an opt-in switch.
- *Seventh pass* --- T4 first, then T1--T3, which is the only workable order:
  until the controls are split by branch the suite cannot finish, and T2's node
  counts fail for reasons that have nothing to do with correctness. The file now
  runs in *2 min 05 s* at `FAIL 0 | WARN 0 | SKIP 0 | PASS 85`, against an
  estimated #sym.tilde 16 h before.

*All behavioural findings are closed, and the suite now runs and asserts.* Every
`influence` #sym.times `splitter` combination runs, the three silent failures are
fixed, and each of the 17 configurations is now asserted to have *actually
split* rather than merely not thrown. What remains is the unit layer that would
*localise* a regression --- T5 and #link(<mixed>)[section 5]'s 5.1--5.3 and 5.5
--- plus packaging cleanup (P1--P6).

== Status at a glance

#table(
  columns: (auto, 1fr, auto),
  inset: 6pt,
  align: (left, left, left),
  table.header([*ID*], [*Finding*], [*Status*]),
  [B1], [`idx_permutation()` / `ndt_dists()` not ported], [#ok],
  [B2], [`splitter` is a function in one branch, a string in the other], [#ok],
  [S1], [Broken split kernel degrades to a stump with no diagnostic], [#ok],
  [S2], [`control$bonferroni` silently ignored on the ctree path], [#ok],
  [S3], [`coef()` / `plot()` empty for any tree that splits], [#ok],
  [T1], [Snapshots stale (`trees_out` vs `trees_mx`)], [#ok],
  [T2], [`length(tree) == 5` is a disguised structural assertion], [#ok],
  [T3], [`expect_no_error()` cannot see the real failure mode], [#ok],
  [T4], [Default control makes the partition tests unrunnable], [#ok],
  [T5], [Collector diagnostics and the scan cache are untested], [#part],
  [P1], [`NAMESPACE` stale], [#part],
  [P2], [`cSEM:::` / `cSEM::` into the package's own namespace], [#open],
  [P3], [Roxygen placeholders; wrong `@return`], [#open],
  [P4], [Argument naming departs from package convention], [#open],
  [P5], [`fit_csem()` hard-codes the estimator], [#open],
  [P6], [No input validation], [#open],
)

The test file now carries two controls chosen per branch, five regenerated
snapshots, and one named `test_that()` per `influence` #sym.times `splitter`
configuration, each asserting the tree grew rather than merely returned. The
remaining test work is a *layer*, not a repair: nothing below the tree level is
covered, so a failure still says "this configuration broke" without saying
where. #link(<mixed>)[Section 5] is that layer, and 5.4 is already shipped.

#pagebreak()

= Resolved

== B1. `idx_permutation()` and `ndt_dists()` are now available <b1>

*Location:* `R/helper_doTrees.R:928`, `R/helper_doTrees.R:959`

Both helpers now exist, and both satisfy the contracts their call sites require.
This was the sole blocker on the partition family; it is cleared.

*`idx_permutation()` --- orientation is correct.* The returned matrix is
`R` #sym.times `n` (`nrow = R, ncol = n`), which is what `permutation_pvalue()`
needs when it indexes `goes_left[perm[r, ]]`. Had it come back transposed, the
result would have been a recycled-length nonsense null distribution rather than
an error, so this is worth having confirmed explicitly. The
sample-without-replacement-within-stratum construction matches
`boot:::permutation.array` as the docstring claims.

*`ndt_dists()` --- names and argument order match the call site.* It returns
`c(DGc, DGi, DLc, DLi)`, and `partition_stat()` reads `unname(ds[stat_kind])` for
`stat_kind` in `"DGi"` / `"DLi"`. The `(Sc_pool, Si_pool, mga_fit)` signature
lines up with `ndt_dists(coll$ndt_pools$Sc, coll$ndt_pools$Si, mga$fit)` at
`R/helper_doTrees.R:635`.

*Verified end to end.* The partition family, previously unrunnable, now grows
trees. At `R_test = 24L, max_cuts = 3L, maxdepth = 1L, minbucket = 200L` (14.8 s):

```
Fitted party:
[1] root
|   [2] z_true in 1: *
|   [3] z_true in 2: *

               z_true      noise_1      noise_2
statistic  0.04315312  0.001707126  0.001593087
p.value    0.04000000  0.280000000  0.083333333
criterion -0.04082199 -0.328504067 -0.087011377
```

The DLi/DGi kernels agree exactly with the reference implementations used to
pre-verify the path in the first pass (`argmax` cutpoint on `noise_1` = `-1.69654`
for both, `-0.231765` for FIT), so the shipped `ndt_dists()` computes what that
analysis assumed.

Two loose ends, both minor:

- *The `strata` argument is dead.* `idx_permutation(n, R, strata = NULL)`
  implements stratified permutation, but `permutation_pvalue()`
  (`R/helper_doTrees.R:699`) calls `idx_permutation(n, R)` with no `strata`. If
  stratified permutation within the node was intended --- which would be the
  natural thing to want when the node is already conditioned on a covariate ---
  it is not wired up. If it was not intended, drop the parameter rather than
  leave an untested branch in the package.
- *The provenance comment is now wrong.* `R/helper_doTrees.R:535` still says
  these come from `R/MGA/csem_test_helpers.R`. They are defined 400 lines below
  it. (`dev/igsca/tree_template/csem_test_helpers.R` is now present but
  untracked.)

== B2. `splitter` is now a single `match.arg()` string <b2>

*Location:* `R/postestimate_doTrees.R:16--29`

Implemented as recommended --- one string-valued argument with an explicit
`"native"` level, resolved once before either branch:

```r
doTrees <- function(data, model, covariates,
                    influence = c("mat", "vec", "FIT", "DLi", "DGi"),
                    splitter  = c("native", "FIT", "DLi", "DGi"),
                    control   = igsca_tree_control()) {
  influence <- match.arg(influence)
  splitter  <- match.arg(splitter)

  split_fn <- switch(splitter,
    native = NULL, FIT = split_max_fitdiff,
    DLi = split_max_dli, DGi = split_max_dgi)
```

`split_fn` is consumed by the ctree branch at `:117` and by `grow_extree()` at
`:221`. Both original failure modes are gone, and the structural hazard behind
the second one is closed: `match.arg()` guarantees a length-1 character, so
`splitter %in% c(...)` can no longer evaluate to `logical(0)` and slip through
`stopifnot()` vacuously.

The retained `stopifnot()` at `:167` is *not* redundant --- it now does real work,
rejecting `"native"` on the partition path, where there is no engine-native scan
and a `NULL` kernel would die inside `argmax_split()`'s `tryCatch` and silently
yield a stump. Both guards produce readable errors:

```r
doTrees(dat, model, covs, influence = "FIT", splitter = "native")
#> Error: splitter should be any one of 'FIT', 'DLi' or 'DGi' when the influence
#>        (selector) function is one of 'FIT', 'DLi' or 'DGi'

doTrees(dat, model, covs, influence = "mat", splitter = "nope")
#> Error: 'arg' should be one of "native", "FIT", "DLi", "DGi"
```

*Verified sweep.* All 17 configurations the test file exercises now run, none
error, and every one produces a split:

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  align: (left, right, right, right),
  table.header([*Configuration*], [*width*], [*n_fail (node/res)*], [*time*]),
  [`mat` #sym.times `native`], [3], [0 / 0], [0.9 s],
  [`mat` #sym.times {`FIT`, `DLi`, `DGi`}], [3], [0 / 2], [1.7--2.0 s],
  [`vec` #sym.times `native`], [3], [0 / 0], [0.3 s],
  [`vec` #sym.times {`FIT`, `DLi`, `DGi`}], [3], [0 / 0], [1.5--1.8 s],
  [{`FIT`, `DLi`, `DGi`} #sym.times {`FIT`, `DLi`, `DGi`}], [2], [0 / 0], [14.0--15.3 s],
)

(ctree rows at `R_test = 50, maxdepth = 2, max_cuts = 8`; partition rows at
`R_test = 24, max_cuts = 3, maxdepth = 1, minbucket = 200`. The
`n_fail_resample = 2` on the `mat` mixed rows is ordinary non-convergence at
extreme cutpoints, not a fault.)

*Residual nits: both since fixed.* The obsolete "add a switch statement" TODO has
been deleted, and the mid-function rebinding of `influence` is gone --- `:45` now
assigns the influence function to its own local, `influence_fn`, read at `:79`:

```r
influence_fn <- switch(influence, mat = influence_mat, vec = influence_vec)
...
h <- influence_fn(E)
```

`influence` therefore stays a character for the whole call, so the `else if
(influence %in% c("FIT","DLi","DGi"))` at `:155` now tests what it appears to
test. (An intermediate revision named this local `selector`, which collided with
the genuine extree `selectfun` of the same name at `:158`/`:219` --- two different
roles under one identifier in one function body. `influence_fn` avoids that:
`selector` is now unambiguously the variable-selection function partykit means by
the term.)

Re-running the full sweep after the rename gave results *byte-identical* to the
pre-rename run --- every `width` and `n_fail` value unchanged, differing only in
wall-clock jitter. As an independent check that the `switch()` still maps to the
right function, `influence = "mat"` and `"vec"` produce different trees, and their
node-2 cutpoints (`-1.09783` on `noise_1`, `1.46613` on `noise_2`) reproduce the
recorded snapshot exactly. That is a real cross-check rather than a coincidence:
on the native path the cutpoint comes from ctree's deterministic maxstat scan, so
it is invariant to `R_test` and the seed.

#pagebreak()

= Silent failures --- all three now resolved

== S1. A broken split kernel degraded to a stump with no diagnostic <s1>

*Location:* `R/helper_doTrees.R:293--345` (`argmax_split`, `warn_dead_splitter`)

*Status: fixed.* `argmax_split()` wraps each kernel call in `tryCatch(..., error
= function(e) NA_real_)`, which is right for a genuinely non-identified candidate
partition but also swallowed *programming* errors: every candidate became `NA`,
`any(is.finite(stats))` was `FALSE`, the loop `next`ed through every covariate,
`argmax_split()` returned `NULL`, and partykit read that as "no admissible
split". A dead kernel and an unsplittable node were indistinguishable in the
returned object --- which is why the original `ndt_dists()` defect had to be found
by hand rather than by the suite.

The `tryCatch` is unchanged (it still has to absorb non-convergent candidates).
What changed is that the two cases are now counted apart. `argmax_split()` tracks
whether the kernel was *evaluated at all* on this invocation, and only a scan
that ran and yielded nothing is recorded as a failure:

```r
    stats <- vapply(cands, function(cc) {
      tryCatch(splitter(model, mf, subset, cc$goes_left, ctrl),
               error = function(e) NA_real_)
    }, numeric(1))
    ## The kernel was actually evaluated on this covariate, so the invocation
    ## counts as a scan regardless of what came back.
    scanned <- TRUE
    if (!any(is.finite(stats))) next
    collector$n_split_scan <- collector$n_split_scan + 1L
    return(cands[[which.max(stats)]]$split)
  }
  if (scanned) {
    collector$n_split_scan <- collector$n_split_scan + 1L
    collector$n_fail_split <- collector$n_fail_split + 1L
  }
  NULL
```

The distinction is load-bearing: a node with no admissible partition exits via
`if (!length(cands)) next` and never sets `scanned`, so it is *not* counted as a
kernel failure. Both counters were added to `new_collector()` and both are
surfaced on the returned tree:

```r
attr(tr, "igsca_info")$n_split_scan   # scans where the kernel actually ran
attr(tr, "igsca_info")$n_fail_split   # ... of which produced nothing usable
```

A totally dead kernel --- `n_fail_split == n_split_scan` --- now also warns, via
a new internal `warn_dead_splitter()` called from both branches of `doTrees()`
just before the attribute is attached:

```
Warning: The 'DLi' split kernel produced no usable statistic at any of the 2
node(s) where it was scanned, so no split point could be chosen from it and the
tree was grown as if no split were admissible. Check that the kernel runs
standalone via partition_stat().
```

The warning deliberately fires only on total failure. A *partial* failure (some
nodes fine, others not) is ordinary non-convergence and would be noise as a
warning --- it is visible in the counters instead, where it can be asserted on.
`splitter = "native"` short-circuits the check and leaves both counters at `0`,
since no kernel is involved.

Tests are in #link(<negctl>)[section 5.4] and now pass.

== S2. `control$bonferroni` was silently ignored on the ctree path <s2>

*Location:* `R/postestimate_doTrees.R:116`

*Status: fixed.* The dropped line has been restored:

```r
cc <- partykit::ctree_control(...)
cc$bonferroni <- isTRUE(control$bonferroni)
```

It was never a no-op. `ctree_control()` derives `bonferroni` from `testtype`, but
the port computes `testtype` from `coin_distribution` first, so with the default
`coin_distribution = "approximate"` it always sent `testtype = "MonteCarlo"` and
`cc$bonferroni` was always `FALSE` --- while on the *partition* path
`grow_extree()` copies the whole control list into `ectrl`, so the same argument
did take effect there. Verified fixed on both distributions
(`R_test = 200`, three covariates):

```
                  z_true     noise_1    noise_2
asymptotic  FALSE  5.31e-57  0.872895   0.810172
asymptotic  TRUE   1.59e-56  0.997947   0.993160
MonteCarlo  FALSE  0.000     0.890      0.845
MonteCarlo  TRUE   0.000     0.998669   0.996276
```

=== Answering the TODO: `testtype` vs `bonferroni` <testtype>

The earlier draft of this section said `testtype` picks "permutation vs
asymptotic". That is half the story and worth stating properly, because
`testtype` conflates *two* dimensions:

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  align: (left, left, left),
  table.header([*`testtype`*], [*p-value from*], [*Multiplicity*]),
  [`"Bonferroni"`], [asymptotic limiting distribution], [adjusted],
  [`"Univariate"`], [asymptotic limiting distribution], [none],
  [`"MonteCarlo"`], [`nresample` permutations], [none],
  [`"Teststatistic"`], [no p-value --- raw statistic is the criterion], [n/a],
)

*Are the p-values not always permutation-based?* No --- and this is the crux. In
the conditional-inference (COIN) framework "permutation test" describes the
*null distribution*, not the *computation*. Strasser and Weber derived the
conditional mean $mu$ and covariance $Sigma$ of the linear statistic $T$ under
that permutation null *in closed form*. Standardising with them gives
$c = (T - mu)^tack.b Sigma^(+) (T - mu)$, which is asymptotically $chi^2$ with
$"rank"(Sigma)$ degrees of freedom (`teststat = "quadratic"`). You can then read
the p-value off that $chi^2$ limit *without ever drawing a permutation*. That is
the asymptotic route, and it is partykit's default.

Two measurements confirm no resampling happens on that path --- same seed,
`influence = "mat"`, root criteria:

```
approximate, R_test = 50 :  z_true 0.00      noise_1 0.74      noise_2 0.82
approximate, R_test = 200:  z_true 0.000     noise_1 0.890     noise_2 0.845
asymptotic,  R_test = 50 :  z_true 5.31e-57  noise_1 0.872895  noise_2 0.810172
asymptotic,  R_test = 200:  z_true 5.31e-57  noise_1 0.872895  noise_2 0.810172
```

First, `5.31e-57` is unreachable for *any* Monte-Carlo test at `R = 50` or `200`.
Second, the asymptotic values are bit-identical across the two `R_test` settings
--- `nresample` is simply unused. So `coin_distribution = "asymptotic"` is a real
setting with the same conditional null and the same $mu$, $Sigma$; it only swaps
an empirical tail for a $chi^2$ tail. It is far cheaper, and it yields continuous
p-values instead of a grid; Monte Carlo is the safer choice in small nodes where
the asymptotics have not bitten.

*Why the fix is still needed.* `.extree_node()` applies the adjustment to
whatever p-values it was handed, independently of how they were produced:

```r
if (ctrl$bonferroni)
    sf$criteria["p.value", ] <- sf$criteria["p.value", ] *
        sum(!is.na(sf$criteria["p.value", ]))
```

So the two knobs *are* orthogonal at the point of use, even though `testtype`
bundles them at the point of configuration --- which is exactly why an explicit
assignment was required.

*Two details worth knowing.* That row holds $log(1 - p)$, not $p$, so multiplying
by $k$ gives $log((1-p)^k)$ --- the adjustment is *Šidák*, $1 - (1-p)^k$, not
Bonferroni's $k p$. It agrees with Bonferroni to first order (hence the exact
$times 3$ on `z_true` above) but not for large p: $1 - (1 - 0.872895)^3 =
0.997947$, matching the table exactly. Second, partykit *does* support asking for
both dimensions at once --- `ctree_control()` accepts a length-2 `testtype`:

```r
testtype <- match.arg(testtype, several.ok = TRUE)
if (length(testtype) > 1) {
  stopifnot(all(testtype %in% c("Bonferroni", "MonteCarlo")))
  ttesttype <- "MonteCarlo"
}
... bonferroni = "Bonferroni" %in% testtype
```

`ctree_control(testtype = c("Bonferroni", "MonteCarlo"))` is identical to the
manual override on every non-function component (verified; the function-valued
slots differ only by closure environment). So the mapping could drop the override
entirely and stay inside the documented API:

```r
testtype <- c(
  if (isTRUE(control$bonferroni))                 "Bonferroni",
  if (control$coin_distribution == "approximate") "MonteCarlo"
)
if (is.null(testtype)) testtype <- "Univariate"
```

which reaches all four combinations --- MC#sym.plus.minus adjustment, asymptotic
#sym.plus.minus adjustment --- with no post-hoc surgery on the control object.
Optional: the current code is correct as it stands.

== S3. `coef()` and `plot()` were empty for any tree that splits <s3>

*Location:* `R/helper_doTrees.R` (`attach_leaf_fits`, `drop_inner_node_objects`,
`coef.igsca_tree`, `plot.igsca_tree`); wired in at
`R/postestimate_doTrees.R:150--160` and `:230--240`

*Status: fixed.* Both methods read `nobs` / `objfun` out of `info_node()` at
*terminal* nodes, but partykit only calls the trafo at nodes it *attempts to
split*, so leaves arrived with `info = NULL`. `coef()` returned `NULL` for every
tree with at least one split (it only appeared to work on stumps, where the root
*is* the leaf) and `plot()` drew unlabelled boxes, because `sprintf("n = %s",
NULL)` yields `character(0)` rather than erroring.

The deeper problem was not the methods: *the tree carried no per-leaf IGSCA fit
at all*. A new `attach_leaf_fits()` refits the model once per leaf after growing
and writes the result into the terminal nodes' own `info`, so both methods now
work through the ordinary `info_node()` route rather than a side table.

Two properties were verified before choosing that design:

- `as.list()` / `as.partynode()` round-trips a grown tree *identically*, so
  rewriting only the leaf entries leaves the rest of the structure untouched.
- Terminal-node `info` is invisible to `print()` --- partykit's
  `formatinfo_node()` falls back to `"*"` for list-valued info --- and that holds
  even when the info carries a full `cSEMResults`. Attaching leaf fits therefore
  does not change printed output or disturb the snapshots.

Measured, `maxdepth = 2` on the ctree path and a `FIT` #sym.times `FIT`
partition tree:

```
ctree      coef()          nobs   objfun        partition  coef()   nobs   objfun
             node 3          68       NA                     node 2   500  3553.173
             node 4         432  3076.803                    node 3   500  3438.003
             node 5         500  3438.003
```

Both previously returned `NULL`. `coef()` also gained row names (node ids) and
now emits an `NA` row instead of silently dropping a node --- `rbind()` discards
`NULL`s, which is precisely how the method used to collapse a whole tree to
`NULL`. `plot()`'s default panel prints `?` rather than an empty box when a value
is missing. `objfun` is the sum of squared casewise GSCA residuals, matching what
the ctree trafo stores on inner nodes; on the partition path it is the *only*
source of an SSR, since that trafo sets `objfun = NA_real_` throughout.

*What the fix exposes.* Leaf refits can fail, and now say so: a new
`n_fail_leaf` counter joins the `igsca_info` attribute, and the failing leaf gets
`converged = FALSE`, `objfun = NA`, `object = NULL`. The `maxdepth = 2` tree
above hit exactly that at node 3 (n = 68). This is worth knowing when choosing
`minbucket`: with 13 indicators the default `30L` admits leaves the model may not
fit. The failure is *sample-dependent rather than a size threshold* --- random
subsets of n = 30 and n = 68 fitted while n = 50 did not --- so it is a rising
probability of inadmissible estimates in small leaves, not a hard floor. Raising
`minbucket` for larger models is the cheap defence; `n_fail_leaf > 0` is the
signal that it is needed.

*Inner-node objects: kept, with an opt-in switch.* Every inner node's info still
holds a full `cSEMResults`, because `saveinfo = TRUE` persists whatever the trafo
returned. `drop_inner_node_objects()` removes those payloads after growing,
leaving criteria, `nobs`, `objfun` and all leaf fits intact. It is *not* called
by default; a commented-out line in each branch of `doTrees()` marks the switch:

```r
## MEMORY OPTIMIZATION (opt-in): each inner node's info holds a full
## cSEMResults object, so a maxdepth = 3 tree keeps up to 7 of them.
## Uncomment the next line to drop them; leaf fits are unaffected.
##   ret <- drop_inner_node_objects(ret)
## Do NOT instead delete `object = ft$fit` from the ytrafo above: a
## mixed-pair splitter reads model$object while the tree is still growing,
## so removing it there silently disables every non-native splitter.
```

That last warning is the point of putting the switch here rather than at the
trafo. Deleting `object = ft$fit` looks like the obvious optimisation and is the
one place it must not be done: on the ctree path `partition_stat()` reads
`model$object` for every candidate partition of a mixed-pair splitter, and on the
partition path the selectors read it too --- removing it there would take out
every non-native configuration, silently, exactly the way S1 used to fail.
Measured saving from the safe switch: 1.9 #sym.arrow.r 1.1 Mb (ctree,
2 inner nodes) and 1.2 #sym.arrow.r 0.7 Mb (partition, 1 inner node), with
`coef()` unchanged afterwards.

#pagebreak()

= Test-suite findings

The test file is the right shape --- one file per R file, matching the
conventions in `MEMORY.md` --- and T1--T4 are now fixed. This section records
what changed and what was measured along the way; T5 is what is left.

== T1. Snapshots were stale --- deleted and regenerated

*Status: fixed.* `_snaps/postestimate_doTrees.md` recorded the printed object
under the name `trees_out` while the tests bind `trees_mx` / `trees_vec` /
`trees_NPT_FIT` / ...; `expect_snapshot()` captures the *expression* as well as
the output, so both existing blocks failed on the `Code` line alone and the three
partition blocks had no snapshot at all. The stale file was removed (commit
`5e268c95`) and regenerated under the new controls --- five snapshots, one per
`influence` value.

Two properties were checked rather than assumed:

- *Reproducible across processes.* A second `devtools::test()` run in a fresh
  session gives `WARN 0` --- no "Adding new snapshot", no diff --- so the
  `set.seed()` at the head of each block pins both the libcoin Monte-Carlo draw
  and `idx_permutation()`.
- *The reference tree is unchanged by T4.* The `mat` snapshot still splits at
  `noise_1 <= -1.09783`, the cutpoint recorded in #link(<b2>)[B2] long before any
  of this, because the two `native` snapshot trees were deliberately left at the
  package defaults. The three partition snapshots are new and are all the same
  2-node split on `z_true`.

One packaging detail: removing the stale file took the last tracked file out of
`tests/testthat/_snaps/`, so the regenerated file is *untracked* and needs a
`git add tests/testthat/_snaps/` before it protects anything.

== T2. `length(tree) == 5` --- replaced with assertions that say what they mean

*Status: fixed.* `length()` on a `party` dispatches to `length.party()`, which
returns the *number of nodes*, not the length of a list, so the five
`expect_true(length(x) == 5)` lines asserted "this tree has exactly 5 nodes" ---
a strictly weaker duplicate of the snapshot immediately below, breaking on any
legitimate change to the control defaults or the RNG stream, for a reason the
`# Dirty substitute for snapshot` comment actively misdirected you about. They
also had to go before T4 could land: the reduced controls produce 2- and 3-node
trees, so all five would have failed for reasons unrelated to correctness.

They are replaced by one helper, applied to all 17 configurations:

```r
expect_grew <- function(tree) {
  expect_s3_class(tree, "igsca_tree")
  info <- attr(tree, "igsca_info")
  expect_identical(info$n_fail_full,  0L)
  expect_identical(info$n_fail_split, 0L)
  # width() returns a double, so expect_identical(..., 1L) would fail on type.
  expect_gt(partykit::width(partykit::node_party(tree)), 1)
}
```

Nothing now asserts a node count, which is what makes the per-branch controls
safe to tune: `expect_grew()` holds at any `maxdepth`, and the exact shape of
each tree is the snapshot's job.

== T3. `expect_no_error()` --- split into one named test per configuration

*Status: fixed.* Each of the five `expect_no_error({...})` blocks wrapped two or
three `doTrees()` calls, so a failure told you that *one* of them errored but not
which --- and, far worse, per S1 the realistic failure is not an error at all: a
broken kernel returns a plausible stump and the block passes. That is exactly
what happened in the first pass, when the DLi and DGi mixed splitters were
entirely non-functional and the suite reported green.

Each block is now a loop that *generates* one `test_that()` per splitter, so the
test name carries the configuration even when the call errors outright:

```r
for (sp in c("FIT", "DLi", "DGi")) {
  test_that(paste0("Matrix of Residuals selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("mat", sp, ctl_mixed))
  })
}
```

`grow_tree(influence, splitter, control)` is a one-line local wrapper around
`doTrees()`; it collapses 17 seven-line call sites to 17 one-line ones without
touching what `expect_snapshot()` captures.

*These assertions have teeth, verified.* #link(<negctl>)[Section 5.4]'s mocked
dead kernel makes exactly the `expect_grew()` assertions fail --- `width == 1`
and `n_fail_split > 0` --- which is the regression the old blocks could not see.

*Two counters were deliberately left out of `expect_grew()`*, both for measured
reasons rather than caution:

- `n_split_scan > 0` does *not* hold on the partition path. A matched
  selector/splitter pair reuses the selector's scan through the cache
  (#link(<cache>)[section 5.5]), so `FIT` #sym.times `FIT` legitimately reports
  `0 / 0`; only the mismatched pairs rescan.
- `n_fail_resample == 0L` does not hold either: the `mat` mixed rows reach 2, the
  ordinary non-convergence at extreme cutpoints noted in #link(<b2>)[B2].

== T4. Two controls, one per branch --- the suite finishes in two minutes

*Status: fixed.* The two paths differ in cost by an order of magnitude, because
only the partition *selectors* run permutation p-values; a mixed splitter
re-scans candidates but inherits COIN variable selection from libcoin. Measured
on this fixture:

#table(
  columns: (1fr, auto),
  inset: 6pt,
  align: (left, right),
  table.header([*Configuration*], [*Wall clock*]),
  [one 2-group MGA fit], [0.585 s],
  [ctree `native`, defaults (`maxdepth = 3, R_test = 500`)], [1.4 s],
  [mixed splitter, `maxdepth = 2, max_cuts = 8, R_test = 50`], [1.5--2.6 s],
  [partition selector, `maxdepth = 1, max_cuts = 3, R_test = 24`], [15.9 s],
  [partition selector, `maxdepth = 2, max_cuts = 4, R_test = 24`], [31.6 s],
  [partition selector, one node at `max_cuts = 20, R_test = 500`], [#sym.tilde 15 min],
  [partition selector, one tree at `maxdepth = 3` (default)], [#sym.tilde 1.8 h],
)

All 17 calls used a bare `igsca_tree_control()`. The eight ctree-branch calls
were fine at the defaults; the nine partition-branch calls were not --- at
#sym.tilde 1.8 h each the file was on the order of *16 hours*. What shipped:

```r
ctl_mixed <- igsca_tree_control(R_test = 50L, maxdepth = 2L, max_cuts = 8L)
ctl_part  <- igsca_tree_control(R_test = 24L, max_cuts = 3L, maxdepth = 1L,
                                minbucket = 200L)
```

`ctl_part` is *shallower and cheaper than this section originally recommended*,
on measurement. The proposed `maxdepth = 2L, max_cuts = 4L, minbucket = 100L`
runs in 31.6 s and returns the *same 3-node tree* as the `maxdepth = 1L` control
at 15.9 s: the children never split on this fixture, so the second level buys
runtime and nothing else. Depth-2 recursion is exercised on the ctree branch by
`ctl_mixed` instead, where it is nearly free.

`R_test` matters only for `ctl_part`, and is worth the comment it carries in the
file. `permutation_pvalue()` returns `(1 + k) / (R + 1)`, so the smallest
attainable p-value is `1/(R + 1)`, and the split criterion is a *strict*
inequality (`crit > logmincriterion` #sym.arrow.r.double `p < alpha`). At
`R_test = 19L` the floor is exactly `1/20 = 0.05`, which does *not* clear
`alpha = 0.05`: the partition family needs `R_test >= 20L` or it silently returns
stumps for reasons that have nothing to do with the data. `24L` clears it at
`1/25 = 0.04`.

This does *not* apply to `ctl_mixed`. libcoin's Monte-Carlo p-value is a plain
proportion with no add-one --- measured, `0.74 = 37/50` and `0.890 = 178/200`
exactly --- so `p = 0` is attainable at any `R_test` and there is no such floor on
the ctree branch.

*The two `native` snapshot trees keep the package defaults.* They cost 1.4 s
each, so there is nothing to buy by reducing them, and leaving them alone is what
lets the `mat` snapshot stay comparable to every earlier measurement in this
document (T1). `ctl_mixed` covers the six mixed ctree calls plus
#link(<negctl>)[section 5.4]'s three blocks; `ctl_part` covers all nine partition
calls.

*Result:* `FAIL 0 | WARN 0 | SKIP 0 | PASS 85` in *2 min 05 s*, so no part of
this file needs `skip_on_ci()`.

== T5. Nothing asserts the diagnostics the port exists to collect

Mostly addressed. #link(<negctl>)[Section 5.4]'s three blocks assert
`n_split_scan` and `n_fail_split` on all three paths (dead kernel, working
kernel, native), and T2's `expect_grew()` now asserts `n_fail_full` and
`n_fail_split` on all 17 configurations.

Still untouched: `n_fail_node`, `n_fail_resample`, `n_fail_leaf`,
`root_criteria`, and the node-local scan cache in `argmax_split()` --- the one
piece of genuinely subtle logic in the port (`identical()` on closures as a cache
key). See #link(<cache>)[section 5.5]. `n_fail_leaf` is the one worth doing next:
it is *live* on the default control (measured `n_fail_leaf = 1` on the `mat`
snapshot tree, the n = 68 leaf of #link(<s3>)[S3]), so it is a real diagnostic
firing in the suite today with nothing watching it.

#pagebreak()

= Testing the mixed-splitter configurations <mixed>

The mixed pairs (COIN variable selection + a model-comparison split point) are the
configurations the suite is least able to check, because their failure mode is a
*plausible-looking tree*, not an exception. T3 closed the crudest version of that
gap --- a wholly dead kernel now fails `expect_grew()` --- but a kernel that runs
and picks the *wrong* cutpoint still passes every assertion in the file. Three
measured facts shape what a useful test can look like.

*(a) The splitfun wiring works, but it is an undocumented contract.* `doTrees()`
plants its kernel in `cc$splitfun` and relies on partykit reading it back out of
the control object. A counting wrapper confirms it fires --- 9 kernel calls on a
`maxdepth = 2, max_cuts = 8` tree (1 at the root for `z_true`, 8 at node 2 for
`noise_1`) on partykit 1.2.29. The template's note says "verified against partykit
1.2.28", so this contract has already survived one upgrade unverified. It deserves
an explicit test.

*(b) The fixture cannot falsify cutpoint choice at the root.* `z_true` is a
2-level factor, so `candidate_partitions()` returns exactly *one* candidate:

```r
length(candidate_partitions(1, dat$z_true, dat$z_true, max_cuts = 20L, minbucket = 30L))
#> 1
length(candidate_partitions(2, dat$noise_1, dat$noise_1, max_cuts = 20L, minbucket = 30L))
#> 20
```

Every kernel is forced to return the same split at the root. Any root-level
assertion about *which* cutpoint a mixed splitter chose is vacuous, which is why
no amount of rephrasing the existing `expect_no_error()` would have worked.

*(c) Kernels do disagree, on numeric covariates.* Pinning `whichvar` to `noise_1`
and running `argmax_split()` directly:

```
FIT  cutpoint on noise_1 = -0.231765
DLi  cutpoint on noise_1 = -1.69654
DGi  cutpoint on noise_1 = -1.69654
```

So a differential assertion is meaningful --- but *`DLi` and `DGi` legitimately
agree*, both being model-implied indicator-VCV distances. Do not assert pairwise
distinctness across all three.

*A note on injection.* Now that B2 has made `splitter` a `match.arg()` string, a
test can no longer pass a custom kernel by value. Anything that needs to observe
or subvert a kernel must go through `local_mocked_bindings()`, which rebinds the
symbol in the package namespace where `doTrees()`'s `switch()` resolves it.
Verified working on testthat 3.3.2 (installed): the mocked counting wrapper
below records the same 9 calls the by-value version did.

== Kernel unit tests --- the layer that was missing

Two MGA fits each, and the failure message names the broken kernel. This is the
test that would have caught the missing `ndt_dists()` in seconds:

```r
node_fixture <- function() {
  ind <- parseModel(model)$indicators
  ft  <- try_fit(dat[, ind, drop = FALSE], model)
  list(model = list(object = ft$fit),      # what the trafo hands the kernel
       mf    = dat,
       ctrl  = list(collector = new_collector(), model = model,
                    indicators = ind, max_cuts = 20L, minbucket = 30L))
}

test_that("every split kernel returns a finite statistic", {
  fx <- node_fixture()
  gl <- fx$mf$noise_1 <= stats::median(fx$mf$noise_1)
  for (k in c("FITdiff", "DLi", "DGi")) {
    expect_true(
      is.finite(partition_stat(k, fx$model, fx$mf, seq_len(nrow(fx$mf)), gl, fx$ctrl)),
      info = k
    )
  }
})
```

`partition_stat()` converts a missing-function error into `NA` +
`n_fail_resample++`, so `is.finite()` catches it. Observed values at the median
split are `FITdiff = -0.00398`, `DLi = 0.0439`, `DGi = 0.0136` --- note the
negative FIT difference, which is exactly the case `partition_select()`'s
`.Machine$double.eps` floor exists to handle, and worth a regression test of its
own.

== `argmax_split()` returns an admissible cutpoint

```r
test_that("argmax_split returns an admissible cutpoint for every kernel", {
  fx <- node_fixture(); j <- match("noise_1", names(fx$mf))
  for (k in c("FIT", "DLi", "DGi")) {
    kern <- switch(k, FIT = split_max_fitdiff, DLi = split_max_dli, DGi = split_max_dgi)
    sp <- argmax_split(kern, fx$ctrl$collector, fx$model, fx$mf,
                       seq_len(nrow(fx$mf)), j, fx$ctrl)
    expect_s3_class(sp, "partysplit")
    br <- partykit::breaks_split(sp)
    expect_true(is.finite(br), info = k)
    expect_true(br > min(fx$mf$noise_1) && br < max(fx$mf$noise_1), info = k)
  }
})

test_that("different kernels place different cutpoints", {
  fx <- node_fixture(); j <- match("noise_1", names(fx$mf))
  cut <- function(kern) partykit::breaks_split(
    argmax_split(kern, fx$ctrl$collector, fx$model, fx$mf,
                 seq_len(nrow(fx$mf)), j, fx$ctrl))
  # NB: DLi and DGi legitimately agree -- both are indicator-VCV distances.
  expect_false(isTRUE(all.equal(cut(split_max_fitdiff), cut(split_max_dli))))
})
```

These call the kernels directly rather than through `doTrees()`, so they are
unaffected by B2 and need no mocking. This is where cutpoint choice *must* be
tested, per fact (b) --- at the tree level the root gives the kernel no freedom.

== Assert that partykit actually calls the splitfun

```r
test_that("doTrees installs its splitfun into partykit's split search", {
  calls <- new.env(parent = emptyenv()); calls$n <- 0L
  real  <- split_max_fitdiff
  local_mocked_bindings(
    split_max_fitdiff = function(model, mf, subset, goes_left, ctrl) {
      calls$n <- calls$n + 1L
      real(model, mf, subset, goes_left, ctrl)
    })
  set.seed(11)
  tr <- doTrees(dat, model, covs, influence = "mat",
                splitter = "FIT", control = ctl_mixed)
  expect_gt(calls$n, 0L)                                    # the contract
  expect_gt(partykit::width(partykit::node_party(tr)), 1L)  # it split
})
```

Verified: `calls$n == 9`, `width == 3`. Capturing `real <- split_max_fitdiff`
*before* the mock is what avoids infinite recursion. This test is what breaks
loudly if a partykit upgrade stops honouring `ctrl$splitfun`.

== Negative control --- give the suite teeth <negctl>

*Status: implemented and passing* (`tests/testthat/test-postestimate_doTrees.R`,
10 assertions across three blocks). One test alone would not have been enough:
asserting only that a dead kernel *warns* leaves open whether a working kernel
warns too, so the trio brackets the behaviour from both sides plus the
kernel-free path.

As shipped, after T2--T4 hoisted `ctl_mixed`, `grow_tree()` and `expect_grew()`
to the top of the file:

```r
test_that("a dead split kernel is reported, not silently a stump", {
  local_mocked_bindings(
    split_max_dli = function(model, mf, subset, goes_left, ctrl) {
      stop("kernel is broken")
    }
  )
  set.seed(11)
  expect_warning(
    tr <- grow_tree("mat", "DLi", ctl_mixed),
    "produced no usable statistic"
  )
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_fail_split, 0L)
  expect_identical(info$n_fail_split, info$n_split_scan)
  # The tree really is the stump the diagnostic is warning about
  expect_equal(partykit::width(partykit::node_party(tr)), 1)
})

test_that("a working split kernel records scans but no failures", {
  set.seed(11)
  tr <- expect_no_warning(grow_tree("mat", "DLi", ctl_mixed))
  expect_gt(attr(tr, "igsca_info")$n_split_scan, 0L)
  expect_grew(tr)
})

test_that("the native split path never touches the kernel counters", {
  set.seed(11)
  tr <- grow_tree("mat", "native", ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_split_scan, 0L)
  expect_identical(info$n_fail_split, 0L)
})
```

The middle block is the positive control for `expect_grew()` itself: it runs the
same configuration as the first, unmocked, so between them the pair pins that the
warning and the counters track the kernel rather than the control. Two details
worth carrying: the reduced `ctl_mixed` keeps all three blocks to a few seconds,
which costs nothing in fidelity since mixed splitters run no permutation test
(T4); and `partykit::width()` returns a *double*, not an integer, so
`expect_identical(width(...), 1L)` fails on type alone --- use
`expect_equal(..., 1)` or `expect_gt()`.

== The scan cache <cache>

`argmax_split()`'s cache key is `identical()` on the splitter *closure*
(`R/helper_doTrees.R:294--296`). A matched selector/splitter pair should reuse the
selector's scan; a mismatched pair should rescan with its own kernel. Both are
cheap to check at unit level with a pre-armed collector, and neither is covered
today:

```r
test_that("a matched pair reuses the scan and a mismatched pair rescans", {
  fx <- node_fixture(); j <- match("noise_1", names(fx$mf)); subset <- seq_len(nrow(fx$mf))
  coll <- fx$ctrl$collector
  sc <- scan_covariate("FITdiff", j, fx$model, fx$mf, subset, fx$ctrl)
  coll$scan <- stats::setNames(list(sc), as.character(j))
  coll$scan_subset   <- subset
  coll$scan_splitter <- split_max_fitdiff

  expect_identical(                                   # matched -> cache hit
    argmax_split(split_max_fitdiff, coll, fx$model, fx$mf, subset, j, fx$ctrl),
    sc$split)
  expect_false(isTRUE(all.equal(                      # mismatched -> rescan
    argmax_split(split_max_dli, coll, fx$model, fx$mf, subset, j, fx$ctrl),
    sc$split)))
})
```

Note this test must *not* be combined with `local_mocked_bindings()` on the same
kernels: the cache key is closure identity, and a mock would change it.

== What none of this proves

The five layers above establish that the mixed splitters *run, are invoked, and
discriminate*. They do not establish that a mixed splitter finds the *right*
cutpoint, because this fixture has no numeric covariate with a true threshold ---
`z_true` is binary and `noise_*` are noise. A recovery test needs a fixture where
a numeric covariate drives a known group difference at a known cutpoint, then
asserts the recovered breakpoint falls within a tolerance of it. That is the only
test that checks the mixed configurations are statistically correct rather than
merely alive, and it is worth building if the mixed pairs are a real study arm
rather than a convenience.

#pagebreak()

= Packaging and API

== P1. `document()` has been run --- now prune the export list <p1>

*Status: inverted.* `NAMESPACE` has been regenerated and `man/` populated, so the
first pass's complaint is resolved: `igsca_tree_control` and the `split_max_*`
kernels are now reachable, and `cSEM::doTrees(control = cSEM::igsca_tree_control(...))`
works.

The problem is now the opposite one. 26 new symbols are exported, and most are
internal machinery with no standalone meaning:

```
argmax_split  csem_converged  fit_csem  grow_extree  idx_permutation
ndt_dists  ndt_pools  node_group_data  partition_select  partition_stat
permutation_pvalue  scan_covariate  select_ndt_dgi  select_ndt_dli
select_npt  split_max_dgi  split_max_dli  split_max_fitdiff  try_fit
candidate_partitions  influence_mat  influence_vec  ...
```

Every one of these is now public API you are committing to keep stable, and each
carries a `man/*.Rd` that `R CMD check` will scrutinise (see P3). Replace
`@export` with `@keywords internal` or `@noRd` on all of them. B2's string-valued
`splitter` has already removed the last *user-facing* reason to expose the
`split_max_*` kernels --- but keep them (or their `partition_stat()` backend)
reachable from the test file, since the unit tests in
#link(<mixed>)[section 5] call them directly. A defensible public surface is
`doTrees`, `igsca_tree_control`, `bdiagFit`, `root_criteria`, and the two S3
methods.

== P2. `:::` and `::` into the package's own namespace

*Location:* `R/helper_doTrees.R:572--573` and `R/helper_doTrees.R:608`

Unchanged:

```r
Sc = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "construct"),
Si = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "indicator")
...
cSEM::calculateFIT(mga$fit) - cSEM::calculateFIT(model$object)
```

Both are defined in this package (`bdiagFit` at `R/helper_doTrees.R:26`,
`calculateFIT` at `R/helper_assess.R:2013`). `R CMD check` flags the `:::` form
("a package almost never needs to use `:::` for its own objects"), and both break
under `load_all()` if the installed and sourced versions diverge.

The new code makes the inconsistency starker rather than better: `ndt_dists()`
(`R/helper_doTrees.R:960--966`) calls `bdiagFit()`, `calculateDG()` and
`calculateDL()` bare and correctly, sitting 390 lines below `ndt_pools()` doing
the same thing through `cSEM:::`. Drop the prefixes.

== P3. Roxygen placeholders --- now materialised as `.Rd` files

`document()` has turned the placeholder blocks into committed documentation.
`man/idx_permutation.Rd` in full:

```
\title{Title}
\usage{idx_permutation(n, R, strata = NULL)}
\arguments{
\item{strata}{}
}
\description{
Permutation resampling: sample WITHOUT replacement, within each stratum, ...
}
```

Three `R CMD check` problems in nine lines: the title is literally `Title`, `n`
and `R` are undocumented while `strata` has an empty description, and there is no
`\value{}`. That pattern repeats across the 20-odd new `.Rd` files. P1's pruning
is the cheapest fix --- `@noRd` deletes the file and the warnings with it --- and
whatever survives as public API needs a real title, real `@param` text and a real
`@return`.

Two specific errors in `man/doTrees.Rd`:

- *`\value` is wrong.* It claims class `modelparty` and `party`. Actual class is
  `c("igsca_tree", "constparty", "party")` on both paths; there is no `modelparty`
  anywhere in the package.
- *Half the arguments are undocumented.* `data`, `model` and `covariates` did
  inherit from `csem_arguments` (roxygen matched them to the dot-prefixed
  `.data` / `.model` / `.covariates` entries). `influence`, `splitter` and
  `control` did not --- they appear nowhere in `\arguments{}`, which is an
  `R CMD check` WARNING. They are also the three arguments a reader most needs
  explained, since `influence` selects between two entirely different algorithms
  and `splitter` is constrained by which one you picked (B2's `stopifnot`).

== P4. Argument naming departs from the package convention

Every other exported cSEM entry point takes dot-prefixed arguments (`.data`,
`.model`, `.approach_weights`, ...). `doTrees()` takes bare `data`, `model`,
`covariates`, `influence`, `splitter`, `control`.

(The first pass claimed this also broke `@inheritParams`; that was wrong ---
roxygen matched the dotted entries to the undotted parameters, as P3 shows. The
case for renaming rests on API consistency alone, which is still a good enough
reason while the branch is young. Doing it together with B2's now-settled
signature would keep the churn to one commit.)

== P5. `fit_csem()` hard-codes the estimator with no override

*Location:* `R/helper_doTrees.R:152--164`

```r
fit_csem <- function(.data, .model, .id = NULL) {
  csem(.data = .data, .model = .model, .id = .id,
       .approach_weights = "GSCA", .disattenuate = TRUE,
       .conv_criterion = "sum_diff_absolute", .iter_max = 100,
       .GSCA_modes = "CCMP", .tolerance = 0.0001)
}
```

Reasonable as a study default, but `doTrees()` gives the user no way past it ---
no `...`, no argument. Since the tree machinery is estimator-agnostic (nothing
outside `partition_stat()`'s FIT branch is GSCA-specific), threading `...` from
`doTrees()` through `try_fit()` into `csem()` costs almost nothing and makes
`doTrees()` usable for PLS-PM trees too. At minimum, document that the estimator
is fixed --- especially now that `fit_csem` is exported.

== P6. No input validation

`doTrees()` validates `influence` and `splitter` (via `match.arg`) and nothing
else. Worth adding, since all three cases currently produce errors from deep
inside partykit or `csem()`:

```r
stopifnot(
  "`data` must be a data.frame" = is.data.frame(data),
  "`covariates` must be columns of `data`" =
    length(covariates) > 0L && all(covariates %in% names(data))
)
missing_ind <- setdiff(indicators, names(data))
if (length(missing_ind)) {
  stop2("The following indicators are not columns of `data`: ",
        paste(missing_ind, collapse = ", "))
}
```

(`stop2()` is the package's own error helper, `R/zz_utils.R:114`.)

#pagebreak()

= What the port got right

Worth recording, since these are deliberate improvements over the template rather
than accidents:

- *Indicators from the model, not the data.* `parseModel(model)$indicators`
  (`R/postestimate_doTrees.R:33`) replaces the template's
  `setdiff(names(data), covariates)`. The template silently treated any stray
  column as an indicator; the port does not, and now tolerates extra columns in
  `data`. The two orderings differ (`x31, x32, ..., x11` vs `x11, x12, ...`),
  which is why the recorded snapshot's model formula reads
  `~x31 + x32 + x33 + x11 + ...` --- expected, not a bug.
- *The unified dispatch is the right cut.* Branching on `influence` rather than
  exposing two constructors correctly reflects that the COIN and partition
  families differ only in how the node statistic is produced.
- *`idx_permutation()` and `ndt_dists()` were ported carefully.* The `R` #sym.times `n`
  orientation and the `(Sc_pool, Si_pool, mga_fit)` argument order both match
  their call sites, and `ndt_dists()` uses unprefixed internal calls --- the
  convention the older `ndt_pools()` still violates.
- *B2 was fixed properly rather than patched.* Making `splitter` a `match.arg()`
  string closed the `logical(0)` hole structurally instead of special-casing
  `NULL`, and kept the one `stopifnot()` that still carries meaning. Both guard
  paths now produce errors a user can act on.
- *`influence_ssr` #sym.arrow.r `influence_vec`* reads better next to
  `influence_mat`. The stale "COIN_ssr" references in the comments
  (`R/helper_doTrees.R:113`, `:514`) should follow the rename, as should the
  now-incorrect provenance comment at `:535`.
- The scan cache, minbucket filtering, and `.Machine$double.eps` flooring all
  survived the port intact and behave correctly under test.

= Suggested order of work

#table(
  columns: (auto, 1fr, auto),
  inset: 6pt,
  align: (left, left, left),
  table.header([*\#*], [*Item*], [*Why in this position*]),
  [#strike[1]], [#strike[T4 --- two controls, chosen by branch]], [*Done.* #sym.tilde 16 h #sym.arrow.r 2 min 05 s],
  [#strike[2]], [#strike[T1--T3 --- regenerate snapshots, drop `length()==5`, split `expect_no_error()`]], [*Done.* 85 assertions, `FAIL 0 | WARN 0`],
  [3], [§5.1--5.3 --- kernel, `argmax_split`, and splitfun-wiring tests], [*Now the head of the queue.* Cheap, fast, and they pin the contracts the port depends on],
  [4], [T5 + §5.5 --- scan-cache tests], [Subtlest logic in the port, and still wholly unasserted],
  [5], [P1--P6 --- prune exports, drop `:::`, fix `\value`, validate inputs], [Cleanup; P4's rename pairs naturally with B2's settled signature. Note `attach_leaf_fits()`, `drop_inner_node_objects()` and `warn_dead_splitter()` are already `@noRd`],
  [6], [#link(<s3>)[S3] follow-up --- assert `n_fail_leaf` and raise `minbucket` for large models], [The default `30L` admits unfittable leaves, and the counter fires in the suite today (`n_fail_leaf = 1` on the `mat` snapshot) with nothing watching it],
  [7], [§5.6 --- numeric-cutpoint fixture], [Optional; only if the mixed pairs are a real study arm],
  [8], [#link(<testtype>)[§3.2.1] --- optional: replace the `cc$bonferroni` override with a length-2 `testtype`], [Cosmetic; stays inside the documented partykit API],
)

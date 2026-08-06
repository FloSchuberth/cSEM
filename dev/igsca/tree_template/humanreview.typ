#set document(
  title: "doTrees() Port Human Review",
  author: "Michael",
)
// ===== Code blocks (codly + iridis paren colouring) =====
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.10": *
#import "@preview/iridis:0.2.0" // For colored parentheses
#show: codly-init.with()
  #codly(languages: codly-languages, stroke: 1pt + red)
  // iridis-show also installs a `show math.equation` rule that is broken on Typst
  // 0.14 (it reads a removed `accent.size` field, dying on accents like overline()).
  // Apply only its code-parenthesis colorizer, leaving math untouched.
  #show raw: iridis.internals.colorize-code(
    state("parenthesis", 0),
    ("(", "[", "{"), (")", "]", "}"),
    iridis.iridis-palette,
  )

#set page(numbering: "1", margin: 2.4cm)
#set text(font: "New Computer Modern", size: 10.5pt)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")


= On-Going Questions

+ Does it matter whether or not `coin_distribution = "asymptotic"` or `"approximate"`?
  - What is the default in ctree?
  - Do all the different options work as intended?
  - What about when I use the partition family path instead of conditional-test?
+ Why does the documentation seem catered to someone who already understands the codebase, as opposed to someone reading it for the first time?
+ Does igscaTree handle convergence failures in a similar way as semtree?
  - Within a resample?
  - Within a particular tested multigroup model?
  - Separation in convergence failure for variable selection versus split point selection
+ What is a kernel?
+ What helper functions or ways of interacting with the fitted tree need to be added/documented?
  - Split-point + p-value?
+ What should I make the default arguments for `igsca_tree_control()`?

= helper_doTrees.R

== igsca_tree_control

+ How much does this overlap or contradict with `partykit::ctree_control()`? Is the functionality consistent? The documentation?
+ Do the default arguments make sense?

= test-postestimate_doTrees.R

+ What are reasonable defaults for `ctl_mixed` and `ctl_part`?

= postestimate_doTrees.R

== `doTrees()`

+ What do each of the attributes of igsca_info mean?
+ What does the 'native' splitter using coin really mean?


=== Outline of function

Preliminary:

```r
split_fn <- switch(
    splitter,
    native = NULL,
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )
csem_tree_args(...) # helper_doTrees.R
validate_tree_input(...) # helper_doTrees.R
fml <- paste(...) |> stats::as.formula()
collector <- new_collector() # helper_doTrees.R
```

Functions to investigate in their order:
+ `csem_tree_args()`
+ `validate_tree_input()`
+ `new_collector()`

Switch between: (1) conditional inference or (2) multigroup approach

```r
if (influence %in% c("mat", "vec")) {
  # conditional inference
} else if (influence %in% c("FIT", "DLi", "DGi")) {
  # multigroup approach
}
```

==== Conditional Inference

```r
influence_fn <- switch(influence, mat = influence_mat, vec = influence_vec)
ytrafo <- function(data, weights, control) {
  ft <- try_fit(...)
  E <- calculateGSCAErrors(ft$fit)
  h <- influence_fn(E)
  ef <- matrix(0, ...) 
  ...
  return(list(
          estfun = ef,
          converged = TRUE,
          objfun = sum(E^2), # TODO: Should I make this 1 - calculateFIT(...), instead?
          object = ft$fit,
          nobs = length(subset)
        ))
}
```

Functions to investigate in their order:
+ `try_fit()`


Preparing the control options

```r
testtype <- c(
  if (isTRUE(control$bonferroni)) "Bonferroni",
  if (control$coin_distribution == "approximate") "MonteCarlo"
)
if (is.null(testtype)) {
  testtype <- "Univariate"
}
cc <- partykit::ctree_control(
  teststat = "quadratic",
  splitstat = "quadratic",
  testtype = testtype,
  nresample = control$R_test,
  alpha = control$alpha,
  minsplit = control$minsplit,
  minbucket = control$minbucket,
  maxdepth = control$maxdepth,
  maxsurrogate = 0L,
  nmax = c(yx = Inf, z = Inf),
  saveinfo = TRUE
)
```

Over-riding the default splitter if something other than `"native"` was passed to `doTrees(.splitter = ...)`, (i.e., `FIT`, `DLi`, `DGi`)

```r
if (!is.null(split_fn)) {
  # Sets the splitter to one of the three functions that we're looking for.
  cc$args <- args
  cc$indicators <- indicators
  cc$collector <- collector
  cc$max_cuts <- control$max_cuts
  cc$splitfun <- function(
    model,
    trafo, # TODO: Why is this trafo instead of ytrafo?
    data,
    subset,
    weights,
    whichvar,
    ctrl
  ) {
    argmax_split(
      split_fn,
      collector,
      model,
      model.frame(data),
      subset,
      whichvar,
      ctrl
    )
  }
  cc$svsplitfun <- cc$splitfun # never called (maxsurrogate = 0)
}
```

+ When is the `svsplitfun()` ever used? Should I force `maxsurrogate = 0` internally of `igsca_tree_control()`?
+ What is `argmax_split(...)` up to?

Launching IGSCA conditional tree

```r
ret <- partykit::ctree(fml, data = data, ytrafo = ytrafo, control = cc)
```

Return Formatting

```r
class(ret) <- c("igsca_tree", class(ret))
warn_dead_splitter(collector, splitter)

## partykit never calls the trafo at a node it does not try to split, so
## leaves arrive with info = NULL and carry no IGSCA fit. Refit them.
ret <- attach_leaf_fits(ret, ret$data, args, indicators, collector)

## MEMORY OPTIMIZATION (opt-in): each inner node's info holds a full
## cSEMResults object, so a maxdepth = 3 tree keeps up to 7 of them.
## Uncomment the next line to drop them; leaf fits are unaffected.
##   ret <- drop_inner_node_objects(ret)
## Do NOT instead delete `object = ft$fit` from the ytrafo above: a
## mixed-pair splitter reads model$object while the tree is still growing,
## so removing it there silently disables every non-native splitter.

attr(ret, "igsca_info") <- list(
  n_fail_full = collector$n_fail_full,
  n_fail_node = collector$n_fail_node,
  ## 0 on the native path (COIN resampling lives inside libcoin); counts
  ## failed candidate MGA fits when a mixed-pair splitter is in play.
  n_fail_resample = collector$n_fail_resample,
  ## Split-kernel scans: n_fail_split counts the scans that produced no
  ## finite statistic. Both stay 0 when splitter = "native".
  n_split_scan = collector$n_split_scan,
  n_fail_split = collector$n_fail_split,
  ## Failed IGSCA refits at terminal nodes (attach_leaf_fits).
  n_fail_leaf = collector$n_fail_leaf,
  root_criteria = root_criteria(ret)
)
return(ret)
```

+ How does `warn_dead_splitter()` work?
+ What does `attach_leaf_fits()` do?
+ What do each of the above collector fields mean?

===== `partykit::ctree()`

The entry is

```r
ret <- partykit::ctree(fml, data = data, ytrafo = ytrafo, control = cc)
```

What do the following unused arguments of `ctree()` do? Should I care about them?

+ `converged`
  - No, because we verify convergence internally and have our own methods for handling them.
  - May want to revisit this --- this is entirely due to the use of the `collector`
+ `scores`
+ `doFit`


Creation of data

```r
mf <- match.call(expand.dots = FALSE) # TODO: What does this do? Store the calling function of `partykit::ctree(...) into mf?
m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "scores"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$yx <- "none" # TODO: What does mf$yx do and why does partykit look in ytrafo for it?
if (is.function(ytrafo)) {
    if (all(c("y", "x") %in% names(formals(ytrafo)))) 
        mf$yx <- "matrix"
}
mf$nmax <- control$nmax
mf[[1L]] <- quote(partykit::extree_data) # Overwrite the argument to use `partykit::extree_data` instead
d <- eval(mf, parent.frame()) # Create data

subset <- .start_subset(d)
weights <- model.weights(model.frame(d))
```

+ `mf$yx` is "none" for `igsca_ctree`, does this matter?
+ What is `control$nmax`? Why does it have two entries of `yx` and `z`?
+ What is so important about `partykit::extree_data`?

Validating and manipulating `ytrafo`

```r
if (is.function(ytrafo)) {
    if (is.null(control$update)) # is TRUE
        control$update <- TRUE
    nf <- names(formals(ytrafo))
    if (all(c("data", "weights", "control") %in% nf)) # is TRUE
        ytrafo <- ytrafo(data = d, weights = weights, control = control) # TODO: What is this doing? Is it overwriting the function...?
    nf <- names(formals(ytrafo))
    stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) || all(c("y", "x", "weights", "offset", "start") %in% nf))
}
```

+ What does `control$update` have to do with anything? And does it matter that it's set to `TRUE` in my context?
+ What is happening on Line 6? Is `ytrafo` being overwritten to an object? A different function? Why?


Validating and manipulating `converged`

```r
# Not relevant because don't pass a converged function, but it's interesting what exactly is happening here...
if (is.function(converged)) {
    stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
    converged <- converged(d, weights, control = control)
}
else {
    converged <- TRUE
}
```

Create an `update` function for building the actual tree

```r
update <- function(subset, weights, control, doFit = TRUE) {
    extree_fit(
        data = d,
        trafo = ytrafo,
        converged = converged,
        partyvars = d$variables$z,
        subset = subset,
        weights = weights,
        ctrl = control,
        doFit = doFit
    )
}
tree <- update(subset = subset, weights = weights, control = control)
```

Output formatting

```r
trafo <- tree$trafo
tree <- tree$nodes
mf <- model.frame(d)
if (is.null(weights)) 
    weights <- rep(1, nrow(mf))
fitted <- data.frame(`(fitted)` = fitted_node(tree, mf), `(weights)` = weights, check.names = FALSE)
fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
names(fitted)[3] <- "(response)"
ret <- party(tree, data = mf, fitted = fitted, info = list(call = match.call(), control = control))
ret$update <- update
ret$trafo <- trafo
class(ret) <- c("constparty", class(ret))
ret$terms <- d$terms$all
return(ret)
```

+ What does `party(...)` do?

====== `extree_fit(...)`


====== `party(...)`

==== Multigroup Approach

```r

```
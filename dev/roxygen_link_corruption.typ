#set document(
  title: "Silent link repointing in roxygen2",
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

#align(center)[
  #text(size: 18pt, weight: "bold")[`document()` can silently repoint cSEM's own help links]\
  #text(size: 10pt, style: "italic")[
    Why `[fit()]` became `generics::fit()` #sym.dash.em and what to do about it
  ]
]

= Summary

Running `devtools::document()` rewrote twelve `.Rd` files that had nothing to do
with the change being documented. Every rewrite was of the same shape: a link to
one of cSEM's *own* functions was repointed at a different package.

```
-\code{\link[=fit]{fit()}}
+\code{\link[generics:fit]{generics::fit()}}

-\code{\link[=summarize]{summarize()}}
+\code{\link[dplyr:summarize]{dplyr::summarize()}}
```

This is not a roxygen bug and not a source-code bug. It is a *lookup* that
depends on the state of your R library rather than on the contents of your
package, and it fails open --- rewriting links rather than erroring. The
corruption is easy to miss because it appears in files you did not edit, and it
is only visible in `git diff`.

Two things follow. The exposure is *small and known*: exactly two topics in the
whole package can be affected, at 16 call sites. And the fix is *local*: qualify
those 16 links and the fragile lookup is never reached.

= The mechanism

`roxygen2:::resolve_link_package()` decides what a markdown link `[foo()]`
becomes:

```r
me <- me %||% roxy_meta_get("current_package")
if (has_topic(topic, me)) return(NA_character_)      # (1) local link
pkgs <- local_pkg_deps(pkgdir)                       # (2) otherwise, scan deps
pkg_has_topic <- pkgs[map_lgl(pkgs, has_topic, topic = topic)]
...
} else if (length(pkg_has_topic) == 1) {
  if (pkg_has_topic %in% base) return(NA_character_) # (3) base packages exempt
  else return(pkg_has_topic)                         # (4) qualify to that package
}
...
warn_roxy_tag(... "Could not resolve link to topic ...")   # (5) nothing owns it
```

Step (1) is the one that fails. `has_topic()` is a `help()` lookup wrapped in a
`tryCatch` that swallows *any* error into `FALSE`:

```r
has_topic <- function(topic, package) {
  tryCatch({
    out <- exec("help", topic, package, .env = global_env())
    inherits(out, "dev_topic") ||
      (inherits(out, "help_files_with_topic") && length(out) == 1)
  }, error = function(c) FALSE)
}
```

Note `.env = global_env()`: `help` is resolved from the *global* environment, so
which lookup runs depends on the session, not on your sources.

- With `pkgload`'s shims attached (the usual `load_all()` case) it searches the
  *source* `man/` directory and returns a `dev_topic`.
- Without them it searches the *installed* library and must return exactly one
  match.

If neither succeeds --- the package is not installed in the active library, the
shims are absent, the topic is ambiguous, or `help()` errors for any other
reason --- `has_topic()` returns `FALSE`, and control falls to step (2), which is
where the damage happens.

*`local_pkg_deps()` searches `Suggests` too:*

```r
deps <- deps[deps$type %in% c("Depends", "Imports", "Suggests"), ]
```

This is the detail that makes the diagnosis certain. `dplyr` is only a *Suggests*
dependency of cSEM. Nothing but the fallback branch could have produced
`dplyr::summarize()`, so the observed corruption proves step (1) returned `FALSE`
for `summarize` during that run.

= Reproduction

The trigger is a failing topic lookup *for the package being documented* while
the installed dependencies still resolve normally. Shadowing `help()` for that
one package reproduces it deterministically, without having to recreate a broken
library:

```r
deps <- roxygen2:::local_pkg_deps(".")
base <- roxygen2:::base_packages()

resolve <- function(t) {
  if (isTRUE(roxygen2:::has_topic(t, "cSEM")))
    return(paste0("\\link[=", t, "]{", t, "()}   <- local"))
  owners <- setdiff(deps[vapply(deps, function(p)
    isTRUE(roxygen2:::has_topic(t, p)), logical(1))], base)
  if (length(owners) == 1)
    paste0("\\link[", owners, ":", t, "]{", owners, "::", t, "()}   <- REPOINTED")
  else if (length(owners) == 0) "unresolved (warning)"
  else paste0("ambiguous: ", paste(owners, collapse = ", "))
}

topics <- c("fit", "summarize", "csem", "assess")
for (t in topics) cat(sprintf("  %-10s %s\n", t, resolve(t)))

## Fail the lookup for cSEM only -- what an uninstalled dev package, or a
## shimmed help() that cannot see man/, does to has_topic().
assign("help", function(topic, package, ...) {
  if (identical(package, "cSEM")) stop("topic lookup unavailable")
  do.call(utils::help, list(topic, package))
}, envir = globalenv())

for (t in topics) cat(sprintf("  %-10s %s\n", t, resolve(t)))
```

Measured output:

```
healthy session:
  fit        \link[=fit]{fit()}   <- local
  summarize  \link[=summarize]{summarize()}   <- local
  csem       \link[=csem]{csem()}   <- local
  assess     \link[=assess]{assess()}   <- local

lookup fails for cSEM only:
  fit        \link[generics:fit]{generics::fit()}   <- REPOINTED
  summarize  \link[dplyr:summarize]{dplyr::summarize()}   <- REPOINTED
  csem       unresolved (warning)
  assess     unresolved (warning)
```

The first two lines of the second block are byte-identical to what appeared in
`git diff`, which is what identifies this as the mechanism rather than a
plausible story about it.

*The natural triggers*, none of which involve editing anything:

- Documenting a package that is *not installed* in the active library --- a fresh
  clone, a CI job, a new machine, or a library switched by `renv`.
- A session where `pkgload`'s `help` shim is not attached, combined with an
  installed copy that lacks the topic (an older install, or one built before the
  function existed).
- Any `help()` error at all: the `tryCatch` does not distinguish "not found" from
  "something went wrong".

*What was not established.* The condition that produced the original corruption
in this repository is not reproducible here: `has_topic("fit", "cSEM")` and
`has_topic("summarize", "cSEM")` both return `TRUE` in this session under both
lookup modes, and re-running `document()` now touches only the files it should.
The code path is identified with certainty; the environmental trigger on that
particular run is not.

= Exposure: exactly two topics, sixteen sites

Scanning every `[foo()]` markdown link in `R/*.R` and asking which of them a
non-base dependency also owns:

```
distinct bare [foo()] targets: 68     total occurrences: 296

AT RISK (topic also owned by a non-base dependency):
  fit          x4    -> generics
  summarize    x12   -> dplyr
```

The other 66 targets are safe, for one of two reasons. Most are owned by no
dependency at all, so step (4) is unreachable. The rest --- `predict` (stats),
`plot` (graphics) --- are owned only by *base* packages, and step (3) exempts
those explicitly, returning a local link.

#table(
  columns: (auto, 1fr, auto),
  inset: 6pt,
  align: (left, left, right),
  table.header([*Topic*], [*File*], [*Sites*]),
  [`fit`],       [`R/helper_doTrees.R`],       [1],
  [`fit`],       [`R/zz_arguments.R`],         [3],
  [`summarize`], [`R/00_csem.R`],              [2],
  [`summarize`], [`R/csem_resample.R`],        [1],
  [`summarize`], [`R/exportToExcel.R`],        [2],
  [`summarize`], [`R/helper_summarize.R`],     [2],
  [`summarize`], [`R/postestimate_infer.R`],   [2],
  [`summarize`], [`R/postestimate_verify.R`],  [1],
  [`summarize`], [`R/print.cSEMSummarize.R`],  [1],
  [`summarize`], [`R/zz_arguments.R`],         [1],
)

Both are cSEM's own exported functions with their own `.Rd` files
(`R/csem_fit.R:48` and `R/postestimate_summarize.R:41`), which is precisely why
the repointing is wrong rather than merely ugly: the reader is sent to
`generics::fit`, a stub generic that documents none of cSEM's arguments.

= Fixes

== Recommended: qualify the sixteen links

A qualified link is handled by `resolve_qualified_link()` and never consults
`has_topic()`, so the fragile branch becomes unreachable. Use the two-part form,
which pins the target while keeping the displayed text:

```r
#' @seealso [fit()][cSEM::fit], [summarize()][cSEM::summarize]
```

renders as `fit()` and `summarize()` but generates `\link[cSEM:fit]{fit()}`. The
one-part form `[cSEM::fit()]` also works, but renders the qualification visibly
as `cSEM::fit()`, which is noisier in prose.

This is roxygen's own advice --- its unresolved-link warning ends with
"Alternatively, you can fully qualify the link with a package name" --- and it is
a 16-line change confined to `@seealso` and `@details` text.

== Mitigation: install before documenting

`R CMD INSTALL .` (or `devtools::install()`) before `document()` makes the
installed-library lookup succeed, which keeps step (1) on the local branch. This
works, but it is a property of the machine rather than of the repository: it has
to be remembered by every contributor, on every machine, forever. Prefer it as a
habit, not as the fix.

== Detection: read `git status` after documenting

Whatever else is done, `document()` should never modify a `.Rd` file whose source
you did not touch. Any such file in `git status` is a signal to inspect rather
than stage:

```
git checkout -- $(git diff --name-only man/ | grep -v <files you edited>)
```

That is how the twelve corrupted files were caught here.

== Not a fix: adding a wordlist entry or a NAMESPACE import

Both were considered and neither addresses this. The links are wrong regardless
of spelling, and `importFrom(generics, fit)` would make the situation *worse* by
adding a genuine competing binding for the topic name.

= Related, but different

`devtools::check_man()` also reports unresolved links from `R/00_csem.R` to
`parseModel`, `cSEMModel`, `calculateWeightsPLS` and `calculateWeightsGSCA`.
These emit the same warning from step (5), but the cause is not the same: those
topics are genuinely undocumented (internal functions with no `.Rd`), so no
lookup could ever succeed. The fix there is either to document them or to stop
linking to them --- not qualification.

The distinction matters when reading a `document()` run: a link to a *documented,
exported* function that resolves elsewhere is the failure described here; a link
to an internal function that resolves nowhere is an ordinary missing topic.

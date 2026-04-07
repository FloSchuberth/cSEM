# Plan: Reducing cSEM Package Imports

**Date:** 2026-04-07
**Branch:** igsca-finetuning
**Current import count:** 28 packages in `Imports:`

## Current Imports Audit

| Package | License | Functions Used | Usage Scope | Verdict |
|---------|---------|---------------|-------------|---------|
| `alabama` | GPL (>= 2) | `auglag()` | 1 call in `estimators_weights.R` | **Keep** — constrained optimization, hard to replace |
| `cli` | MIT | `boxx()`, `rule()` | ~15 calls across print files | **Keep** — actively used for formatted output |
| `crayon` | MIT | `col_align()`, `col_nchar()`, `bold()`, `%+%` | Heavily used via `import(crayon)` in all print files | See note below |
| `expm` | GPL (>= 2) | `sqrtm()` | 3 calls in 2 files | **Keep** — matrix square root is non-trivial |
| `future` | LGPL (>= 2.1) | `plan()` | Core parallelization infra | **Keep** — essential for parallel resampling |
| `future.apply` | GPL (>= 2) | `future_lapply()` | Core parallelization infra | **Keep** — essential for parallel resampling |
| `generics` | MIT | `tidy`, `glance` | Re-exported S3 generics | **Keep** — needed for broom-style tidiers |
| `lavaan` | GPL (>= 2) | `cfa()`, `lavPredict()`, model syntax | Core dependency for FSR estimation | **Keep** — deeply integrated |
| `lifecycle` | MIT | `deprecate_soft` (imported but never called) | Only `\lifecycle{}` Rd badges | **REMOVE candidate** |
| `magrittr` | MIT | `%>%` | ~152 uses across 11 files | See note below |
| `MASS` | GPL-2 \| GPL-3 | `ginv()`, `cov.rob()` | ~15 calls across 3 files | **Keep** — generalized inverse is critical |
| `Matrix` | GPL (>= 2) | `bdiag()` | 3 calls in 2 files | **Keep** — block diagonal matrices |
| `matrixcalc` | GPL (>= 2) | `matrix.trace()`, `is.positive.semi.definite()`, `is.symmetric.matrix()`, `is.square.matrix()` | ~11 calls in 5 files | **REPLACE candidate** |
| `matrixStats` | Artistic-2.0 | `rowProds()`, `colSds()`, `colQuantiles()` | ~25 calls in 4 files | See note below |
| `polycor` | GPL (>= 2) | `polychor()`, `polyserial()` | 3 calls in `helper_foreman.R` | **Keep** — specialized correlation types |
| `progressr` | GPL (>= 3) | `with_progress()`, `progressor()` | 4 files for progress bars | **Keep** — user-facing progress reporting |
| `psych` | GPL (>= 2) | `principal()` | 2 calls in 2 files | **REPLACE candidate** |
| `purrr` | MIT | `map()`, `transpose()`, `modify_depth()` | ~25 calls, mostly in `helper_test_MGD.R` | See note below |
| `Rdpack` | GPL (>= 2) | `\insertRef{}`, `\insertCite{}` | Hundreds of Rd references | **Keep** — essential for documentation |
| `rlang` | MIT | `.data` pronoun | 1 import + conditional require | **Keep** — `.data` pronoun used in tidy eval |
| `symmoments` | GPL | `callmultmoments()` | 2 calls in `helper_estimators_paths.R` | **Keep** — multivariate moment computation |
| `TruncatedNormal` | GPL-3 | `mvrandn()` | 1 call in `postestimate_predict.R` | **Keep** — specialized distribution sampling |
| `partykit` | GPL-2 \| GPL-3 | `mob()`, `mob_control()` | `postestimate_doTrees.R` | **Keep** — core to SEM trees |
| `Formula` | GPL-2 \| GPL-3 | `as.Formula()` (imported but apparently unused in code) | Only `@importFrom` tag exists | **REMOVE candidate** |
| `boot` | Unlimited | `boot()` | 1 call in `postestimate_doTrees.R` | See note below |
| `tibble` | MIT | `as_tibble()`, `tibble()` | Tidier output + MGD helpers | **Keep** — used in exported S3 methods |
| `stats` | (base R) | — | — | N/A (base) |
| `utils` | (base R) | — | — | N/A (base) |

## Overlapping Function Analysis

When multiple imported packages provide similar functionality, prefer consolidating onto the more stable/major upstream package rather than importing from several.

### `crayon` vs `cli` (Terminal Formatting)
- **Both imported.** `cli` is the modern successor to `crayon` (same maintainer, Gábor Csárdi).
- `crayon` provides: `col_align()`, `col_nchar()`, `bold()`, `%+%` — used heavily in print methods.
- `cli` provides: `boxx()`, `rule()` — used in ~15 places for hypothesis box formatting.
- `cli` also has `style_bold()`, `ansi_nchar()`, `ansi_align()` which can replace crayon equivalents.
- **Recommendation:** Consolidate onto `cli` (the actively maintained package). `crayon` is in maintenance mode.

### `matrixcalc` vs base R (Matrix Properties)
- `matrixcalc` provides: `matrix.trace()`, `is.symmetric.matrix()`, `is.square.matrix()`, `is.positive.semi.definite()`.
- Base R already has: `isSymmetric()`, and the rest are trivial one-liners (`sum(diag(A))`, `nrow(A)==ncol(A)`, eigenvalue check).
- **No overlap with `Matrix` or `MASS`** — those provide different functionality (`bdiag()` and `ginv()`/`cov.rob()` respectively).
- **Recommendation:** Replace with base R helpers. No need for a separate package for these.

### `purrr` vs base R (List Manipulation)
- `purrr::map()` → direct equivalent: `lapply()`
- `purrr::transpose()` → no base equivalent, but a ~5 line internal helper suffices
- `purrr::modify_depth()` → no base equivalent, ~10 line recursive helper
- **Recommendation:** Replace with base R. `purrr` is a large tidyverse dependency for just 3 functions.

### `magrittr` vs base R pipe (Piping)
- `magrittr` provides `%>%`. Base R `|>` is available since R 4.1.0 (cSEM's minimum).
- **Caveat:** Some uses may rely on magrittr's `.` placeholder. The base pipe `|>` uses `_` as placeholder (R 4.2.0+) or anonymous functions.
- **Recommendation:** Can be replaced, but scope is large (~152 uses). Medium priority.

### `psych` vs base R (PCA)
- `psych::principal()` with `nfactors=1` is equivalent to extracting the first eigenvector from `eigen()`.
- **Recommendation:** Replace with base R eigen decomposition. `psych` is a very large package for one function call.

### No Overlap Found
- `MASS` (`ginv`, `cov.rob`) — unique functionality, no overlap with other imports
- `Matrix` (`bdiag`) — unique block diagonal construction
- `expm` (`sqrtm`) — unique matrix square root
- `lavaan` — unique SEM functionality
- `polycor` — unique polychoric/polyserial correlations
- `symmoments` — unique multivariate moments
- `TruncatedNormal` — unique truncated MVN sampling
- `partykit` / `Formula` / `boot` — unique tree/model infrastructure

## Recommended Changes (Ranked by Simplicity)

### Tier 1: Easy Removals (no code changes or trivial changes)

#### 1. Remove `lifecycle`
- **Effort:** Trivial
- **Rationale:** `deprecate_soft` is imported but **never actually called** anywhere in R code. The `\lifecycle{}` Rd badges are just documentation macros that work through the `RdMacros` field — lifecycle just needs to be in `RdMacros:` (and possibly `Suggests:`), not `Imports:`.
- **Action:** Move `lifecycle` from `Imports:` to `Suggests:`. Keep in `RdMacros:`. Remove `@importFrom lifecycle deprecate_soft` from `zz_package.R`.
- **Risk:** Low. The Rd badges will still render. If `deprecate_soft` is ever needed in future, it can be called with `lifecycle::deprecate_soft()`.
- **License:** MIT — not a factor.

#### 2. Remove `Formula`
- **Effort:** Trivial
- **Rationale:** `as.Formula` is declared via `@importFrom` in `postestimate_doTrees.R` but **never actually called** in any R code. It may have been needed by `partykit::mob()` at one point but partykit handles this internally.
- **Action:** Remove `Formula` from `Imports:` in DESCRIPTION. Remove `@importFrom Formula as.Formula` from `postestimate_doTrees.R`. Run `devtools::document()`.
- **Risk:** Low. Test `doTrees()` to confirm partykit doesn't rely on Formula being loaded by cSEM.
- **License:** GPL-2 | GPL-3 — removal reduces GPL surface.

### Tier 2: Replace with Base R (small code changes)

#### 3. Replace `matrixcalc` with internal helpers
- **Effort:** Small (write 4 small functions)
- **Rationale:** Only 4 functions used, all trivially implementable in base R:
  - `matrix.trace(A)` → `sum(diag(A))`
  - `is.symmetric.matrix(A)` → `isSymmetric(A)` (base R)
  - `is.square.matrix(A)` → `nrow(A) == ncol(A)`
  - `is.positive.semi.definite(A)` → `all(eigen(A, symmetric = TRUE, only.values = TRUE)$values >= -sqrt(.Machine$double.eps))`
- **Action:** Create internal helper functions (or inline). Replace ~11 call sites across 5 files.
- **License:** GPL (>= 2) — compatible, but one less dependency is cleaner.
- **Risk:** Low. The base R replacements are mathematically equivalent. Need careful tolerance handling for PSD check.

#### 4. Replace `psych::principal()` with base R eigen decomposition
- **Effort:** Small (write 1 helper)
- **Rationale:** `psych::principal()` is called with `nfactors = 1`, which is just extracting the first principal component. This is `eigen(S)$vectors[,1]` scaled appropriately.
- **Action:** Write an internal `first_principal_component()` function. Replace 2 call sites.
- **License:** GPL (>= 2) — compatible, but psych is a large package to depend on for one function.
- **Risk:** Medium. Need to verify the exact scaling/sign convention matches `psych::principal()`.

### Tier 3: Moderate Effort Consolidations

#### 5. Replace `magrittr` pipe with base R pipe `|>`
- **Effort:** Medium (~152 replacements across 11 files)
- **Rationale:** R >= 4.1.0 (already the minimum for cSEM) includes the native pipe `|>`. The magrittr `%>%` could be replaced with `|>` everywhere.
- **Caveats:**
  - `|>` doesn't support the `.` placeholder the same way. Uses like `x %>% foo(., bar)` need rewriting as `x |> (\(z) foo(z, bar))()` or restructured.
  - The `%+%` operator from crayon (string concatenation) is NOT from magrittr, so that stays.
  - `helper_test_MGD.R` has 104 pipe uses — largest single file.
- **Action:** Systematic replacement. Can be done file-by-file.
- **Risk:** Medium. Some pipe chains may use `.` placeholder in non-trivial ways (especially in `helper_test_MGD.R` which also uses `dplyr`/`tidyr`).
- **License:** MIT — not a licensing concern.

#### 6. Reduce `purrr` usage
- **Effort:** Medium
- **Rationale:** Only 3 functions used: `map()`, `transpose()`, `modify_depth()`.
  - `purrr::map()` → `lapply()`
  - `purrr::transpose()` → write a small internal `list_transpose()` helper
  - `purrr::modify_depth()` → write a recursive internal helper or restructure
- **Caveats:** `modify_depth()` is used in complex nested list operations in `helper_test_MGD.R`. Replacing these requires careful testing.
- **Action:** Replace `map()` calls with `lapply()`. Write internal `list_transpose()`. Carefully rewrite `modify_depth()` calls.
- **License:** MIT — copyable with attribution if needed.

#### 7. Consider `crayon` + `cli` consolidation
- **Effort:** Large
- **Rationale:** `crayon` is used via full `import(crayon)` for `col_align()`, `col_nchar()`, `bold()`, and `%+%`. `cli` is used for `boxx()` and `rule()`. These are both terminal formatting packages. `cli` actually supersedes `crayon` — cli can do everything crayon does.
- **Action:** Replace crayon functions with cli equivalents. `bold()` → `cli::style_bold()`, `%+%` → `paste0()`, etc. Then remove crayon.
- **Risk:** Medium-high. Crayon is deeply used throughout all print methods. The `col_align()` and `col_nchar()` functions handle ANSI-aware string width which would need equivalent cli calls.
- **License:** Both MIT.

### Tier 4: Not Recommended to Remove

These packages are either deeply integrated, provide critical functionality, or would require significant rewrites:

- **`lavaan`** — Core dependency for CFA/FSR estimation and model syntax
- **`MASS`** — `ginv()` (generalized inverse) used ~15 times; `cov.rob()` for robust correlation
- **`Matrix`** — `bdiag()` for block diagonal construction
- **`expm`** — `sqrtm()` for matrix square root (no simple base R equivalent)
- **`alabama`** — `auglag()` for augmented Lagrangian optimization
- **`future` / `future.apply`** — Entire parallelization infrastructure
- **`progressr`** — Progress bar system integrated with future
- **`polycor`** — Polychoric/polyserial correlations (specialized)
- **`symmoments`** — Multivariate moment calculations (specialized)
- **`TruncatedNormal`** — Truncated multivariate normal sampling (specialized)
- **`Rdpack`** — Bibliography management in Rd files (hundreds of references)
- **`rlang`** — `.data` pronoun for tidy evaluation
- **`generics`** — Provides `tidy()` and `glance()` generics
- **`tibble`** — Used in exported S3 methods (`tidy.cSEMResults`, `glance.cSEMResults`)
- **`partykit`** — Core to `doTrees()` functionality
- **`boot`** — Used in `doTrees()` for parametric bootstrap of fit indices. License is "Unlimited" so no concern, and it ships with R. Low priority for removal.

### `matrixStats` Note
- **License:** Artistic-2.0 (permissive, compatible with GPL-3)
- `rowProds()` is used ~20 times and is performance-critical (vectorized C code). A base R replacement (`apply(x, 1, prod)`) would be slower.
- `colSds()` and `colQuantiles()` are used a few times and could be replaced with `apply()` variants.
- **Verdict:** Keep for now. The performance benefit of `rowProds()` matters for resampling loops.

## Summary: Achievable Import Reduction

| Change | Packages Removed | Effort | Risk |
|--------|-----------------|--------|------|
| Remove `lifecycle` from Imports | 1 | Trivial | Low |
| Remove `Formula` from Imports | 1 | Trivial | Low |
| Replace `matrixcalc` | 1 | Small | Low |
| Replace `psych` | 1 | Small | Medium |
| Replace `magrittr` with `|>` | 1 | Medium | Medium |
| Replace `purrr` | 1 | Medium | Medium |
| Consolidate `crayon` into `cli` | 1 | Large | Medium-high |

**Quick wins (Tier 1):** Remove 2 imports with essentially no code changes.
**Small effort (Tier 1+2):** Remove 4 imports with minor code changes.
**Full effort (Tier 1-3):** Remove up to 7 imports (28 → 21).

## Recommended Execution Order

1. `lifecycle` — trivial, do first
2. `Formula` — trivial, do second
3. `matrixcalc` — small, well-contained replacements
4. `psych` — small but needs careful validation
5. `magrittr` → `|>` — mechanical but large scope
6. `purrr` — needs careful testing of nested list operations
7. `crayon` → `cli` — largest scope, do last

Each step should be a separate commit with tests run after.

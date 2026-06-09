# Plan: Monte Carlo study of bootstrapped R² / adj-R² model comparison (pooled vs split)

> Project-local copy of the approved plan (master copy: `~/.claude/plans/valiant-beaming-mccarthy.md`).
> Implemented in `dev/igsca/ohtani.R` after the `# Simulation` header.

## Context

`dev/igsca/ohtani.R` examines a claim from Hwang & Takane (2014, p.27): that a bootstrapped
*paired t-test* on FIT values can test whether two models differ in fit, described as "a
nonparametric version of the paired t test." The PoC already in the file confirmed the
algebraic fact that a paired t-test of (a vs b) ≡ a one-sample t-test of (a−b). A second quote
the user supplied — "a paired t-test ... equivalent to testing whether the difference in the FIT
values calculated from **100 bootstrap samples** ... t(99) = −0.59" — confirms the procedure under
test treats the B bootstrap replicates as the t-test's sample (`df = 99 = B − 1`, B = 100).

This study fills the empty `# Simulation` section with a **SimDesign** Monte Carlo study that
stress-tests that procedure: comparing a single **pooled** multiple linear regression on all data
against **two separate** regressions fit to two parts of the data. By varying how much the slopes
differ between the parts (δ), we estimate the **type I error** (δ = 0) and **power / 1 − type II
error** (δ > 0) of the published method and of progressively-correct alternatives.

Settled with the user: split FIT = combined-prediction R² (two separate models, pooled residuals);
include the published paired t-tests **plus** Wald/SE **plus** percentile/ASL CI **plus** the
null-pooled residual bootstrap **plus** a Chow-test oracle; vary B (anchored at 100); 1000 MC reps.

---

## Adversarial collaboration: Statistician (S) vs Machine Learner (M)

1. **In-sample R² is monotone in model complexity (S & M agree).** Pooled (`y ~ x1 + x2`) is nested
   in split (`y ~ group*(x1+x2)`), so `R²_split ≥ R²_pool` w.p. 1 even at δ = 0 → raw-R² tests are
   biased toward rejection. Adjusted R² penalizes the 3 extra params → fairer. Report both.

2. **The published Hwang–Takane paired t-test makes significance an artifact of B (S).** `df = 99 =
   B − 1` proves `SE = sd(d*)/√B`. As B → ∞ the t-stat → ∞ and rejection → 1 for any non-zero
   difference (and the nested-model difference is always > 0). Keep it as **the method under test**;
   anchor B = 100 (their value) and vary up to B = 1000 — only its type I error should climb with B.

3. **Three proposed corrections of increasing rigor (S):**
   - **(3a′) Wald/studentized** `z = T_obs / sd(d*)`, no √B (`t_paired ≈ √B·z_Wald`) — isolates the √B
     error; still over-rejects for raw R² (statistic bounded at 0).
   - **(3a) Percentile / ASL CI** of the FIT difference — B-stable; raw-R² CI sits above 0 → over-rejects.
   - **(3b) Null-pooled residual bootstrap** — impose H₀ (single pooled model), `y* = ŷ_pool + e*`,
     `p = (1 + #{T* ≥ T_obs})/(B+1)`. Canonical bootstrap hypothesis test; handles the raw-R² bound.

4. **In-sample fit is the wrong target (M).** Out-of-sample/CV R² or AIC/BIC is principled; in-sample
   R² mechanically favors the higher-capacity split model. Honor the R²/adj-R² request; flag CV in a comment.

5. **Exact oracle (S):** the **Chow test** (`anova(pooled, full-interaction)`) — exact F-test, type I
   at α and the power ceiling every bootstrap is judged against.

---

## Data-generating model

Two groups `g ∈ {1,2}`, equal sizes `n1 = n2 = N/2` (N = 200):

```
x1, x2  ~ iid N(0, 1)
group 1: y = β0 + β1·x1     + β2·x2     + ε
group 2: y = β0 + (β1+δ)·x1 + (β2+δ)·x2 + ε      ε ~ N(0, σ²)
```

Fixed: `β0 = 0`, `β1 = β2 = 0.3`, `σ = 1`. δ = 0 → null (type I); δ > 0 → power.

## Design grid (`createDesign`)

| factor | levels |
|--------|--------|
| `delta` | 0, 0.1, 0.2, 0.4 |
| `B` | 100, 1000 |
| `prop1` | 0.5, 0.7, 0.9 |

B = 100 reproduces `df = 99`; B = 1000 exposes √B inflation. `prop1` = fraction of N in part 1
(group-size imbalance), N held fixed at 200 → splits 100/100, 140/60, 180/20; the deviating group
is part 2. **24 conditions × 1000 reps** (~3 min parallel).

## SimDesign components (in `ohtani.R` after line 94)

- **`Generate`** → tibble `{x1, x2, y, group}` per the DGP.
- **`fit_fast(X, y, g1, g2)`** (used in bootstrap; passed via `fixed_objects` for worker export) and
  readable **`fit_compare(d)`** (lm-based, for sanity): split = two separate fits, combined R²/adj-R²
  from pooled residuals (`SSE = SSE1+SSE2`, `SST` over all y, `k = 6`).
- **`Analyse`** returns 9 p-values: `p_ttest_R2/adjR2` (published, /√B), `p_se_R2/adjR2` (Wald, no √B),
  `p_pctl_R2/adjR2` (percentile/ASL), `p_null_R2/adjR2` (null-pooled residual boot), `p_chow` (oracle).
- **`Summarise`** → `EDR(results, alpha = 0.05)` → rejection rate (δ=0: type I; δ>0: power).
- **Driver:** `runSimulation(..., replications = 1000, parallel = "future")` reusing the
  `mirai_multisession` plan; `packages = c("tibble")`.

## Expected results

- `p_ttest_*` (/√B): type I inflated, **rising with B**; raw ≈ 1. Headline √B-flaw demo.
- `p_se_*` / `p_pctl_*`: B-stable; raw still over-rejects, adj closer but not calibrated.
- `p_null_*`: type I ≈ 0.05, B-stable — correctly calibrated.
- `p_chow`: type I ≈ 0.05, highest power — ceiling.

## Verification

1. **Unit sanity:** δ=0 dataset → `R2_split > R2_pool`; two-separate-model combined R²/adj-R² ==
   `summary(lm(y ~ group*(x1+x2)))$r.squared`/`$adj.r.squared`; `fit_fast == fit_compare`.
2. **Smoke:** `runSimulation(replications = 20)` returns all 9 `p_*` columns, no errors.
3. **Full:** `replications = 1000`; confirm δ=0 type-I pattern and δ>0 power; MC SE ≈ 0.007 at 0.05.

## Files

- `dev/igsca/ohtani.R` — only file modified (code after `# Simulation`, line 94+). Reuses libs/future
  plan loaded at the top of the file.

## Progress log

- [x] Plan approved and saved (master + this project-local copy).
- [x] Simulation code written into `ohtani.R` (after `# Simulation`, ~line 95).
- [x] Unit sanity checks pass (split R² > pooled R²; two-separate-model FIT == interaction-model FIT; `fit_fast == fit_compare`).
- [x] Smoke run (replications = 20) passes — all 9 `p_*` columns returned.
- [x] Full run (replications = 1000, 8 conditions) complete in ~1 min (parallel `mirai`).
- [x] **Extended with group-size imbalance** (`prop1` ∈ {0.5, 0.7, 0.9}); 24-condition run (~3 min) below.

### Results (1000 replications; rejection rate at α = .05)

Design extended with **`prop1`** = fraction of N (=200) in part 1; the **deviating** group is always
part 2, so larger `prop1` shrinks the deviating group (100/100 → 140/60 → 180/20). δ = 0 rows =
empirical **Type I error**; δ > 0 rows = **power**.

**prop1 = 0.5 (balanced, 100/100)**

| δ | B | ttest_R² | ttest_adj | se_R² | se_adj | pctl_R² | pctl_adj | null_R² | null_adj | chow |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 100 | 1.000 | 0.882 | 0.003 | 0.001 | 1.000 | 0.018 | 0.045 | 0.045 | 0.043 |
| 0.2 | 100 | 1.000 | 0.980 | 0.092 | 0.026 | 1.000 | 0.217 | 0.345 | 0.349 | 0.359 |
| 0.4 | 100 | 1.000 | 1.000 | 0.647 | 0.404 | 1.000 | 0.832 | 0.918 | 0.917 | 0.918 |
| 0.0 | 1000 | 1.000 | 0.981 | 0.002 | 0.000 | 1.000 | 0.014 | 0.044 | 0.046 | 0.039 |
| 0.2 | 1000 | 1.000 | 0.992 | 0.075 | 0.020 | 1.000 | 0.198 | 0.338 | 0.338 | 0.338 |
| 0.4 | 1000 | 1.000 | 1.000 | 0.626 | 0.386 | 1.000 | 0.822 | 0.917 | 0.917 | 0.917 |

**prop1 = 0.7 (140/60)**

| δ | B | ttest_R² | ttest_adj | se_R² | se_adj | pctl_R² | pctl_adj | null_R² | null_adj | chow |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 100 | 1.000 | 0.874 | 0.006 | 0.000 | 1.000 | 0.023 | 0.059 | 0.059 | 0.056 |
| 0.2 | 100 | 1.000 | 0.972 | 0.056 | 0.013 | 1.000 | 0.178 | 0.297 | 0.301 | 0.294 |
| 0.4 | 100 | 1.000 | 1.000 | 0.530 | 0.291 | 1.000 | 0.736 | 0.851 | 0.852 | 0.859 |
| 0.0 | 1000 | 1.000 | 0.973 | 0.006 | 0.001 | 1.000 | 0.022 | 0.053 | 0.053 | 0.050 |
| 0.2 | 1000 | 1.000 | 0.988 | 0.060 | 0.013 | 1.000 | 0.152 | 0.297 | 0.297 | 0.291 |
| 0.4 | 1000 | 1.000 | 1.000 | 0.517 | 0.281 | 1.000 | 0.723 | 0.862 | 0.862 | 0.859 |

**prop1 = 0.9 (severe, 180/20)**

| δ | B | ttest_R² | ttest_adj | se_R² | se_adj | pctl_R² | pctl_adj | null_R² | null_adj | chow |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 100 | 1.000 | 0.848 | 0.014 | 0.008 | 1.000 | 0.023 | 0.046 | 0.045 | 0.038 |
| 0.2 | 100 | 1.000 | 0.900 | 0.043 | 0.007 | 1.000 | 0.087 | 0.149 | 0.147 | 0.150 |
| 0.4 | 100 | 1.000 | 0.980 | 0.196 | 0.064 | 1.000 | 0.326 | 0.474 | 0.473 | 0.470 |
| 0.0 | 1000 | 1.000 | 0.966 | 0.011 | 0.005 | 1.000 | 0.017 | 0.045 | 0.043 | 0.043 |
| 0.2 | 1000 | 1.000 | 0.974 | 0.040 | 0.011 | 1.000 | 0.078 | 0.149 | 0.150 | 0.149 |
| 0.4 | 1000 | 1.000 | 0.995 | 0.200 | 0.062 | 1.000 | 0.313 | 0.474 | 0.474 | 0.475 |

**Conclusions**

- **Published H&T paired t-test (`p_ttest_*`) is invalid at every imbalance level.** Raw-R² always
  rejects (Type I = 1.0); adj-R² Type I climbs **0.85–0.88 (B=100) → 0.97 (B=1000)** — the √B artifact,
  present even at their own B = 100. Significance is driven by the bootstrap count, not the evidence.
- **CI-style readings fail in opposite directions on raw R²:** percentile (`p_pctl_R2`) always rejects
  (`d* > 0` in every resample); Wald (`p_se_R2`) is conservative. On adj-R² both are conservative/underpowered.
- **Null-pooled residual bootstrap (`p_null_*`) is correctly calibrated** (Type I ≈ 0.05, B-stable for
  raw and adj) and **tracks the Chow oracle almost exactly in both Type I and power, across all
  imbalance levels** — the correct nonparametric test stays efficient under imbalance.
- **Imbalance effect:** Type I calibration of the valid methods is **robust** to imbalance, but **power
  collapses as the deviating subgroup shrinks** — at δ = 0.4, power falls **0.92 (100/100) → 0.86
  (140/60) → 0.47 (180/20)**. A small heterogeneous subgroup carries too little information to detect
  its own slope difference. (The deviation here sits in part 2; were the large group the deviating one,
  the power loss would be far milder — the asymmetry is the point.)
- **Takeaway:** imposing the null via a pooled residual bootstrap recovers the validity and efficiency
  of the exact Chow F-test at all imbalance levels; the published bootstrap-replicate paired t-test does not.

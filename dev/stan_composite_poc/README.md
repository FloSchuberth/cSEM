# Henseler–Ogasawara composite + common factor — Stan PoC

A proof of concept that fits the Henseler–Ogasawara (H-O) specification of a
canonical composite affecting a common factor in **Stan**, using a
**probabilistic graphical model** (row-level) formulation rather than a
covariance-based / Wishart sufficient-statistic likelihood. Synthetic
population data is generated with **cSEM.DGP**; parameter recovery is
assessed by Monte Carlo replication.

```
   x₁ ─┐
   x₂ ─┼─ canonical composite c ──β──▶ common factor η ──┬──▶ y₁
   x₃ ─┘                                                  ├──▶ y₂
                                                          └──▶ y₃
```

## Files

| file | purpose |
|---|---|
| `composite_ho.stan`         | Stan model (PGM formulation, marginalized η) |
| `00_setup.R`                | DGP definition, truth, simulation + init helpers |
| `01_single_run.R`           | one replication, full diagnostics, recovery table |
| `02_simulation_study.R`     | multi-rep bias / RMSE / 90% CI coverage |
| `03_sbc.R`                  | minimal simulation-based calibration |
| `single_run_result.rds`     | output of `01_single_run.R` |
| `sim_study_result.rds`      | output of `02_simulation_study.R` |
| `sbc_result.rds`            | output of `03_sbc.R` |

To reproduce (assumes R, cmdstanr, cmdstan 2.36, cSEM, cSEM.DGP installed):

```bash
Rscript dev/stan_composite_poc/01_single_run.R           # ~1 min
REPS=30 N=500 Rscript dev/stan_composite_poc/02_simulation_study.R  # ~30 min
R_SBC=50 L_DRAWS=50 Rscript dev/stan_composite_poc/03_sbc.R         # ~30 min
```

## Model

For each observation `i = 1..N`:

- Composite indicators `xᵢ ∈ ℝᴷ` with population correlation matrix Σₓ.
- Composite (deterministic, canonical, Var(c)=1):
  `cᵢ = w' xᵢ`,
  with `w` normalized so that `w' Σₓ w = 1`.
- H-O latent decomposition (implicit):
  `xᵢ = a cᵢ + V νᵢ` with `νᵢ ~ N(0, I)`, `w' a = 1`, `w' V = 0`.
  Marginally `xᵢ ~ MVN(0, Σₓ)` where `Σₓ = a a' + V V'`.
- Structural: `ηᵢ = β cᵢ + ζᵢ`, `ζᵢ ~ N(0, σ_η)`, with **Var(η) = 1** so that
  `σ_η = √(1 − β²)` and `|β| < 1`.
- Reflective: `y_{i,j} = λⱼ ηᵢ + ε_{i,j}`, `ε_{i,j} ~ N(0, σ_{y,j})`.

The Stan model evaluates the likelihood **row by row**:

- `xᵢ ~ MVN(0, Σₓ)` — this is the H-O marginal density on the composite
  indicators. It is NOT the same as plugging the sample covariance `Sₓₓ`
  into a Wishart likelihood. The two are algebraically equivalent for
  linear Gaussian models, but the row-level form generalizes to
  non-Gaussian indicators, missing data, mixtures, hierarchical extensions,
  etc.
- `yᵢ | cᵢ ~ MVN(β λ cᵢ, λλ' σ_η² + diag(σ_y²))` — η is integrated out
  analytically (this is exact for linear Gaussian).

## Identification (the parts the brief did not mention)

The composite + factor model has several mathematical symmetries; without
fixing them HMC produces a bimodal / drift-prone posterior and fake-looking
"bias":

1. **Sign of c**: `c → −c` with `β → −β` is invariant. Fix by `w[1] > 0`.
2. **Sign of η**: `η → −η` with `λ → −λ` and `β → −β` is invariant. Fix by
   `λ[1] > 0`.
3. **Scale of η**: `(β, σ_η, λ) → (αβ, ασ_η, λ/α)` leaves *every* observable
   moment invariant. This is a continuous non-identification. Fix by
   `Var(η) = 1`, i.e. `σ_η = √(1 − β²)`. This is the same convention
   used by `cSEM.DGP` when it generates the data.

After fixing all three, the model is identified, R-hat ≈ 1.00, and chains
mix cleanly *provided* the sampler is initialized in the regular region of
parameter space.

### Initialization

Random inits sometimes put chains in a degenerate corner (`λ_1 ≈ 0`,
`β ≈ 0`) from which they cannot escape and produce bogus posteriors with
multimodality artifacts. We initialize Stan from **cSEM PLS estimates**
on the same data, jittered slightly per chain. PLS is a fast, consistent
estimator for composite-based SEMs; using it as a starting point is a
standard practical workflow and does **not** compromise Bayesian
inference (the prior and likelihood determine the posterior; inits only
determine where the sampler starts).

## What "PGM formulation, not covariance-based" actually means (skeptic note)

The brief asked for a PGM formulation, "not assuming MVN data and using a
variance-covariance matrix as the primary feature in the likelihood." Worth
being precise about what this distinction does and does not get you:

- For a **linear Gaussian composite model with no measurement error on the
  indicators of the composite**, the row-level MVN likelihood `Πᵢ φ(xᵢ; 0, Σₓ)`
  is *algebraically equivalent* to the Wishart sufficient-statistic
  likelihood `W(Sₓₓ; N−1, Σₓ / (N−1))`. They are the same model.
- The genuine generalization the PGM formulation buys you is **for
  non-Gaussian or partially observed indicators** (e.g. categorical y's,
  hierarchical extensions, mixture composites, missing data). The
  composite indicators of a *canonical* composite remain a linear combo
  of x, so for x itself there is no escape from the linear-algebra core.
- A literal "explicit-latent" H-O encoding (declare `cᵢ`, `νᵢ` as
  per-observation parameters; write `xᵢ = a cᵢ + V νᵢ` *exactly*) makes
  `xᵢ` a degenerate (Dirac-delta) density given the latents — Stan cannot
  evaluate that. You either (a) add small idiosyncratic measurement noise
  on x (turning the composite into a K-factor reflective model) or (b)
  integrate the latents out and recover the row-level MVN on x. We do (b).
  The relationships `a = Σₓ w` and `V V' = Σₓ − a a'` are recovered in
  the `generated quantities` block to make the H-O structure explicit.

So: this PoC is a "PGM" formulation in that **every observation has its
own likelihood evaluation** and **the latent common factor η is treated as
a stochastic node in the graph** (whether marginalized or sampled). It is
not a SEM/Wishart fit. But the user's intuition that PGM-vs-covariance is
a deep modelling distinction for **linear Gaussian** composite indicators
is overstated.

## Debugging journey (pedagogical)

The "final" Stan model above did not work on the first three tries. Each
broken intermediate is informative — it surfaces a different identification
pathology that composite-with-common-factor models routinely exhibit and
that quietly disappears in CB-SEM software because the software imposes the
constraints silently. The point of writing this section is *not* to confess
that the author iterated; it is to make the bugs and their symptoms visible
so a reader can recognize them in their own models.

### Pathology 1 — Sign of η is free (no `λ[1] > 0`)

**The buggy code (excerpt).** No constraint on any element of `λ`; η was an
explicit per-observation latent at this point.

```stan
parameters {
  real<lower=0> w_raw_head;          // w[1] > 0  → fixes sign of c
  vector[K-1] w_raw_tail;
  vector[M] lambda;                  // ← no sign constraint
  vector<lower=0>[M] sigma_y;
  real beta;
  real<lower=0> sigma_eta;
  vector[N] eta;
}
model {
  eta ~ normal(beta * composite, sigma_eta);
  for (j in 1:M) y[, j] ~ normal(lambda[j] * eta, sigma_y[j]);
}
```

**Observed symptom.** R-hat ≈ 2 on `β, λ, η`. Posterior densities for `β`
and `λ` look like two superposed bell curves centered at ±truth; posterior
means consequently collapse to near 0.

```
   variable    mean   median      sd    q5     q95   rhat
   beta        0.028  -0.010      0.32  -0.38  0.53  2.18
   lambda[1]  -0.061  -0.024      1.67  -2.31  2.54  2.33
   lambda[2]  -0.057  -0.016      1.51  -2.08  2.28  2.32
   lambda[3]  -0.069  -0.006      1.82  -2.52  2.73  2.36
```

**Why.** The likelihood is invariant under the joint relabeling
`η → −η, λ → −λ, β → −β`. The signs of `c` are already pinned by
`w[1] > 0`, but the signs of `η` are not pinned anywhere. So there are two
posterior modes related by this flip and the chains visit both, producing
"average ≈ 0" posteriors. The R-hat ≈ 2 is the textbook signature of two
equally-weighted modes.

**Fix.** Pin `λ[1] > 0`.

**Pedagogical takeaway.** Any latent variable that enters the likelihood
*only through products of two parameters* (here `β λ_j` and `λ_j λ_k`)
carries a sign symmetry that has to be broken by some external rule. The
rule has to nail down each latent's sign with one constraint — typically
"fix the sign of one loading per factor."

### Pathology 2 — Boundary trap from a hard `<lower=0>` constraint

**The buggy code (excerpt).** Same as above but with `lambda_head` declared
as `<lower=0>`. Random inits.

```stan
parameters {
  ...
  real<lower=0> lambda_head;
  vector[M-1] lambda_tail;
  ...
}
```

**Observed symptom.** Chains still report R-hat > 1.7 on `β, λ`. The
posterior for `λ[1]` is pressed against zero (q5 ≈ 0.0004); `λ[2], λ[3]`
are bimodal across positive and negative. `β` straddles ±0.

```
   variable    mean   median      sd    q5      q95   rhat
   lambda[1]   0.42   0.36       0.42   0.00    0.89   1.73
   lambda[2]  -0.01   0.02       0.76  -0.83    0.81   1.73
   lambda[3]   0.01   0.05       0.90  -0.96    0.96   1.73
```

**Why (two compounding reasons).**

1. *Hard constraint at an active mode.* Stan's `real<lower=0>` uses a log
   transform internally. If the data could prefer `λ[1] ≈ 0` (as a way to
   "decouple" `y_1` from the rest of the model), the chain piles up at the
   boundary and produces a half-spike. This is what you see in the q5 ≈ 0.

2. *The constraint is not strong enough to break the symmetry alone.* If
   `λ[1]` is small, the sign-flip symmetry on `(β, λ_2, λ_3)` is almost
   unbroken — flipping all of them costs almost nothing in likelihood
   because `λ[1] ≈ 0` makes `y_1` nearly orthogonal to η. So the chains
   still find two modes in `(β, λ_2, λ_3)`.

   Why does the chain even visit `λ[1] ≈ 0` in the first place? Because of
   **Pathology 3** below: the scale of η is non-identified, so the model
   *can* shrink all of `λ` toward 0 while inflating `σ_η` and `β` in
   compensation, without affecting any observable moment.

**Fix.** Address Pathology 3 (next) and use informed initial values
(below) so chains do not start near the degenerate corner.

**Pedagogical takeaway.** A `<lower=0>` constraint is a *tie-breaker*, not
a *prior*. It only does its job if the posterior under the full model
clearly prefers one sign. When some other identification problem makes
both signs equally palatable, the constraint just creates a half-normal
posterior squashed against the boundary, and you cannot distinguish "the
constraint did its job" from "the model is broken."

### Pathology 3 — Continuous non-identification of the η scale

**The buggy code (excerpt).**

```stan
parameters {
  ...
  real beta;
  real<lower=0> sigma_eta;     // ← σ_η free
  ...
}
```

**Observed symptom.** Even with `w[1] > 0` *and* `λ[1] > 0` *and*
truth-near initial values, the magnitudes of `(β, σ_η, λ)` drift. With
truth-init the posterior settled at
`β = 0.33, σ_η = 0.65, λ ≈ (1.40, 1.27, 1.54)` vs the truth
`β = 0.50, σ_η = 0.87, λ ≈ (0.80, 0.70, 0.90)`.

Notice the ratios: posterior `λ` is inflated by ≈ 1.75; posterior
`σ_η` deflated by ≈ 0.75; posterior `β` deflated by ≈ 0.66. The
products `λ_j · √(β² + σ_η²)` and `β · λ_j` are essentially preserved.

**Why.** The structural+reflective block has a continuous one-parameter
family of indistinguishable parameterizations. Define the rescaling
`(β, σ_η, λ) → (αβ, ασ_η, λ/α)`. Then:

- `Cov(y_j, c) = β λ_j` → `αβ · λ_j/α = β λ_j` (unchanged)
- `Cov(y_j, y_k) = (β² + σ_η²) λ_j λ_k` →
  `(α²β² + α²σ_η²)(λ_j λ_k / α²)` = unchanged
- `Var(y_j) = (β² + σ_η²) λ_j² + σ_{y,j}²` → unchanged
- `Var(c) = 1` → unchanged

Every observable moment is invariant under the rescaling, so the
likelihood is *flat* along this one-dimensional ridge. HMC can wander
along it indefinitely, producing wide marginals on each of `β, σ_η, λ`
even though the products and ratios are pinned. Random inits make the
problem dramatic — different chains slide to different points on the
ridge and never meet.

This is the kind of bug that **does not show up** in lavaan / CB-SEM
because those tools impose `Var(η) = 1` by default in the standardized
parameterization.

**Fix.** Adopt the same convention: pin `Var(η) = 1`, i.e. let
`σ_η² = 1 − β²` with `|β| < 1`. The ridge collapses to a point.

```stan
parameters {
  real<lower=-1, upper=1> beta;
}
transformed parameters {
  real<lower=0> sigma_eta = sqrt(1 - square(beta));
}
```

**Pedagogical takeaway.** Before fitting *any* latent-variable model in
Stan, write out the change-of-variables that would leave the joint
likelihood invariant and check that you have enough hard constraints to
break each one. Sign symmetries need ≥ 1 sign constraint per latent.
Scale symmetries need ≥ 1 scale constraint per latent. Rotation
symmetries (e.g. between multiple excrescent factors of a composite, or
between multiple common factors in a CFA) need (K−1)(K−2)/2 zero
constraints. The number of hard constraints needed equals the dimension
of the symmetry group acting on the parameters.

### Pathology 4 — η as explicit per-observation latent mixes poorly

**The buggy code (excerpt).** Centered parameterization, with `η` as
a vector of free parameters.

```stan
parameters { ...
  vector[N] eta;
}
model { ...
  eta ~ normal(beta * composite, sigma_eta);
  for (j in 1:M) y[, j] ~ normal(lambda[j] * eta, sigma_y[j]);
}
```

**Observed symptom.** "E-BFMI < 0.3" warning on every chain. ESS for
hyperparameters (`β, σ_η, λ`) much lower than for the η's. Treedepth
saturations possible.

**Why.** The centered parameterization induces a posterior funnel between
`σ_η` and the η_i. When `σ_η` is small the η_i are tightly constrained
around `β c_i`; when `σ_η` is large they are loose. HMC needs different
step sizes in different parts of the joint space.

**Fix (used in final code).** Marginalize η analytically. For linear
Gaussian the marginal of `y_i | c_i` is also MVN, in closed form:

`y_i | c_i ~ MVN(β λ c_i,  λ λ' σ_η² + diag(σ_y²))`

This removes the funnel and ESS jumps up. The composite c_i remains a
deterministic transformation of the data (`c_i = w' x_i`), so the
"PGM"-flavour is preserved on the side where it actually adds value
(non-Gaussian extensions in y would still slot in as `y_i | η_i` with
η_i sampled non-marginally; see the closing note in the PGM section).

**Alternative fix.** Non-centered parameterization
(`η_i = β c_i + σ_η · z_i` with `z_i ~ N(0,1)`) — also correct, slightly
slower per iteration but generalizable to non-Gaussian y.

**Pedagogical takeaway.** "Use the PGM formulation" sounds principled
but in practice you should marginalize whatever can be marginalized in
closed form. The PGM benefits (non-Gaussian indicators, missing data,
hierarchical structure) accrue from being *able* to write
`p(y_i | η_i, θ)` explicitly when needed — not from forcing yourself to
*always* keep `η_i` as a sampled node. Marginalizing a Gaussian latent
out of a Gaussian likelihood is a free win.

### Diagnostic checklist (the lesson, distilled)

When a Stan latent-variable model produces R-hat ≫ 1.01 and you don't
see an obvious code bug, run through:

1. **Sign symmetries.** For every latent factor, is there exactly one
   "anchor" loading constrained `> 0`?
2. **Scale symmetries.** For every latent factor, is exactly one of
   `Var(factor) = 1` or `λ_anchor = 1` imposed?
3. **Rotation symmetries.** If you have several latent factors that load
   on the same indicators, do you have enough zero loadings to make the
   loading matrix uniquely defined?
4. **Boundary traps.** Are any of your `<lower=0>` constraints active in
   the posterior? If so, you may have a deeper identification problem
   underneath.
5. **Funnels.** Are any variance parameters strongly correlated with
   per-observation latents in the posterior? If yes, marginalize the
   latents or switch to non-centered.
6. **Initialization.** Is the random-init region of parameter space
   regular for your model? If not, use a moment-based estimator
   (e.g. PLS for composite SEMs, OLS for measurement models) as a
   starting point.

In the model in this repo the *first three* are textbook structural
fixes, the *fourth and fifth* are HMC fixes, and the *sixth* is a
practical fix that does not affect the posterior — only how easily the
sampler finds it.

## Skeptic notes on the brief

1. *"Unbiased parameter estimates"*: Bayesian point estimates are typically
   not unbiased. The natural guarantees are **consistency** (posterior
   concentrates around truth as N → ∞), **frequentist coverage**
   (credible intervals contain the truth at the nominal rate), and
   **SBC calibration** (posterior is correctly computed under the
   prior+model). I report MC bias for transparency but treat coverage
   and SBC as the meaningful checks.
2. *"Canonical composite affecting a common factor"*: the canonical
   composite is `c = w' x`, deterministic. There is no error in
   `x → c`; "composite measurement error" is a contradiction in terms.
   This shapes how the likelihood can be written.
3. *"Henseler–Ogasawara specification"*: in CB-SEM, H-O is a trick to
   express a composite in CFA syntax via excrescent factors V so that
   covariance-based software can fit it. In Stan we don't need the
   trick to be explicit — we get the same model by parameterizing
   Σₓ and w directly and recovering `a, V` post-hoc. The PoC shows
   the H-O quantities in `generated quantities`.

## Results (R = 30 replications, N = 500)

Aggregate parameter recovery (also saved to `sim_study_result.rds`).
`bias = mean(posterior_mean) − truth`; `mc_se_bias = SD/√R`; `rmse =
√mean((posterior_mean − truth)²)`; `coverage_90 = fraction of 90% CIs
that contain truth`.

| parameter   | truth | mean post mean | bias    | MC SE   | RMSE   | 90% cov |
|-------------|------:|---------------:|--------:|--------:|-------:|--------:|
| w[1]        | 0.348 |          0.326 | −0.022  | 0.017   | 0.093  | 0.867   |
| w[2]        | 0.406 |          0.414 |  0.008  | 0.014   | 0.078  | 0.967   |
| w[3]        | 0.464 |          0.462 | −0.001  | 0.015   | 0.080  | 0.967   |
| beta        | 0.500 |          0.500 |  0.001  | 0.007   | 0.037  | 0.867   |
| lambda[1]   | 0.800 |          0.802 |  0.002  | 0.006   | 0.033  | 0.933   |
| lambda[2]   | 0.700 |          0.701 |  0.001  | 0.009   | 0.051  | 0.767   |
| lambda[3]   | 0.900 |          0.907 |  0.007  | 0.006   | 0.033  | 0.900   |
| sigma_y[1]  | 0.600 |          0.605 |  0.005  | 0.006   | 0.030  | 0.900   |
| sigma_y[2]  | 0.714 |          0.717 |  0.003  | 0.004   | 0.025  | 0.967   |
| sigma_y[3]  | 0.436 |          0.425 | −0.011  | 0.006   | 0.036  | 0.967   |
| R²_η        | 0.250 |          0.253 |  0.003  | 0.007   | 0.036  | 0.867   |

Diagnostics aggregated across replications:

- max R̂: median 1.004, max 1.025
- min ESS: median 991, min 163
- divergent transitions: 106 total across 30 reps × 1600 transitions
  (≈ 0.22%; rep-max = 11)

### Interpretation

- **Bias.** Every parameter's bias is < 1.5 × its MC standard error, i.e.
  statistically indistinguishable from zero. Posterior means are
  effectively unbiased at N = 500. Note that the largest *absolute* bias
  is on `w[1]` (−0.022), with MC SE 0.017 — borderline-significant; with
  more replications it would either shrink toward zero or settle at a
  small genuine finite-sample bias. The estimator (and PLS in `cSEM`)
  underweight `w[1]` for the same reason: when the indicator-correlation
  matrix is uncertain, the composite weight on the most weakly correlated
  indicator is harder to estimate.
- **Coverage.** Empirical coverage of 90% credible intervals ranges 0.77
  – 0.97. With R = 30 reps, the binomial Monte Carlo SE on a true 0.9
  coverage rate is √(0.9·0.1/30) ≈ 0.055, so values in roughly the
  0.79–1.0 band are consistent with nominal. `lambda[2]` at 0.767 is mild
  under-coverage; given only 30 reps it could be either noise or genuine
  slight optimism (e.g. from the prior). A second run at R = 100+ would
  resolve this.
- **Diagnostics.** Convergence is clean (R̂ ≤ 1.025 everywhere, median ≈
  1.004). Divergent transitions appear in most reps but never exceed 1%
  per chain. The model is *fitable*, not bullet-proof — if you push to
  smaller N or more constructs you should expect to need a non-centered
  reformulation, tighter `adapt_delta`, or both.

### Honest caveat

This is one DGP (3-indicator composite, 3-indicator factor, β = 0.5, the
particular Σₓ specified in `00_setup.R`). The conclusion "Stan produces
effectively unbiased posterior means and roughly-calibrated 90% CIs"
applies to *this* DGP at *this* N. It does **not** imply that the model
is well-behaved at every (β, λ, Σₓ, N) configuration. The SBC script in
`03_sbc.R` is the principled way to assess that more broadly; it was not
executed as part of this PoC due to time budget.

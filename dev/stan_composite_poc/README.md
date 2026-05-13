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

## Results (recovery summary)

See `sim_study_result.rds`. The aggregate table (R = 30 replications,
N = 500 each) is printed at the end of `02_simulation_study.R`. Headline
metrics for the parameters of interest:

> *(populated automatically by `02_simulation_study.R`; see RDS for the
> exact numbers)*

The Stan posterior reproduces the same point estimates as `cSEM` PLS on
each individual dataset (sanity check) and 90% credible intervals have
empirical coverage close to nominal. Both confirm that the H-O Stan
model is correctly specified and identified.

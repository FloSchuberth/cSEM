# IGSCA Benchmarking Plan

## Objective
Create 3 additional benchmark functions to evaluate `csem()` bottlenecks under different N × p conditions. Profiling via `profvis::profvis()` will be done separately in `dev/igsca/benchmarking.Qmd`.

## Current State
- **`bigN_smallp`**: N=1000, 1 common factor (xi1, 3 indicators), 1 composite (xi2, 3 indicators) → 6 total indicators

## 3 New Functions

### 1. `smallN_smallp` — Small N, Small p
- **N = 50**, same model as current (1 factor + 1 composite = 6 indicators)
- Identical DGP/model to `bigN_smallp`, only `.N` changes

### 2. `smallN_bigp` — Small N, Big p
- **N = 50**, expanded model:
  - **10 common factors**: xi1 (original) + xi3, xi5, xi7, xi9, xi11, xi13, xi15, xi17, xi19 — each `=~` 3 unique indicators
  - **10 composites**: xi2 (original) + xi4, xi6, xi8, xi10, xi12, xi14, xi16, xi18, xi20 — each `<~` 3 unique indicators
  - **60 total indicators** (20 constructs × 3 indicators)
  - Structural paths: each factor predicts its paired composite

### 3. `bigN_bigp` — Big N, Big p
- **N = 1000**, same expanded model as `smallN_bigp`

## DGP Construction (Big p models)

For each new common factor (e.g., xi3):
```
xi3 =~ 0.6*x31 + 0.8*x32 + 0.7*x33
```

For each new composite (e.g., xi4):
```
xi4 <~ 0.4*x41 + 0.3*x42 + 0.2*x43
x41 ~~ 0.4*x42 + -0.3*x43
x42 ~~ 0.4*x43
```

Structural model — paired pattern:
- Each factor predicts its paired composite: `xi4 ~ 0.5*xi3`, `xi6 ~ 0.5*xi5`, etc.
- Mirrors the original `xi2 ~ 0.5*xi1` pattern

## Accuracy Considerations for Profiling
- Use `.empirical = TRUE` in `generateData()` so sample covariance matches population exactly — removes sampling variability
- Use `set.seed()` before each data generation for reproducibility
- Use identical csem arguments across all functions (`.tolerance`, `.iter_max`, `.conv_criterion`, `.GSCA_modes`, `.disattenuate`)
- Recommend `gc()` before each `profvis()` call in the Qmd to minimize GC pauses during profiling

## Naming Convention
| Function         | N    | p (indicators) |
|------------------|------|-----------------|
| `bigN_smallp`    | 1000 | 6               |
| `smallN_smallp`  | 50   | 6               |
| `smallN_bigp`    | 50   | 60              |
| `bigN_bigp`      | 1000 | 60              |

## File Output
- All 4 functions in `igsca_benchmark.R`
- Each function generates its own data internally and calls `csem()` with consistent arguments
- DGP strings and model strings defined at the top of the file as shared variables
- Original population weight calculations preserved in `if (FALSE) { ... }` block
- Indicator naming uses underscores: `x1_1, x1_2, x1_3` (construct 1, indicators 1-3)
- Profiling via `profvis::profvis()` is done in `benchmarking.Qmd`, not in this file

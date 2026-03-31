# Optimizing Step 3: Unique Scores (U) Update in IGSCA

## Problem Statement

The `updateUD()` function in `R/helper_igsca.R` implements Step 3 of the IGSCA
alternating least squares algorithm as described in Hwang et al. (2017). This step
updates the unique scores matrix U by solving an orthogonal Procrustes problem
constrained so that U is orthogonal to the construct scores Eta.

The current implementation is extremely computationally expensive for large sample
sizes (N), making IGSCA impractical for large datasets.

---

## Mathematical Formulation (from gscam_original.pdf, Step 3)

**Objective:** Maximize `tr(Ũ' Γ_⊥' Z D)` subject to `Ũ'Ũ = I_J`

Where:
- `Z` (N × J): standardized indicator data matrix
- `Γ` (equivalently `Eta_normed`, N × P): normalized construct scores
- `Γ_⊥` (N × (N-P)): orthonormal basis for the null space of Γ
- `D` (J × J): diagonal matrix of unique loadings
- `U` (N × J): unique scores matrix
- N = sample size, P = number of constructs, J = number of indicators

**Solution via SVD:**
1. QR decompose Γ = [Q₁, Q₂]R where Q₂ = Γ_⊥
2. Compute SVD(Γ_⊥' Z D) = P Δ T'
3. Set Ũ = P T' (orthogonal Procrustes solution)
4. Recover U = Γ_⊥ Ũ

---

## Why the Current Implementation Is Expensive

### Current code (R/helper_igsca.R, lines 541-548):

```r
Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
    (n_constructs + 1):n_case, drop = FALSE]
svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
Utilde <- svd_mx$v %*% t(svd_mx$u)
U <- Eta_Q2 %*% Utilde
```

### Bottleneck 1: Complete QR decomposition — O(N²P) time, O(N²) memory

`qr.Q(qr(Eta_normed), complete = TRUE)` computes the **full** N × N orthogonal
matrix Q. When Eta_normed is N × P with P << N, we only need P columns for the
column space, but `complete = TRUE` computes all N columns to extract the (N-P)
null space columns.

| N (cases) | P (constructs) | Q matrix size | Memory (doubles) |
|-----------|----------------|---------------|------------------|
| 1,000     | 5              | 1000 × 1000   | ~8 MB            |
| 10,000    | 5              | 10000 × 10000 | ~800 MB          |
| 50,000    | 5              | 50000 × 50000 | ~20 GB           |
| 100,000   | 5              | 100000×100000  | ~80 GB           |

### Bottleneck 2: Eta_Q2 storage — O(N(N-P)) ≈ O(N²) memory

`Eta_Q2` is the extracted N × (N-P) submatrix. Since P << N, this is nearly N × N.

### Bottleneck 3: Matrix products involving Eta_Q2 — O(N²J) computation

- `t(Z_normed) %*% Eta_Q2`: (J × N)(N × (N-P)) → O(JN(N-P)) ≈ O(JN²)
- `D %*% ...`: O(J(N-P)) (cheap, D is diagonal)
- `svd(...)`: SVD of J × (N-P) matrix → O(J(N-P)²) ≈ O(JN²)
- `Eta_Q2 %*% Utilde`: (N × (N-P))((N-P) × J) → O(NJ(N-P)) ≈ O(NJN) = O(N²J)

### Total complexity of current approach

| Resource   | Cost                  | For N=10000, J=20, P=5 |
|------------|-----------------------|------------------------|
| Memory     | O(N²)                | ~800 MB                |
| Computation| O(N²J + N²P + JN²)   | ~2 × 10⁹ FLOPs        |

This runs **inside the ALS loop** (typically 50-500 iterations), multiplying the
total cost enormously.

---

## Proposed Efficient Approach: Thin QR Projection

### Key Insight

We never need the explicit null space basis Γ_⊥. Instead, we can work with the
**orthogonal projector** P_⊥ = I - Q_thin Q_thin', where Q_thin is the N × P
thin QR factor of Eta_normed.

**Claim:** `SVD(P_⊥ Z D)` yields the same Procrustes solution as `SVD(Γ_⊥' Z D)`.

### Proof of Equivalence

1. **Projection identity:** Since Q_thin and Γ_⊥ together form a complete
   orthonormal basis for R^N:
   ```
   P_⊥ = I - Q_thin Q_thin' = Γ_⊥ Γ_⊥'
   ```

2. **SVD relationship:** Let M = Γ_⊥' Z D (the (N-P) × J matrix). Then:
   ```
   P_⊥ Z D = Γ_⊥ Γ_⊥' Z D = Γ_⊥ M
   ```
   Since Γ_⊥ has orthonormal columns (Γ_⊥' Γ_⊥ = I), premultiplying by Γ_⊥
   preserves singular values and right singular vectors.

3. **SVD correspondence:** If M = P Δ T', then:
   ```
   Γ_⊥ M = (Γ_⊥ P) Δ T'
   ```
   So `SVD(P_⊥ Z D)` has:
   - Left singular vectors: Γ_⊥ P (already in the null space of Eta)
   - Right singular vectors: T (identical to original)
   - Singular values: Δ (identical to original)

4. **Procrustes solution:**
   ```
   U = (left s.v.) × (right s.v.)' = (Γ_⊥ P) T' = Γ_⊥ (P T') = Γ_⊥ Ũ
   ```
   This is identical to the original algorithm.

5. **Constraint verification:**
   - `U' Eta_normed = 0` ✓ (columns of U lie in null space of Eta')
   - `U' U = I_J` ✓ (left and right singular vectors are orthonormal)

### Proposed code:

```r
# Efficient projection method: O(NJ² + NPJ) instead of O(N²J)
# Thin QR of Eta_normed: N × P (not N × N)
Q_thin <- qr.Q(qr(Eta_normed))
# Project: M_proj = (I - QQ')ZD = P_⊥ ZD, all N × J matrices
ZD <- Z_normed %*% D
M_proj <- ZD - Q_thin %*% crossprod(Q_thin, ZD)
# SVD of N × J matrix (J << N, so fast) and Procrustes solution
svd_mx <- svd(M_proj)
U <- svd_mx$u %*% t(svd_mx$v)
```

### Complexity comparison

| Resource    | Old (complete QR)     | New (thin QR projection) | Speedup ratio |
|-------------|-----------------------|--------------------------|---------------|
| Memory      | O(N²)                | O(NJ + NP)              | N / (J+P)     |
| Computation | O(N²J + N²P)         | O(NJ² + NPJ)            | N / J         |

**Concrete example (N=10000, J=20, P=5):**

| Resource    | Old           | New          | Speedup |
|-------------|---------------|--------------|---------|
| Memory      | ~800 MB       | ~1.6 MB      | 500×    |
| Computation | ~2×10⁹ FLOPs | ~4×10⁶ FLOPs | 500×    |

---

## Analysis of the Commented-Out Hwang 2021 Alternative (lines 502-527)

### What it does

The commented-out code from Hwang (private communication, 2021) uses a different
computational path to achieve the same mathematical result:

```r
svd_etaprod <- svd(t(Eta_normed) %*% Eta_normed)   # SVD of P×P matrix
gd2 <- diag(svd_etaprod$d)
gv <- svd_etaprod$v
GU <- Eta_normed %*% gv %*% solve(sqrt(gd2))       # Orthonormalized Eta
M3 <- Z_normed %*% D - GU %*% (t(GU) %*% Z_normed) %*% t(D)   # = P_⊥ Z D
```

This IS mathematically equivalent to our proposed approach — `M3 = P_⊥ Z D`
because `GU GU' = Eta (Eta'Eta)^{-1} Eta' = P_Eta`, so `I - GU GU' = P_⊥`.

### Why it produces different numerical results

1. **`solve(sqrt(gd2))` instability:** This inverts the square root of the
   eigenvalues of Eta'Eta. If Eta has near-collinear columns, these eigenvalues
   can be very small, and their inverse amplifies floating-point errors in GU.

2. **Condition number squaring in the `n_case > n_indicators` branch:**
   ```r
   svd_M3prod <- svd(t(M3) %*% M3)    # forms M3'M3
   u <- M3 %*% v %*% solve(sqrt(d2))   # another inversion
   ```
   Computing eigendecomposition via M3'M3 squares the condition number of M3.
   Combined with the already-amplified errors from step 1, this can cause
   significant numerical drift.

3. **Our proposed approach avoids both issues:**
   - Thin QR decomposition is backward-stable (no condition number issues)
   - Direct SVD of M_proj avoids condition number squaring
   - No matrix inversions required (the Procrustes solution u·v' uses only
     orthogonal factors, not singular values)

---

## Summary

The thin QR projection approach is:
1. **Mathematically equivalent** to the original Hwang et al. (2017) method
2. **~500× faster** and uses **~500× less memory** for typical problem sizes
3. **More numerically stable** than the Hwang 2021 alternative
4. **Simpler code** — fewer lines, no complete QR, no null space extraction

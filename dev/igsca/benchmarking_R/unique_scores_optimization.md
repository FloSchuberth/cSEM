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
2. Compute SVD(D Z' Γ_⊥) = u Δ v' (equivalently, SVD of Γ_⊥' Z D transposed)
3. Set Ũ = v u' (orthogonal Procrustes solution)
4. Recover U = Γ_⊥ Ũ

---

## Why the Original Implementation Is Expensive

### Original code:

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
- `svd(...)`: SVD of J × (N-P) matrix → O(J(N-P)²) ≈ O(JN²)
- `Eta_Q2 %*% Utilde`: (N × (N-P))((N-P) × J) → O(NJ(N-P)) ≈ O(N²J)

### Total complexity of original approach

| Resource   | Cost                  | For N=10000, J=20, P=5 |
|------------|-----------------------|------------------------|
| Memory     | O(N²)                | ~800 MB                |
| Computation| O(N²J + N²P + JN²)   | ~2 × 10⁹ FLOPs        |

This runs **inside the ALS loop** (typically 50-500 iterations), multiplying the
total cost enormously.

---

## Implemented Approach: Implicit Householder

### Key Insight

R's `qr()` stores the Householder reflections that define Q without ever forming
the N × N matrix. The functions `qr.qty()` and `qr.qy()` apply Q' and Q to
arbitrary matrices using these reflections — in O(NPJ) time and O(NJ) memory.

This lets us compute Q₂' Z and Q₂ Ũ without forming Q₂ explicitly, while
**structurally guaranteeing** that U = Q₂ Ũ lies in the null space of Eta.

### Why the projection approach (SVD of P_⊥ Z D) fails

An earlier attempt used `SVD(P_⊥ Z D)` where `P_⊥ = I - Q Q'` is computed via
thin QR. While mathematically equivalent, this approach has a **structural flaw**
for rank-deficient cases (e.g., single-indicator common factors):

- When indicator j is the sole indicator of construct k, then Z[,j] ∝ Eta[,k],
  so the j-th column of `P_⊥ Z D` is approximately zero.
- The SVD of M_proj (N × J) assigns an **arbitrary** left singular vector to
  the near-zero singular value. This vector must be orthogonal to the other
  J-1 left singular vectors, but need NOT lie in null(Eta').
- This arbitrary vector contaminates U[,j] via the Procrustes product `u v'`,
  causing U[,j]' Eta ≠ 0 and D[j,j] ≠ 0, which makes the loading c_j < 1
  instead of the correct c_j = 1.

In contrast, the original code computes U = Γ_⊥ Ũ, which **structurally
guarantees** that every column of U is in null(Eta'), regardless of the SVD
behavior on rank-deficient inputs. The implicit Householder approach preserves
this structural guarantee.

### Implemented code:

```r
qr_eta <- qr(Eta_normed)
# Q2' Z via implicit Householder: O(NPJ), result is N×J, extract rows (P+1):N
QtZ_null <- qr.qty(qr_eta, Z_normed)[(n_constructs + 1):n_case, , drop = FALSE]
# SVD of D Z' Q2 = D t(Q2' Z): J × (N-P), same input as original approach
svd_mx <- svd(D %*% t(QtZ_null))
Utilde <- svd_mx$v %*% t(svd_mx$u)  # (N-P) × J
# Recover U = Q2 Utilde via implicit Householder: O(NPJ)
U <- qr.qy(qr_eta, rbind(matrix(0, n_constructs, n_indicators), Utilde))
```

### How it works

1. **`qr(Eta_normed)`** — stores Householder reflections defining the full
   Q = [Q₁ | Q₂], without forming Q. Memory: O(NP).

2. **`qr.qty(qr_eta, Z_normed)`** — computes Q' Z = [Q₁'Z; Q₂'Z] by applying
   Householder reflections to Z. Cost: O(NPJ). Result: N × J.
   We extract rows (P+1):N to get Q₂'Z, which is (N-P) × J.

3. **`svd(D %*% t(QtZ_null))`** — SVD of D Z' Q₂ (J × (N-P)). This is the
   **exact same SVD input** as the original code, so the Procrustes solution
   Ũ = v u' is numerically equivalent.

4. **`qr.qy(qr_eta, [0; Utilde])`** — computes Q [0; Ũ] = Q₂ Ũ by applying
   Householder reflections. The zero block in the first P rows ensures only Q₂
   (the null space basis) acts on Ũ. Cost: O(NPJ). Result: N × J.

### Complexity comparison

| Resource    | Original (complete QR) | Implicit Householder   | Speedup ratio |
|-------------|------------------------|------------------------|---------------|
| Memory      | O(N²)                 | O(NJ + NP)             | N / (J+P)     |
| Computation | O(N²J + N²P)          | O(NJ² + NPJ)           | N / J         |

**Concrete example (N=10000, J=20, P=5):**

| Resource    | Original      | Implicit Householder | Speedup |
|-------------|---------------|----------------------|---------|
| Memory      | ~800 MB       | ~1.6 MB              | 500×    |
| Computation | ~2×10⁹ FLOPs | ~4×10⁶ FLOPs         | 500×    |

### Properties preserved

- **Same SVD input** as original → numerically equivalent Procrustes solution
- **Structural null-space guarantee** → U = Q₂ Ũ, so U' Eta = 0 always holds
- **Single-indicator constructs** → D[j,j] ≈ 0, loading c_j = 1 (verified)
- **Orthonormality** → U'U = Ũ'Q₂'Q₂Ũ = Ũ'Ũ = I_J

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

This IS mathematically equivalent — `M3 = P_⊥ Z D` because
`GU GU' = Eta (Eta'Eta)^{-1} Eta' = P_Eta`, so `I - GU GU' = P_⊥`.

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

3. **Same structural flaw as the thin QR projection approach:** It computes
   SVD(M3) where M3 = P_⊥ Z D, so it suffers from the same issue with
   arbitrary singular vectors escaping null(Eta') for rank-deficient cases.

---

## Summary

The implicit Householder approach is:
1. **Numerically equivalent** to the original Hwang et al. (2017) method
   (same SVD input, structural null-space guarantee)
2. **~500× faster** and uses **~500× less memory** for typical problem sizes
3. **Correctly handles single-indicator constructs** (loading = 1)
4. **More numerically stable** than the Hwang 2021 projection alternative

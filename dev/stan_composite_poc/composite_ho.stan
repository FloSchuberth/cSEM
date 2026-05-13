// ===========================================================================
// Henseler-Ogasawara (H-O) composite + common factor in Stan.
// Probabilistic graphical model formulation (row-level likelihoods,
// NOT a sample-covariance / Wishart likelihood).
//
// MODEL
//   Composite indicators x_i ∈ R^K, reflective indicators y_i ∈ R^M.
//   Composite (canonical, deterministic): c_i = w' x_i, with w'Σ_x w = 1.
//   H-O parameterization (latents):  x_i = a c_i + V ν_i,  a = Σ_x w.
//   Structural:                       η_i = β c_i + ζ_i,  ζ_i ~ N(0, σ_η).
//   Reflective:                       y_{i,j} = λ_j η_i + ε_{i,j}.
//
//   Indicators x_i: row-level MVN(0, Σ_x) — algebraically equivalent to
//     marginalizing (c_i, ν_i) under H-O.
//   Reflective y_i: η_i is marginalized analytically -- given c_i,
//     y_i ~ MVN(β λ c_i, λ λ' σ_η^2 + diag(σ_y^2)).
//   This is still a row-level likelihood (each i contributes independently).
//
// IDENTIFICATION
//   * w_raw[1] > 0  → fixes sign of c (else c → −c, β → −β symmetry).
//   * λ[1]  > 0     → fixes sign of η (else η → −η, λ → −λ, β → −β symmetry
//                     leaves the model invariant when c is signed).
//   * Var(c) = 1 set by post-hoc normalization w = w_raw / sqrt(w_raw' Σ_x w_raw).
//   * Var(η) = 1 imposed by setting σ_η = sqrt(1 − β^2), |β| < 1. This is the
//     cSEM/lavaan convention. WITHOUT this constraint the model has a
//     continuous non-identification: (β, σ_η, λ) → (αβ, ασ_η, λ/α) leaves
//     every observable moment invariant.
//   * Σ_x parameterized as a correlation matrix via LKJ-Cholesky; indicators
//     assumed unit-variance (the standardization used by cSEM.DGP).
//
// WHY MARGINALIZE η?
//   With η_i as explicit latents the posterior is mildly hierarchical and a
//   centered parameterization mixes poorly when σ_η is sizeable (and the
//   marginal η_i are highly correlated with the (β,λ) signs in a way that
//   produces multimodality before sign constraints fully take effect).
//   Marginalizing η is a closed-form, exact operation for linear Gaussian
//   models. The composite indicators x and the H-O latent structure (a, V,
//   ν, c) are *not* marginalized in any non-trivial way: c_i = w' x_i is
//   recovered exactly as a per-observation transformed parameter.
// ===========================================================================

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> M;
  matrix[N, K] x;
  matrix[N, M] y;
}

transformed data {
  vector[K] zeros_K = rep_vector(0.0, K);
}

parameters {
  real<lower=0> w_raw_head;
  vector[K - 1] w_raw_tail;
  cholesky_factor_corr[K] L_corr;

  real<lower=0> lambda_head;
  vector[M - 1] lambda_tail;
  vector<lower=0>[M] sigma_y;

  real<lower=-1, upper=1> beta;    // Var(c)=Var(η)=1 ⇒ |β| < 1
}

transformed parameters {
  vector[K] w_raw;
  vector[K] w;
  vector[M] lambda;
  real<lower=0> sigma_eta = sqrt(1 - square(beta));   // Var(η) = 1
  matrix[K, K] Sigma_x;
  vector[N] composite;

  w_raw[1] = w_raw_head;
  for (k in 2:K) w_raw[k] = w_raw_tail[k - 1];

  lambda[1] = lambda_head;
  for (j in 2:M) lambda[j] = lambda_tail[j - 1];

  Sigma_x = multiply_lower_tri_self_transpose(L_corr);

  {
    real wSw = quad_form(Sigma_x, w_raw);
    w = w_raw / sqrt(wSw);
  }
  composite = x * w;
}

model {
  // ---- Priors ----
  w_raw_head ~ normal(0, 1);
  w_raw_tail ~ normal(0, 1);
  L_corr     ~ lkj_corr_cholesky(2.0);

  lambda_head ~ normal(0, 2);
  lambda_tail ~ normal(0, 2);
  sigma_y     ~ exponential(1);

  // beta has support (-1, 1); uniform prior over this range
  // (sigma_eta is a deterministic function of beta, so no separate prior)

  // ---- Likelihood: composite indicators (row-level MVN, H-O marginal) ----
  for (n in 1:N)
    x[n] ~ multi_normal_cholesky(zeros_K, L_corr);

  // ---- Likelihood: reflective indicators given composite ----
  // y_i | c_i ~ MVN(β λ c_i, Σ_y), where Σ_y = λ λ' σ_η^2 + diag(σ_y^2).
  // Build Cholesky of Σ_y once; it does not depend on i.
  {
    matrix[M, M] Sigma_yy =
        diag_matrix(square(sigma_y)) + square(sigma_eta) * (lambda * lambda');
    matrix[M, M] L_yy = cholesky_decompose(Sigma_yy);
    matrix[N, M] mu_y;
    for (n in 1:N)
      mu_y[n] = (beta * composite[n]) * lambda';
    for (n in 1:N)
      y[n] ~ multi_normal_cholesky(mu_y[n]', L_yy);
  }
}

generated quantities {
  // H-O composite loadings (a = Σ_x w; satisfies w' a = 1, w' V = 0)
  vector[K] a_loading = Sigma_x * w;
  // Excrescent residual covariance V V' = Σ_x − a a'
  matrix[K, K] VVt = Sigma_x - a_loading * a_loading';
  // R^2 of structural equation
  real R2_eta = square(beta) / (square(beta) + square(sigma_eta));
}

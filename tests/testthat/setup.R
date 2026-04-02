if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  # On certain computers using BLAS, all available threads will be used for large computations, instead of only a few. This helps mitigate this issue.
  #   See this github issue for an example in the context of lme4 https://github.com/lme4/lme4/issues/492
  RhpcBLASctl::blas_set_num_threads(2)
}
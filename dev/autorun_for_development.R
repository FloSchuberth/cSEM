Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(1)

# Pretty error happening
library(rlang)
rlang::global_handle()

# rlang::last_trace()
# Could always reread https://adv-r.hadley.nz/debugging.html#debugging and https://rstudio.github.io/r-manuals/r-exts/Debugging.html#debugging-r-code
library(lobstr)

# lobstr::tree()


library(cli)
cli::pretty_print_code()

# Test by typing `lm` in the console

#' library(tibble)
# For printing something easier to see
#' tibble::as_tibble()

#' devtools::test_active_file()

# How to see source code of method
# getAnywhere('t.test.default')


# For replacing a function quickly 
# library(testthat)
# local_mocked_bindings(..., .package = NULL, .env = caller_env())
# https://testthat.r-lib.org/reference/local_mocked_bindings.html

library(testthat)
library(withr)
# withr::deferred_run()
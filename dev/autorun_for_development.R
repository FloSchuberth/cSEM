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

# ---- pkgnet: function network as a large, standalone picture ----
# pkgnet graphs the INSTALLED package, never the source tree. `pkg_path` is only
# used for the covr coverage overlay, so reinstall first or the network still
# shows functions you already deleted:
#   R CMD INSTALL --no-multiarch --no-docs .
#
# The report embeds the widget at a hardcoded height:480px, unreadable at cSEM's
# ~237 nodes. Build the widget yourself and resize it before saving.
if (!interactive()) {
  devtools::install(reload = FALSE)
  r <- pkgnet::FunctionReporter$new()
  r$set_package("cSEM", getwd())
  g <- r$graph_viz
  g$height <- "1200px"
  g$width <- "100%"
  # Optional: let dragged nodes push their neighbours around. Must be set at the
  # node level -- visPhysics(enabled = TRUE) alone leaves visIgraphLayout's
  # physics = FALSE on each node, and that wins. Slow to stabilise at ~237 nodes,
  # and the simulation discards igraph's layout, so leave it off for reading.
  # g <- visNetwork::visNodes(g, physics = TRUE)
  g <- visNetwork::visExport(g, type = "png", name = "cSEM_function_network")
  out <- here::here("dev/igsca/function_network.html")
  # selfcontained = TRUE genuinely inlines everything, but saveWidget creates
  # "<name>_files" next to `out` and then unlinks it relative to getwd(), so the
  # scaffolding is orphaned whenever those differ. Saving from the target
  # directory makes the cleanup land. The html is identical either way.
  withr::with_dir(dirname(out), {
    visNetwork::visSave(g, file = basename(out), selfcontained = TRUE)
  })
}

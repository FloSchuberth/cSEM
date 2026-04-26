# Pretty error happening
library(rlang)
rlang::global_handle()

# rlang::last_trace()
# Could always reread https://adv-r.hadley.nz/debugging.html#debugging and https://rstudio.github.io/r-manuals/r-exts/Debugging.html#debugging-r-code
library(lobstr)

# lobstr::tree()


# library(cli)
cli::pretty_print_code()

# Test by typing `lm` in the console

#' library(tibble)
# For printing something easier to see
#' tibble::as_tibble()
source("test-main.R")
# source("tests/testthat/test-main.R")

# ==============================================================================
# Tests 
# ==============================================================================

summarize(res_single_linear)
summarize(res_single_nonlinear)
summarize(res_single_2ndorder)
summarize(res_single_nonlinear_2ndorder)

summarize(res_single_linear_boot)
summarize(res_single_linear_boot)
summarize(res_single_2ndorder_boot)
summarize(res_single_nonlinear_2ndorder_boot)


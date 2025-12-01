source("test-main.R")
# 
# # ==============================================================================
# # Tests 
# # ==============================================================================

infer(res_single_linear_boot)

infer(res_single_linear_boot, .quantity = "mean")
infer(res_single_linear_boot, .quantity = "sd")
infer(res_single_linear_boot, .quantity = "bias")
infer(res_single_linear_boot, .quantity = "CI_standard_z")
infer(res_single_linear_boot, .quantity = "CI_standard_t")
infer(res_single_linear_boot, .quantity = "CI_percentile")
infer(res_single_linear_boot, .quantity = "CI_basic")
infer(res_single_linear_boot, .quantity = "CI_bc")
infer(res_single_linear_boot, .quantity = "CI_bca")

infer(res_single_2ndorder_boot, .quantity = "mean")
infer(res_single_2ndorder_boot, .quantity = "sd")
infer(res_single_2ndorder_boot, .quantity = "bias")
infer(res_single_2ndorder_boot, .quantity = "CI_standard_z")
infer(res_single_2ndorder_boot, .quantity = "CI_standard_t")
infer(res_single_2ndorder_boot, .quantity = "CI_percentile")
infer(res_single_2ndorder_boot, .quantity = "CI_basic")
infer(res_single_2ndorder_boot, .quantity = "CI_bc")
infer(res_single_2ndorder_boot, .quantity = "CI_bca")

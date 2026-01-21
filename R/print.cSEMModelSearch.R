#' `cSEMModelSearch` method for `print()`
#'
#' The `cSEMModelSearch` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [doModelSearch()]
#'
#' @export
#' @keywords internal
print.cSEMModelSearch <- function(x, ...) {
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview")
  )
  
  cat2(
    col_align("\n\n\tModel selection criterion used", width = 36), "= ", 
    x$model_selection_criterion,
    col_align("\n\tPenalize value (fbar)", width = 36), "= ", 
    x$fbar,
    col_align("\n\tFitness value of the original model", width = 36), "= ", 
    x$original_fitness
    # col_align("\n\tPermutation seed", width = 36), "= ", 
    # info$Information_permutation$Permutation_seed,
    # col_align("\n\n\tTotal bootstrap runs", width = 37), "= ", 
    # info$Information_bootstrap$Total_runs[[1]],
    # "\n\tAdmissible bootstrap results:"
  )
}
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
    rule2("Overview AGAS-PLS")
  )
  
  cat2(
    col_align("\n\n\tPopulation size", width = 36), "= ", 
    x$Information$pop_size,
    col_align("\n\n\tNumber of generations", width = 36), "= ", 
    x$Information$n_generation,
    col_align("\n\n\tMutation probability", width = 36), "= ", 
    x$Information$prob_mutation,
    col_align("\n\n\tCrossover probability", width = 36), "= ", 
    x$Information$prob_crossover,
    col_align("\n\tPenalize value (fbar)", width = 36), "= ", 
    format(x$Information$fbar, scientific = FALSE),
    col_align("\n\n\tModel selection criterion used", width = 36), "= ", 
    x$Information$model_selection_criterion,
    col_align("\n\tFitness value of the original model", width = 36), "= ", 
    x$original_fitness
  )
}
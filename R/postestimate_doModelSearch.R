#' Do a model specification search
#'
#' \lifecycle{"experimental"}
#'
#' Perform a model search 
#'
#' @usage doModelSearch(.object = NULL, 
#' .popsize = 20,
#' .iter_max= 20,
#' .mutation_prob = 0.5,
#' .cross_prob = 0.8,
#' .only_structural = FALSE,
#' .ms_criterion = 'bic')
#'
#' @return A list ???
#' 
#' @inheritParams csem_arguments
#' 
#' @seealso [cSEMResults]
#' 
#' @references \insertAllCited{}
#' 
#' @export
doModelSearch <- function(.object = NULL, 
                          .popsize = 20, 
                          .iter_max= 20,
                          .mutation_prob = 0.5,
                          .cross_prob = 0.8,
                          .only_structural = FALSE,
                          .ms_criterion = 'bic') {
  
  if(.object$Information$Model$model_type != 'Linear'){
    stop2('Currently, `doModelSearch()` supports only linear models.')
  }
  
  if(.object$Information$Arguments$.approach_weights != "PLS-PM"){
    stop2("Currently, `doModelSearch()` supports only PLS-PM.")
  }
  
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop2(
      "Package `GA` required. Use `install.packages(\"GA\")` and rerun.")
  }
  
  model_criteria <- calculateModelSelectionCriteria(
    .object = .object,
    .by_equation = FALSE,
    .only_structural = .only_structural,
    .ms_criterion = .ms_criterion
  )
  
#   HIER NOCH ABFRAGE EINBAUEN FUER DIE VERSCHIEDNEN MS CRITERIA
  BIC = model_criteria$BIC
  print(BIC)
  
  ga_control <- GA::ga(
    type = "binary",
    nBits = length(.object$Information$Model$structural),
    popSize = .popsize,
    maxiter = .iter_max,
    pmutation = .mutation_prob,
    pcrossover = .cross_prob,
    fitness = function(x) agas_fitness1(.matrix_vector = x,
                                        .data = .object$Information$Data,
                                        .model_org = .object$Information$Model,
                                        .only_structural = .only_structural,
                                        .ms_criterion = .ms_criterion),
    elitism = TRUE,
    parallel = FALSE,
    # mutation = function(object, parent) .agas_mutation(object, parent,
    #                                                    .n_variables = nr_c,
    #                                                    .n_exogenous = .n_exogenous),

    mutation = function(object, parent) .agas_mutation1(.object = object,
                                                        .parent = parent,
                                                       .model_org = .object$Information$Model),
    keepBest = TRUE
  )
  
#   HERE IS STILL A PROBLEM. HOW IS THE @bestSol STRUCTURED?
  
  best <- ga_control@bestSol[[.iter_max]]
  best_matrix <- matrix(best, nrow = nrow(.object$Information$Model$structural), 
                        ncol = ncol(.object$Information$Model$structural), 
                        byrow = TRUE,
                        dimnames = dimnames(.object$Information$Model$structural))
  
  best_fitness <- - ga_control@fitnessValue
  
  out <- list(
    best_matrix  = best_matrix,
    best_fitness = best_fitness
  )
  
  out
  
}
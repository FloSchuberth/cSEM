#' Do a model specification search
#'
#' \lifecycle{experimental}
#'
#' Perform a model search 
#'
#' @usage doModelSearch(.object = NULL, 
#' .popsize = 20,
#' .iter_max= 20,
#' .mutation_prob = 0.5,
#' .cross_prob = 0.8,
#' .only_structural = FALSE,
#' .ms_criterion = 'bic',
#' .seed = args_default()$.seed,
#' .fbar = -10000)
#'
#' @return A list containing a list of model(s) showing the 'best' value of the model 
#' selection criterion, the value of the model selection criteria, and a list of cSEM_models  
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
                          .ms_criterion = 'bic',
                          .seed = args_default()$.seed,
                          .fbar = -10000) {
  
  
  
  # if(any(class(.object) %in% "cSEMResults_default")) {
  
  model_original <- .object$Information$Model
  
  #   ADD ADDITIONAL CHECKs CSEM OBKECT
  
  if(model_original$model_type != 'Linear'){
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
    nBits = length(model_original$structural),
    popSize = .popsize,
    maxiter = .iter_max,
    pmutation = .mutation_prob,
    pcrossover = .cross_prob,
    fitness = function(x) calculatefitness(.matrix_vector = x,
                                           .data = .object$Information$Data,
                                           .model_org = model_original,
                                           .only_structural = .only_structural,
                                           .ms_criterion = .ms_criterion,
                                           .fbar = .fbar),
    elitism = TRUE,
    parallel = FALSE,
    mutation = function(object, parent) mutateVector(.object = object,
                                                      .parent = parent,
                                                      .model_org = model_original),
    keepBest = TRUE,
    seed = .seed
  )
  
  
  best <- ga_control@bestSol[[.iter_max]]
  
  best_list <- apply(best, 1, function(x) {
    matrix(
      x,
      nrow = nrow(model_original$structural),
      ncol = ncol(model_original$structural),
      byrow = TRUE,
      dimnames = dimnames(model_original$structural)
    )
  },simplify = FALSE)
  
  best_fitness <- - ga_control@fitnessValue
  
  #   Convert best_list into models that can be estimated by csem
  best_model <- lapply(best_list,function(x){
    model <- model_original
    model$structural <- x
    model
  })
  
  out <- list(
    best_matrix  = best_list,
    best_fitness = best_fitness,
    best_model = best_model
  )
  
  out
  
}
#' Automated model specification search
#'
#' \lifecycle{experimental}
#'
#' Performs an automated model specification search using a genetic algorithm \insertCite{Holland1992}{cSEM} in combination 
#' with partial least squares path modeling (PLS-PM).
#' Specifically, the function implements the approach presented in \insertCite{Trinchera2026;textual}{cSEM}. 
#' The model search is limited to the structural model.
#' Thereby it is ensured that the resulting structural model is recursive, i.e., feedback loops and non-recursive structures are prohibited.  
#' Moreover, the set of exogenous constructs is taken from the original model provided to `doModelSearch()` and remains exogenous during the search.
#' If a model violates the restrictions or does not pass [verify()] during the model search, a penalty fitness value `.fbar` is assigned to the model.
#' Therefore, ensure that `.fbar` is sufficiently small. 
#' The fitness of the model is determined by the criterion provided to `.ms_criterion`.
#' 
#' To estimate the found model `Inputcsem` results can be used, see the example.
#' 
#'  
#'
#' @usage doModelSearch(
#'  .object          = NULL, 
#'  .pop_size        = 20,
#'  .n_generations   = 20,
#'  .prob_mutation   = 0.5,
#'  .prob_crossover  = 0.8,
#'  .fbar            = -100000,
#'  .ms_criterion    = c('bic','aic','hq'),
#'  .seed            = NULL,
#'  .only_structural = FALSE
#' )
#'
#' @return A list of class `cSEMModelSearch`. The list contains the follwing elements:
#' \describe{
#'   \item{`$Information`}{A list of input parameters to `doModelSearch`}
#'   \item{`$Results`}{A list with elements, `original_fitness`, `best_fitness`, and `best_matrix`}
#'   \item{`$Inputcsem`}{A list with elements, `data`, and `model`, where `data` 
#'   comprises the dataset used to estimate the original model and `model` is a list of [cSEMModel] lists.}
#' }
#' 
#' 
#' @inheritParams csem_arguments
#' @param .ms_criterion Character string. A single character string naming the model selection criterion to compute. 
#' One of `"bic"`, `"aic"`, or `"hq"`. Defaults to `"bic"`
#' 
#' @references \insertAllCited{}
#' 
#' @seealso [cSEMResults] [calculateModelSelectionCriteria]
#' 
#' @example inst/examples/example_doModelSearch.R
#' 
#' @export
doModelSearch <- function(.object = NULL, 
                          .pop_size = 20, 
                          .n_generations= 20,
                          .prob_mutation = 0.5,
                          .prob_crossover = 0.8,
                          .fbar = -100000,
                          .ms_criterion = c('bic','aic','hq'),
                          .seed = NULL,
                          .only_structural = FALSE) {
  
  
# Check arguments
  
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop2(
      "Package `GA` required. Use `install.packages(\"GA\")` and rerun.")
  }
  
  if(any(length(.pop_size)!=1 | !is.numeric(.pop_size)| .pop_size<0)){
    stop2("`.pop_size` must be a scalar.")
  }
  
  if(any(length(.n_generations)!=1 | !is.numeric(.n_generations) | .n_generations<0)){
    stop2("`.n_generations` must be a scalar.")
  }
  
  if(any(!is.numeric(.prob_mutation)|.prob_mutation>1|.prob_mutation<0|length(.prob_mutation)>1)){
    stop2("`.prob_mutation` must be a scalar rangning between 0 and 1.")
  }

  if(any(!is.numeric(.prob_crossover)|.prob_crossover>1|.prob_crossover<0|length(.prob_crossover)>1)){
    stop2("`.prob_crossover` must be a scalar rangning between 0 and 1.")
  }
  
  if(any(length(.fbar)>1|!is.numeric(.fbar))){
    stop2("`.fbar` must be a scalar.")
  }
  
  .ms_criterion <- match.arg(.ms_criterion,several.ok = FALSE) 

  if(inherits(.object, "cSEMResults_multi")) {
    
    stop2("Multigroup objects are currently not supported.") 
    
  }else if(inherits(.object, "cSEMResults_2ndorder")) {
    
    stop2("Models containing second-order constructs are currently not supported.") 
    
  }else if(inherits(.object, "cSEMResults_default")) {
    
    model_original <- .object$Information$Model    
  }else{
    
    stop2(
      "The following error occured in the assess() function:\n",
      "`.object` must be a `cSEMResults` object."
    )
  }
  

  if(model_original$model_type != 'Linear'){
    stop2('Currently, `doModelSearch()` supports only linear models.')
  }
  
  if(.object$Information$Arguments$.approach_weights != "PLS-PM"){
    stop2("Currently, `doModelSearch()` supports only PLS-PM.")
  }
  
  
  fitness <- calculateModelSelectionCriteria(
    .object = .object,
    .by_equation = FALSE,
    .only_structural = .only_structural,
    .ms_criterion = .ms_criterion
  )[[1]]
  
  ga_control <- GA::ga(
    type = "binary",
    nBits = length(model_original$structural),
    popSize = .pop_size,
    maxiter = .n_generations,
    pmutation = .prob_mutation,
    pcrossover = .prob_crossover,
    fitness = function(x) calculateFitness(.matrix_vector = x,
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
  
  
  best <- ga_control@bestSol[[.n_generations]]
  
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
    "Information" = list(
      model_selection_criterion = .ms_criterion,
      fbar = .fbar,
      pop_size = .pop_size,
      n_generations = .n_generations,
      prob_mutation = .prob_mutation,
      prob_crossover= .prob_crossover, 
      seed = .seed,
      only_structural=.only_structural
    ),
    "Results" = list(
    original_fitness = fitness,  
    best_fitness = best_fitness,
    best_matrix  = best_list
    ),
    "Inputcsem"=list(
    data = .object$Information$Data,
    model = best_model
    )
  )
  
  class(out) <- 'cSEMModelSearch'
  
  out
  
}
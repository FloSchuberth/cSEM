# --- Internal helpers (same file, not exported) ------------------------------

# --- Internal state (no globals leaked) --------------------------------------
.pkg_state <- new.env(parent = emptyenv())
.pkg_state$best_individuals_all <- list()
.pkg_state$best_individual      <- NULL
.pkg_state$best_fitness         <- -Inf

# Function to transform the measurement model into a list
.transform_measurement_model <- function(.mes_mod) {
  measurement_model_list <- list()
  for (i in 1:nrow(.mes_mod)) {
    latent_var <- rownames(.mes_mod)[i]
    indicators <- colnames(.mes_mod)[which(.mes_mod[i, ] == 1)]
    measurement_model_list[[latent_var]] <- indicators
  }
  return(measurement_model_list)
}

.agas_mutation <- function(.object, .parent, .n_variables, .mutation_prob, .n_exogenous) {
  mutate <- .parent <- as.vector(.object@population[.parent,])
  mutate_matrix <- matrix(mutate, nrow = .n_variables, byrow = TRUE)
  
  diag(mutate_matrix) <- 0  
  
  # Ensure the first three rows are all zeros
  for (i in 1:.n_exogenous) {
    mutate_matrix[i, ] <- 0
  }
  
  mutate_vector <- as.vector(t(mutate_matrix))
  
  # Create list indices that i want to mutate
  indices <- which(!diag(.n_variables), arr.ind = TRUE)
  indices <- indices[indices[, 1] > .n_exogenous, ]  # Exclude first three rows
  
  # Convert row and column indices to vector indices
  if (length(indices) > 0) {
    vector_indices <- (indices[, 1] - 1) * .n_variables + indices[, 2]
    
    # Select a random index from the sub-diagonal indices 
    if (length(vector_indices) > 0 && runif(1) <= .mutation_prob) {  # Use mutation_prob here
      available_indices <- vector_indices  # Keep track of available indices to flip
      while (length(available_indices) > 1) {
        j <- sample(available_indices, size = 1)
        mutate_vector[j] <- abs(mutate_vector[j] - 1)
        
        # Check if the mutation creates a cycle
        adj_matrix <- matrix(mutate_vector, nrow = .n_variables, byrow = TRUE)
        has_cycle <- .has_cycle_matrix(adj_matrix)
        
        if (has_cycle) {
          mutate_vector[j] <- abs(mutate_vector[j] - 1)
          available_indices <- setdiff(available_indices, j)
        } else {
          break  
        }
      }
    }
  }
  
  return(mutate_vector)
}


# --- DFS cycle routine -------------------------------------------------------
#' @keywords internal
.has_cycle_matrix <- function(.adj) {
  n <- nrow(.adj)
  # 0 = unvisited, 1 = visiting, 2 = done
  state <- integer(n)
  
  dfs <- function(u) {
    if (state[u] == 1) return(TRUE)   # back-edge → cycle
    if (state[u] == 2) return(FALSE)  # already processed, no cycle here
    
    state[u] <<- 1  # mark as visiting
    
    neighbors <- which(.adj[u, ] != 0)
    for (v in neighbors) {
      if (dfs(v)) return(TRUE)
    }
    
    state[u] <<- 2  # done
    FALSE
  }
  
  for (i in seq_len(n)) {
    if (state[i] == 0 && dfs(i)) {
      return(TRUE)
    }
  }
  
  FALSE
}

# --- Matrix utils ------------------------------

#' @keywords internal
.matrix_to_string <- function(.m) {
  row_strs <- apply(.m, 1, paste, collapse = ",")
  paste(row_strs, collapse = ";")
}

#' @keywords internal
.check_matrix_criteria <- function(.adj_matrix, .n_exogenous) {
  if (any(rowSums(.adj_matrix[1:.n_exogenous, , drop = FALSE]) != 0)) return(FALSE) # exogenous all zeros
  if (any(diag(.adj_matrix) != 0)) return(FALSE)                       # diag all zeros
  TRUE
}


#' @keywords internal
.check_matrix <- function(.mat) {
  n <- nrow(.mat)
  for (i in seq_len(n)) if (all(.mat[i, ] == 0) && all(.mat[, i] == 0)) return(FALSE)
  TRUE
}

# --- Create model string from matrix -----------------------------------------

#' @keywords internal
.create_sem_model_string_from_matrix <- function(.adj_matrix, .measurement_model, .variables) {
  .structural_coefficients <- list()
  .type_of_variable <- setNames(rep("composite", length(.variables)), .variables)
  
  model_string <- "# Composite model\n" # Initialize the model string
  
  # Include the measurement model specified in the input for composite types
  for (var in names(.measurement_model)) {
    if (.type_of_variable[var] == "composite") {
      items <- paste(.measurement_model[[var]], collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s <~ %s\n", var, items), sep = "")
    }
  }
  
  # Adding reflective measurement models
  model_string <- paste(model_string, "\n# Reflective measurement model\n")
  for (var in names(.measurement_model)) {
    if (.type_of_variable[var] == "reflective") {
      items <- paste(.measurement_model[[var]], collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s =~ %s\n", var, items), sep = "")
    }
  }
  
  # Structural model section
  model_string <- paste(model_string, "\n# Structural model\n")
  for (i in seq_along(.variables)) {
    dependent <- .variables[i]
    predictors <- .variables[.adj_matrix[i, ] == 1]
    if (length(predictors) > 0) {
      relationship_str <- paste(predictors, collapse = " + ")
      model_string <- paste(model_string, sprintf("  %s ~ %s\n", dependent, relationship_str), sep = "")
    }
  }
  
  model_string <- gsub("\\n\\s+$", "", model_string) # Cleanup the string
  model_string <- trimws(model_string)  
  
  return(model_string)
}

# --- Fitness ---------------------------------------------------

#' @keywords internal
.agas_fitness <- function(.matrix_vector, .dataset_generated, .n_exogenous, .measurement_model, .variables, .only_structural) {
  n_variables <- length(.variables)
  adj_matrix  <- matrix(.matrix_vector, nrow = n_variables, byrow = TRUE)
  
  if (!.check_matrix(adj_matrix)) return(-100000)
    # adj_matrix <- .repair_individual_unused(adj_matrix)
  
  if (!.check_matrix_criteria(adj_matrix, .n_exogenous)) return(-100000)
  
  if (.has_cycle_matrix(adj_matrix)) return(-100000)
  
  model_string <- .create_sem_model_string_from_matrix(adj_matrix, .measurement_model, .variables)
  out <- csem(.data = .dataset_generated, .model = model_string)
  ver <- verify(out)
  if (!sum(ver) == 0) return(-100000)
  
  model_criteria <- calculateModelSelectionCriteria(
    out, .by_equation = FALSE, .only_structural = .only_structural
  )
  
  sem_fitness <- -model_criteria$BIC
  
  
  if (is.na(sem_fitness)) return(-100000)
  
  if (sem_fitness > .pkg_state$best_fitness) {
    .pkg_state$best_individual      <- adj_matrix
    .pkg_state$best_fitness         <- sem_fitness
    .pkg_state$best_individuals_all <- append(.pkg_state$best_individuals_all, list(adj_matrix))
  }
  
  sem_fitness
}

# --- Public function ---------------------------------------------------------

#' Do a model search
#'
#' \lifecycle{stable}
#'
#' Perform a model search \insertCite{Hair2016;textual}{cSEM}
#'
#' @usage doModelSearch(.object = NULL)
#'
#' @return The mean matrix generated by the AGAS-PLS runs
#' @inheritParams csem_arguments
#' @seealso [cSEMResults]
#' @references \insertAllCited{}
#' @export
doModelSearch <- function(.object = NULL, 
                          .n_exogenous = 3,
                          .popsize = 20, 
                          .maxiter=20,
                          .mutation_prob = 0.5,
                          .seeds = c(2, 4, 5), 
                          .only_structural = TRUE) {
  if (is.null(.object)) stop("`.object` must be a cSEM results object.")
  
  # Retrieve information on the model (number variables, name variables)
  path_estimates <- .object$Estimates$Path_estimates 
  .variables <- rownames(path_estimates)
  .n_variables <- length(.variables)
  .measurement_model <- .transform_measurement_model(.object$Information$Model$measurement)
  
  # compute criteria once
  model_criteria <- calculateModelSelectionCriteria(
    .object,
    .by_equation = FALSE,
    .only_structural = .only_structural
  )
  
  BIC = model_criteria$BIC
  print(BIC)
  
  .agas_dataset <- .object$Information$Data
  
  seeds <- as.integer(.seeds)
  mat_list <- vector("list", length(.seeds))
  names(mat_list) <- as.character(.seeds)
  
  for (i in seq_along(.seeds)) {
    
    .pkg_state$best_individuals_all <- list()
    .pkg_state$best_individual      <- NULL
    .pkg_state$best_fitness         <- -Inf
    
    ga_control <- GA::ga(
        type = "binary",
        nBits = .n_variables * .n_variables,
        popSize = .popsize,
        maxiter = .maxiter,
        pmutation = 1.0,
        pcrossover = 0.8,
        fitness = function(x) .agas_fitness(x, .agas_dataset, .n_exogenous, .measurement_model, .variables, .only_structural),
        elitism = TRUE,
        parallel = FALSE,
        seed = i,
        mutation = function(object, parent) .agas_mutation(object, parent, .n_variables, .mutation_prob, .n_exogenous)
      )
    mat_list[[i]] <- .pkg_state$best_individual
  }
  
  arr <- array(unlist(mat_list),
               dim = c(.n_variables, .n_variables, length(mat_list)))
  mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  rownames(mean_mat) <- .variables
  colnames(mean_mat) <- .variables
  
  print(mean_mat)
  out <- .matrix_to_string(mean_mat)
  return(out)
  
}

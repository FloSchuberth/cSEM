#' Internal: AGAS MUTATION
.agas_mutation <- function(.object, .parent, .n_variables, .n_exogenous) {
  
#   Would be nice to label the matrix. Perhaps additional argument is needed. 
  mutate <- .parent <- as.vector(.object@population[.parent,])
  mutate_matrix <- matrix(mutate, nrow = .n_variables, byrow = TRUE)
  
  diag(mutate_matrix) <- 0  
  
  # Ensure the first three rows are all zeros
  for (i in 1:.n_exogenous) {
    mutate_matrix[i, ] <- 0
  }
  
  # Comment in later
  # mutate_matrix[,.n_exogenous,drop=FALSE] <- 0
  
#   Try to work here with matices instead of vectors
  mutate_vector <- as.vector(t(mutate_matrix))
  
  # Create list indices that i want to mutate
  indices <- which(!diag(.n_variables), arr.ind = TRUE)
  indices <- indices[indices[, 1] > .n_exogenous, ]  # Exclude first three rows
  
  # Convert row and column indices to vector indices
  if (length(indices) > 0) {
    vector_indices <- (indices[, 1] - 1) * .n_variables + indices[, 2]
    
    # Select a random index from the sub-diagonal indices 
    if (length(vector_indices) > 0) {  # Use mutation_prob here
      available_indices <- vector_indices  # Keep track of available indices to flip
      while (length(available_indices) > 1) {
        j <- sample(available_indices, size = 1)
        mutate_vector[j] <- abs(mutate_vector[j] - 1)
        
        # Check if the mutation creates a cycle
        struc_model <- matrix(mutate_vector, nrow = .n_variables, byrow = TRUE)
        has_cycle <- has_cycle_matrix(.matrix = struc_model)
        
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


#' Internal: AGAS MUTATION
.agas_mutation1 <- function(.object, 
                            .parent, 
                            .model_org) {
  cons_exo <- .model_org$cons_exo

  mutate <- .parent <- as.vector(.object@population[.parent,])
  
  
  mutate_matrix <- matrix(mutate, 
                          nrow = nrow(.model_org$structural),
                          ncol = ncol(.model_org$structural),
                          byrow = TRUE,
                          dimnames = dimnames(.model_org$structural))
  
  diag(mutate_matrix) <- 0  
  
  # Ensure that the rows of the exogenous constructs are zero
  mutate_matrix[cons_exo,] <- 0

#   Create matrix indicating the elements in the mutate matrix that are allowed to mutate
ind_matrix <- mutate_matrix  
ind_matrix[] <- 1
diag(ind_matrix) <- 0
ind_matrix[cons_exo,] <- 0  
indices <- which(ind_matrix == 1, arr.ind = TRUE)

while(nrow(indices)>0){
  randn <- sample(seq(nrow(indices)),1)

# Flip element
  mutate_matrix[indices[randn,,drop=FALSE]] <- abs(mutate_matrix[indices[randn,,drop=FALSE]] - 1)

has_cycle <- has_cycle_matrix(.matrix = mutate_matrix)

if(has_cycle){
  mutate_matrix[indices[randn,,drop=FALSE]] <- abs(mutate_matrix[indices[randn,,drop=FALSE]] - 1)
  indices <- indices[-randn,,drop=FALSE]
}else{
  break
}
}

return(as.vector(t(mutate_matrix)))
}

#' Internal: DFS cycle routine 
#' @keywords internal
has_cycle_matrix <- function(.matrix) {
  n <- nrow(.matrix)
  # 0 = unvisited, 1 = visiting, 2 = done
  state <- integer(n)
  
  dfs <- function(u) {
    if (state[u] == 1) return(TRUE)   # back-edge → cycle
    if (state[u] == 2) return(FALSE)  # already processed, no cycle here
    
    state[u] <<- 1  # mark as visiting
    
    neighbors <- which(.matrix[u, ] != 0)
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

#' Internal:  --- Matrix utils ------------------------------

#' @keywords internal 
.matrix_to_string <- function(.matrix) { #IS THIS FUNCTION NEEDED AT ALL?
  row_strs <- apply(.matrix, 1, paste, collapse = ",")
  paste(row_strs, collapse = ";")
}

#' Internal:
#' @keywords internal
check_matrix_criteria <- function(.adj_matrix, .n_exogenous) {
  if (any(rowSums(.adj_matrix[1:.n_exogenous, , drop = FALSE]) != 0)) return(FALSE) # exogenous all zeros
  if (any(diag(.adj_matrix) != 0)) return(FALSE)                       # diag all zeros
  TRUE
}

check_matrix_criteria1 <- function(.matrix, .model_org) {
  
  cons_exo <- .model_org$cons_exo
  
  if (any(rowSums(.matrix[cons_exo, , drop = FALSE]) != 0)) {
    return(FALSE)
    }
  if (any(diag(.matrix) != 0)) {
    return(FALSE)
    }
  TRUE
}



#' @keywords internal
check_matrix <- function(.matrix) {
  n <- nrow(.matrix)
  for (i in seq_len(n)) if (all(.matrix[i, ] == 0) && all(.matrix[, i] == 0)) return(FALSE)
  TRUE
}

# --- Create model string from matrix -----------------------------------------

#' @keywords internal
create_sem_model_string_from_matrix <- function(.adj_matrix, .measurement_model, .variables) {
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

#' Internal:  --- Fitness ---------------------------------------------------

#' @keywords internal
agas_fitness <- function(.matrix_vector, .data, .n_exogenous, .measurement_model, .variables, .only_structural) {
  n_variables <- length(.variables)
  adj_matrix  <- matrix(.matrix_vector, nrow = n_variables, byrow = TRUE)
  
  if (!check_matrix(.matrix = adj_matrix)) return(-100000)
  # adj_matrix <- .repair_individual_unused(adj_matrix)
  
  if (!check_matrix_criteria(adj_matrix, .n_exogenous)) return(-100000)
  
  if (has_cycle_matrix(.matrix = adj_matrix)) return(-100000)
  
  model_string <- create_sem_model_string_from_matrix(adj_matrix, .measurement_model, .variables)
  out <- csem(.data = .data, .model = model_string)
  ver <- verify(out)
  if (!sum(ver) == 0) return(-100000)
  
  model_criteria <- calculateModelSelectionCriteria(
    out, .by_equation = FALSE, .only_structural = .only_structural
  )
  
  sem_fitness <- -model_criteria$BIC
  
  if (is.na(sem_fitness)) return(-100000)
  
  sem_fitness
}

agas_fitness1 <- function(.matrix_vector, .data, .model_org , .only_structural, .ms_criterion) {
  
  # cons_exo <- .model_org$cons_exo
  # cons_endo <- .model_org$cons_endo
  
  # n_cons <- length(unique(c(cons_exo,cons_endo)))
  
  # n_cons <- length(.variables)
  struc_model  <- matrix(.matrix_vector, 
                         nrow = nrow(.model_org$structural), 
                         ncol = ncol(.model_org$structural),
                         byrow = TRUE,
                         dimnames=dimnames(.model_org$structural))
  
  if (!check_matrix(.matrix = struc_model)) {
    return(-100000)
  }
  # struc_model <- .repair_individual_unused(struc_model)
  
  if (!check_matrix_criteria1(.matrix = struc_model, .model_org = .model_org)){
    return(-100000)
  } 
  
  if (has_cycle_matrix(.matrix = struc_model)){
    return(-100000)
    } 
  
  model <- .model_org
  model$structural <- struc_model
  
  # model_string <- create_sem_model_string_from_matrix(struc_model, .measurement_model, .variables)
  out <- csem(.data = .data, .model = model)
  ver <- verify(out)
  if (!sum(ver) == 0){
    return(-100000)
  }
  
  model_criteria <- calculateModelSelectionCriteria(.object = out,
                                                    .ms_criterion = .ms_criterion,
                                                    .by_equation = FALSE,
                                                    .only_structural = .only_structural)
  
  
#   ADJUST HERE TO ALLOW FOR DIFFEREN MS CRITERIA
  sem_fitness <- -model_criteria$BIC
  
  
#   WHAT IS THIS CASE?
  if (is.na(sem_fitness)){
    sem_fitness <- -100000
  }
  
  
  
  sem_fitness
}

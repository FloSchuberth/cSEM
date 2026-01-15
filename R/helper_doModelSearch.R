#' Internal: AGAS MUTATION
.agas_mutation <- function(.object, 
                            .parent, 
                            .model_org) {
  cons_exo <- .model_org$cons_exo

# REMOVE AS.VECTOR?
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
check_matrix_criteria <- function(.matrix, .model_org) {
  
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

#' @keywords internal
agas_fitness <- function(.matrix_vector, .data, .model_org , .only_structural, .ms_criterion) {
  
  struc_model  <- matrix(.matrix_vector, 
                         nrow = nrow(.model_org$structural), 
                         ncol = ncol(.model_org$structural),
                         byrow = TRUE,
                         dimnames=dimnames(.model_org$structural))
  
  if (!check_matrix(.matrix = struc_model)) {
    return(-100000)
  }

  if (!check_matrix_criteria(.matrix = struc_model, .model_org = .model_org)){
    return(-100000)
  } 
  
  if (has_cycle_matrix(.matrix = struc_model)){
    return(-100000)
    } 
  
  model <- .model_org
  model$structural <- struc_model
  
#   MUST STRUCTURAL MODEL ALWAYS BE LOWER TRIANGUALR MATRIX (FLO CHECK THIS)
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
  
  sem_fitness
}

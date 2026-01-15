#' Internal: Function that mutates a vector
#' 
#' Function that randomly mutates one of the elements of the vectorized structural model by flipping it
#'  from 0 to 1 or vice versa) taking into account that not every element is allowed to be mutated, The elements in the rows of exogenous constructs 
#' and the elements on the diagonal are not mutated. Similarly, elements are not mutated when this would result in 
#' feedback loops.
#' This function adheres to the necessary requirements  for a function to be supplied to \code{\link[GA]{ga}}.
#' 
#' @usage mutateVector(.object,
#'                    .parent,
#'                    .model_org)
#'                      
#' @return A numeric vector; Vectorized version of the structural model in which one element has been flipped from 0 to 1 or vice versa
#' 
#' @keywords internal
#' 
mutateVector <- function(.object, 
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
    
    has_cycle <- checkCycles(.matrix = mutate_matrix)
    
    if(has_cycle){
      mutate_matrix[indices[randn,,drop=FALSE]] <- abs(mutate_matrix[indices[randn,,drop=FALSE]] - 1)
      indices <- indices[-randn,,drop=FALSE]
    }else{
      break
    }
  }
  
  return(as.vector(t(mutate_matrix)))
}


#' Internal: Function that calculates the fitness value 
#' 
#' Function that calculates the fitness value of a model. 
#' This function adheres to the necessary requirements  for a function to be supplied to \code{\link[GA]{ga}}.
#' 
#' @return numeric
#' 
#' @keywords internal
calculatefitness <- function(.matrix_vector, 
                             .data, 
                             .model_org,
                             .only_structural,
                             .ms_criterion,
                             .fbar){
  
  struc_model  <- matrix(.matrix_vector, 
                         nrow = nrow(.model_org$structural), 
                         ncol = ncol(.model_org$structural),
                         byrow = TRUE,
                         dimnames=dimnames(.model_org$structural))
  
  
  # if (!.check_matrix(adj_matrix)) return(-100000)
  # 
  # if (!.check_matrix_criteria(adj_matrix, .n_exogenous)) return(-100000)
  # 
  # if (.has_cycle_matrix(adj_matrix)) return(-100000)
  
  if (!checkIsolatedConstruct(.matrix = struc_model)|
      any(rowSums(struc_model[.model_org$cons_exo, , drop = FALSE]) != 0)|
      any(diag(struc_model)!=0)|
      checkCycles(.matrix = struc_model)) {
    
    return(.fbar)
  }
  
  model <- .model_org
  model$structural <- struc_model
  
  #   MUST STRUCTURAL MODEL ALWAYS BE LOWER TRIANGUALR MATRIX (FLO CHECK THIS)
  out <- csem(.data = .data, .model = model)
  if (sum(verify(out)) != 0){
    return(.fbar)
  }
  
  model_criteria <- calculateModelSelectionCriteria(.object = out,
                                                    .ms_criterion = .ms_criterion,
                                                    .by_equation = FALSE,
                                                    .only_structural = .only_structural)
  
  
  #   ADJUST HERE TO ALLOW FOR DIFFEREN MS CRITERIA
  sem_fitness <- -model_criteria$BIC
  
  sem_fitness
}



#' Internal: Function that search for cycles
#' 
#' Function that searches for cycles in a binary indicator matrix using the depth first search (DFS).
#' 
#' @return Logical vector `'TRUE'`: Cycle detected; `'FALSE'` No cycle detected
#'  
#' @keywords internal
checkCycles <- function(.matrix) {
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

#  check_matrix_criteria <- function(.matrix, .model_org) {
#   
#   cons_exo <- .model_org$cons_exo
#   
#   if (any(rowSums(.matrix[cons_exo, , drop = FALSE]) != 0)) {
#     return(FALSE)
#     }
#   if (any(diag(.matrix) != 0)) {
#     return(FALSE)
#     }
#   TRUE
# }

#'Internal: Check whether construct in the sturtcural model is isolated
#'@return logical: `'TRUE'`: Isolated; `'FALSE'`:connected to at least one variable via a path.

#' @keywords internal
checkIsolatedConstruct <- function(.matrix) {
  n <- nrow(.matrix)
  
  for (i in seq_len(n)){
    if(all(.matrix[i, ] == 0) && all(.matrix[, i] == 0)){
      return(FALSE)
    } 
  } 
  TRUE
}



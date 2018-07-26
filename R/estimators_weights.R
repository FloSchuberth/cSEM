#' Calculate composite weights
#'
#' Calculates weights
#'
#' More description here
#'
#' @usage calculateWeights(
#'   .data                     = NULL,
#'   .model                    = NULL,
#'   .PLS_mode                 = NULL,
#'   .tolerance                = NULL,
#'   .iter_max                 = NULL,
#'   .PLS_weight_scheme_inner  = NULL,
#'   .ignore_structural_model  = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{A (J x J) matrix of inner weights.}
#'   \item{`$Modes`}{A named vector of Modes used for the outer estimation.}
#' }
#' @export
#'

calculateWeightsPLS <- function(
  .data                     = NULL,
  .model                    = NULL,
  .PLS_mode                 = NULL,
  .tolerance                = NULL,
  .iter_max                 = NULL,
  .PLS_weight_scheme_inner  = NULL,
  .ignore_structural_model  = NULL
) {


  ### Make the function more autonomous ========================================
  ## Convert model to "cSEMModel" format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    csem_model <- parseModel(.model)
  } else {
    csem_model <- .model
  }

  ## Prepare, standardize, check, and clean data if not already in this format
  if(!(class(.data) == "cSEMData")) {
    if(is.matrix(.data) && isSymmetric.matrix(.data)) {
      S <- .data
    } else {
      #  If function is used externally, data is not automatically scaled!
      X <- processData(.data = .data, .model = csem_model)
      S <- stats::cov(X)
    }
  }

  ### Preparation ==============================================================
  ## Get/set the modes for the outer estimation

  if(is.null(.PLS_mode)) {
    
    modes <- ifelse(csem_model$construct_type == "Common factor", "ModeA", "ModeB")

  } else if(all(.PLS_mode %in% c("ModeA", "ModeB"))) {
    
    if(setequal(names(.PLS_mode), names(csem_model$construct_type))) {
      modes <- .PLS_mode
      modes <- modes[names(csem_model$construct_type)]
    } else if(length(.PLS_mode) == 1) {
      modes <- rep(.PLS_mode, length(csem_model$construct_type))
      names(modes) <- names(csem_model$construct_type)
    } else if(length(setdiff(names(.PLS_mode), names(csem_model$construct_type))) > 0) {
      stop(paste0("`", setdiff(names(.PLS_mode), names(csem_model$construct_type)), 
                  "`", collapse = ", ")," in `.PLS_mode` is an unknown construct name.", call. = FALSE)
    } else  {
      stop("Mode ", paste0("`", setdiff(names(csem_model$construct_type), names(.PLS_mode)), 
                  "`", collapse = ", ")," in `.PLS_mode` is missing.", call. = FALSE)
    }
  } else {
    stop(paste0("`", setdiff(.PLS_mode, c("ModeA", "ModeB")), "`", collapse = ", "),
         " in `.PLS_mode` is an unknown mode.", call. = FALSE)
  }

  ### Calculation/Iteration ====================================================
  W <- csem_model$measurement
  # Scale weights
  W <- scaleWeights(.S = S, .W = W)

  W_iter       <- W
  tolerance    <- .tolerance
  iter_max     <- .iter_max
  iter_counter <- 0

  repeat {
    # Counter
    iter_counter <- iter_counter + 1

    # Inner estimation
    E <- calculateInnerWeightsPLS(
      .S                       = S,
      .W                       = W_iter,
      .csem_model              = csem_model,
      .PLS_weight_scheme_inner = .PLS_weight_scheme_inner,
      .ignore_structural_model = .ignore_structural_model
    )
    # Outer estimation

    W <- calculateOuterWeightsPLS(
      .S        = S,
      .W        = W_iter,
      .E        = E,
      .modes    = modes
    )

    # Scale weights
    W <- scaleWeights(S, W)

    # Check for convergence
    if(max(abs(W_iter - W)) < .tolerance) {
      # Set convergence status to TRUE as algorithm has converged
      conv_status = TRUE
      break # return iterative PLS weights
    } else if(iter_counter == iter_max & iter_max == 1) {
      # Set convergence status to NULL, NULL is used if no algorithm is used
      conv_status = NULL
      break # return one-step PLS weights
      
    } else if(iter_counter == iter_max & iter_max > 1) {
      # Set convergence status to FALSE, as algorithm has not converged
      conv_status = FALSE
      warning("Iteration did not converge after ", iter_max, " steps. ",
              "Last weights are returned.")
      
    } else {
      W_iter <- W
    }
  }

  # Return
  l <- list("W" = W, "E" = E, "Modes" = modes, "Conv_status" = conv_status,
            "Iterations" = iter_counter)
  return(l)

  ### For maintenance: ### -----------------------------------------------------
  # W (J x K) := (Block-)diagonal matrix of weights with the same dimension as .csem_model$measurement
  # C (J x J) := Empirical composite correlation matrix
  # E (J x J) := Matrix of inner weights
  # D (J x J) := An "adjacency" matrix with elements equal to one, if proxy j and j' are adjacent
} # END calculateWeightsPLS

#' Calculate the inner weights for PLS
#'
#' Calculates inner weights
#'
#' More description here.
#'
#' @usage calculateInnerWeightsPLS(
#'   .S                        = NULL,
#'   .W                        = NULL,
#'   .csem_model               = NULL,
#'   .PLS_weight_scheme_inner  = NULL,
#'   .ignore_structural_model  = NULL
#' )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x J) matrix of inner weights.
#'
calculateInnerWeightsPLS <- function(
  .S                        = NULL,
  .W                        = NULL,
  .csem_model               = NULL,
  .PLS_weight_scheme_inner  = NULL,
  .ignore_structural_model  = NULL
) {

  # Composite correlation matrix (C = V(H))
  C <- .W %*% .S %*% t(.W)

  # Due to floting point errors may not be symmetric anymore. In order
  # prevent that replace the lower triangular elements by the upper
  # triangular elements

  C[lower.tri(C)] <- t(C)[lower.tri(C)]

  # Path model relationship
  tmp <- rownames(.csem_model$measurement)
  E   <- .csem_model$structural[tmp, tmp]
  D   <- E + t(E)

  ## (Inner) weightning scheme:
  if(.PLS_weight_scheme_inner == "path" & .ignore_structural_model) {
    .ignore_structural_model <- FALSE
    warning("Structural model required for the path weighting scheme.\n",
            ".ignore_structural_model = TRUE was changed to FALSE.")
  }

  if(.ignore_structural_model) {
    switch (.PLS_weight_scheme_inner,
            "centroid"  = {E <- sign(C) - diag(1, nrow = nrow(.W))},
            "factorial" = {E <- C - diag(1, nrow = nrow(.W))}
    )
  } else {
    switch (.PLS_weight_scheme_inner,
            "centroid"  = {E <- sign(C) * D},
            "factorial" = {E <- C * D},
            "path"      = {
              ## All construct that have an arrow pointing to construct j
              ## are called predecessors of j; all arrows going from j to other
              ## constructs are called followers of j

              ## Path weighting scheme:
              ## Take construct j:
              #  - For all predecessors of j: set the inner weight of
              #    predescessor i to the correlation of i with j.
              #  - For all followers of j: set the inner weight of follower i
              #    to the coefficient of a multiple regression of j on
              #    all followers i with i = 1,...,I.

              E_temp <- E
              ## Assign predescessor relation
              E[t(E_temp) == 1] <- C[t(E_temp) == 1]

              ## Assign follower relation
              for(j in tmp) {

                followers <- E_temp[j, ] == 1

                if(sum(followers) > 0) {

                  E[j, followers] <-  solve(C[followers, followers, drop = FALSE]) %*%
                    C[j, followers]
                }
              }
            }
    )
  }

  return(E)
} # END calculateInnerWeights

#' Calculate the outer weights for PLS
#'
#' Calculates outer weights
#'
#' More description here..
#'
#' @usage calculateOuterWeightsPLS(
#'    .S      = NULL,
#'    .W      = NULL,
#'    .E      = NULL,
#'    .modes  = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A (J x K) matrix of outer weights.
#'

calculateOuterWeightsPLS <- function(
  .S      = NULL,
  .W      = NULL,
  .E      = NULL,
  .modes  = NULL
) {
  # Covariance/Correlation matrix between each proxy and all indicators (not
  # only the related ones). Note: Cov(H, X) = WS, since H = XW'.
  W <- .W
  proxy_indicator_cor <- .E %*% W %*% .S

  for(i in 1:nrow(W)) {
    block      <- rownames(W[i, , drop = FALSE])
    indicators <- W[block, ] != 0

    if(.modes[block] == "ModeA") {
      ## Mode A - Regression of each indicator on its corresponding proxy
      # Select only
      W[block, indicators] <- proxy_indicator_cor[block, indicators]

    } else if(.modes[block] == "ModeB") {
      ## Mode B - Regression of each proxy on all its indicator
      # W_j = S_jj^{-1} %*% Cov(eta_j, X_j)
      W[block, indicators] <- solve(.S[indicators, indicators]) %*% proxy_indicator_cor[block, indicators]

    } # END ModeB
  }
  return(W)
} # END calculateOuterWeights

#' Calculate weights using one of Kettenrings's approaches
#'
#' Calculates weights according to on of the the five criteria
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*", and "*GENVAR*"
#' suggested by \insertCite{Kettenring1971;textual}{cSEM}.
#'
#' Some more description...
#'
#' @usage calculateWeightsKettenring(.data, .model, .approach)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsKettenring <- function(
  .data           = NULL,
  .model          = NULL,
  .approach       = NULL
) {
  
  ### Make the function more autonomous ========================================
  ## Convert model to "cSEMModel" format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    csem_model <- parseModel(.model)
  } else {
    csem_model <- .model
  }
  
  ## Prepare, standardize, check, and clean data if not already in this format
  if(!(class(.data) == "cSEMData")) {
    if(is.matrix(.data) && isSymmetric.matrix(.data)) {
      S <- .data
    } else {
      #  If function is used externally, data is not automatically scaled!
      X <- processData(.data = .data, .model = csem_model)
      S <- stats::cov(X)
    }
  }
  
  # Check if all constructs are modeled as composite
  if("Common factor" %in% .model$construct_type ){
    stop("Kettenring's is only allowed for pure composite models.")
  }
  
  # read out measurement model 
  W=csem_model$measurement
  
  Construct_names=rownames(W)
  
  # read out names of indicators belonging to one block 
  Blocksqrtmcorrelation=list()
  namesIndicators=list()
  
  for(i in Construct_names){      
    # Collect indicator names that belong to construct i 
    namesIndicators[[i]] = colnames(W)[W[i,]==1]
    
    # selects Correlation matrix of block i and calculate S_{ii}^{-1/2}}
    Blocksqrtmcorrelation[[i]]=solve(expm::sqrtm(as.matrix(S[namesIndicators[[i]],namesIndicators[[i]]])))
    dimnames(Blocksqrtmcorrelation[[i]])=list(namesIndicators[[i]],namesIndicators[[i]])
  }
  
  # restructure the correlation matrix S to get the right order
  Srestructured=S[unlist(namesIndicators),unlist(namesIndicators)]
  
  # Create Blockdiagonal matrix
  H=as.matrix(Matrix::bdiag(Blocksqrtmcorrelation))
  dimnames(H)=list(unlist(namesIndicators),unlist(namesIndicators))
  
  if(.approach == 'MAXVAR'){
    # The MAXVAR implementation is based on a MATLAB code provided by Theo K. Dijkstra
    
    Rs=H%*%Srestructured%*%H
    
    # Calculate eigenvalues and eigenvectors of Rs
    u=as.matrix(eigen(Rs)$vectors[,1],nocl=1)# take eigenvector to the largest eigenvalue
    rownames(u)=  unlist(namesIndicators)
    
    # Calculate weights
    a=matrix(0,ncol=length(unlist(namesIndicators)),nrow=length(Construct_names),
             dimnames=list(Construct_names,unlist(namesIndicators)))
    for(i in Construct_names){
      a[i,namesIndicators[[i]]]=
        as.vector((H[namesIndicators[[i]],namesIndicators[[i]]]%*%
                     u[namesIndicators[[i]],])/sqrt(sum(c(u[namesIndicators[[i]],])^2)))
    }
    for(i in Construct_names){
      if(sum(a[i,])<0) a[i,]=-a[i,]
    }
    
    
    
    W= a[rownames(W),colnames(W)]
    W= scaleWeights(S, W)
    
    l <- list("W" = W, "E" = NULL, "Modes" = rep('MAXVAR',length(Construct_names)),
              "Conv_status" = TRUE, "Iterations" = 0)
    return(l)
    
    
  } else if(.approach=="SUMCORR"){
    
    # maximize the sum of the composite correlations sum(w_i%*%S_ii%*%t(w_i)), see Asendorf (2015, Appendix B.3.2 Empirical, p. 295)
    
    # define function to be minimized:
    # R: empirical correlation matrix
    # RDsqrt: Blockdiagonal matrix with S_{ii}^{-1/2} on its diagonal
    # nameInd: list with indicator names per block
    # name: vector with names of the composite
    fn = function(wtilde,R,RDsqrt,nameInd,nameLV) {
      -(t(wtilde)%*%RDsqrt%*%R%*%RDsqrt%*%t(t(wtilde)))
    }
    
    # constraints; variance of each composite is one
    heq <- function(wtilde,R,RDsqrt,nameInd,nameLV){
      names(wtilde)=unlist(nameInd)
      
      h <- rep(NA,length(nameLV))
      
      for(kk in 1:length(h)){
        h[kk]=t(wtilde[nameInd[[kk]]])%*%t(t(wtilde[nameInd[[kk]]])) -1 
      }
      h
    }
    
    # drawn random starting values for optimization
    p0=runif(length(unlist(namesIndicators)))
    
    # minimze
    res=alabama::auglag(par=p0, fn=fn, R=Srestructured, RDsqrt=H,nameInd=namesIndicators,nameLV=Construct_names ,heq=heq)
    wtilde=res$par
    
    names(wtilde)=unlist(namesIndicators)
    
    # Calculate weights w_i = R_{ii}^{-1/2} wtilde_i
    a=matrix(0,ncol=length(unlist(namesIndicators)),nrow=length(Construct_names),dimnames=list(Construct_names,unlist(namesIndicators)))
    
    for(i in Construct_names){
      a[i,namesIndicators[[i]]]=as.vector(Blocksqrtmcorrelation[[i]]%*%wtilde[namesIndicators[[i]]]   )
    }
    
    for(i in Construct_names){
      if(sum(a[i,])<0) a[i,]=-a[i,]
    }
    
    W = a[rownames(W),colnames(W)]
    
    W = scaleWeights(S, W)
    
    l <- list("W" = W, "E" = NULL, "Modes" = rep('SUMCORR',length(Construct_names)),
              "Conv_status" = res$convergence == 0, "Iterations" = res$counts[1])
    return(l)
    
  } else if(.approach == "GENVAR" ){
    
    # Function to be optimized
    fn = function(wtilde,R,RDsqrt,nameInd,nameLV) {
      det(-(t(wtilde)%*%RDsqrt%*%R%*%RDsqrt%*%t(t(wtilde))))
    } 
    
    # constraints; variance of each composite is one
    heq <- function(wtilde,R,RDsqrt,nameInd,nameLV){
      names(wtilde)=unlist(nameInd)
      
      h <- rep(NA,length(nameLV))
      
      # constraint: t(wtilde_i)%*%t(t(wtilde_i))=1
      for(kk in 1:length(h)){
        h[kk]=t(wtilde[nameInd[[kk]]])%*%t(t(wtilde[nameInd[[kk]]])) -1 
      }
      h
    }
    
    # drawn random starting values 
    p0=runif(length(unlist(namesIndicators)))
    
    # minimze
    res=alabama::auglag(par=p0, fn=fn, R=Srestructured, RDsqrt=H,nameInd=namesIndicators,nameLV=Construct_names ,heq=heq)
    
    wtilde=res$par
    
    names(wtilde)=unlist(namesIndicators)
    
    # Calculate weights w_i = R_{ii}^{-1/2} wtilde_i
    a=matrix(0,ncol=length(unlist(namesIndicators)),nrow=length(Construct_names),dimnames=list(Construct_names,unlist(namesIndicators)))
    
    for(i in Construct_names){
      a[i,namesIndicators[[i]]]=as.vector(Blocksqrtmcorrelation[[i]]%*%wtilde[namesIndicators[[i]]]   )
    }
    
    for(i in Construct_names){
      if(sum(a[i,])<0) a[i,]=-a[i,]
    }
    
    W=a[rownames(W),colnames(W)]
    
    W=scaleWeights(S, W)
    
    l <- list("W" = W, "E" = NULL, "Modes" = rep('GENVAR',length(Construct_names)),
              "Conv_status" = res$convergence == 0, "Iterations" = res$counts[1])
    return(l)
    
  } else if(.approach == "MINVAR"){
    
    stop("Not yet implemented")
  } else if(.approach == "SSQCORR"){
    
    # Function to be optimized
    fn = function(wtilde,R,RDsqrt,nameInd,nameLV) {
      norm(-(t(wtilde)%*%RDsqrt%*%R%*%RDsqrt%*%t(t(wtilde))),"F")
    }
    
    # constraints; variance of each composite is one
    
    heq <- function(wtilde,R,RDsqrt,nameInd,nameLV){
      names(wtilde)=unlist(nameInd)
      
      h <- rep(NA,length(nameLV))
      
      # constraint: t(wtilde_i)%*%t(t(wtilde_i))=1
      for(kk in 1:length(h)){
        h[kk]=t(wtilde[nameInd[[kk]]])%*%t(t(wtilde[nameInd[[kk]]])) -1 
      }
      h
    }
    
    # drawn random starting values 
    p0=runif(length(unlist(namesIndicators)))
    
    # minimze
    res=alabama::auglag(par=p0, fn=fn, R=Srestructured, RDsqrt=H,nameInd=namesIndicators,nameLV=Construct_names ,heq=heq)
    
    wtilde=res$par
    
    names(wtilde)=unlist(namesIndicators)
    
    # Calculate weights w_i = R_{ii}^{-1/2} wtilde_i
    a=matrix(0,ncol=length(unlist(namesIndicators)),nrow=length(Construct_names),dimnames=list(Construct_names,unlist(namesIndicators)))
    
    for(i in Construct_names){
      a[i,namesIndicators[[i]]]=as.vector(Blocksqrtmcorrelation[[i]]%*%wtilde[namesIndicators[[i]]]   )
    }
    
    for(i in Construct_names){
      if(sum(a[i,])<0) a[i,]=-a[i,]
    }
    
    W=a[rownames(W),colnames(W)] 
    
    W=scaleWeights(S, W)
    
    l <- list("W" = W, "E" = NULL, "Modes" = rep('SSQCORR',length(Construct_names)),
              "Conv_status" = res$convergence == 0, "Iterations" = res$counts[1])
    return(l)
    
  } else{stop("Not yet implemented")}
  
} # END calculateWeightsKettenring


#' Calculate weights using GSCA
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsGSCA(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsGSCA <- function(
  .data         = NULL,
  .model        = NULL
) {

  stop("Not yet implemented")

} # END calculateWeightsGSCA

#' Calculate weights using fixed weights
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsFixed(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsFixed <- function(
  .data         = NULL,
  .model        = NULL
) {

  stop("Not yet implemented")

} # END calculateWeightsFixed

#' Calculate unit weights for all blocks, i.e., each indicator of a block is equally weighted 
#' This approach might be useful for Kroon correction.
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsUnit(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'
calculateWeightsUnit = function(
  .data         = NULL,
  .model        = NULL
){
  
  ### Make the function more autonomous ========================================
  ## Convert model to "cSEMModel" format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    csem_model <- parseModel(.model)
  } else {
    csem_model <- .model
  }
  
  ## Prepare, standardize, check, and clean data if not already in this format
  if(!(class(.data) == "cSEMData")) {
    if(is.matrix(.data) && isSymmetric.matrix(.data)) {
      S <- .data
    } else {
      #  If function is used externally, data is not automatically scaled!
      X <- processData(.data = .data, .model = csem_model)
      S <- stats::cov(X)
    }
  }
  
  W=.model$measurement
  W=scaleWeights(S,W)
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = rep('unit',nrow(W)), "Conv_status" = TRUE,
            "Iterations" = 0)
  return(l)  
}

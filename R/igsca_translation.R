#' One-to-One R Translation of igsca_sim.m from Heungsun Hwang
#' 
#' @param Z0 Input data matrix
#' @param W0 Weight indicator: Indicators in rows, LVs in columns
#' @param C0 Loadings indicator. Indicators that load onto factors should have 1s in both W0 and C0.
#' @param B0 Path Coefficients Indicators
#' @param lv_type A boolean vector index for what latent variable is a factor. Length = ncol(W0)
#' @param ov_type A numeric vector index for what latent variable is a composite. Important for D matrix. Length = ncol(W0)
#' @param ind_domi A numeric vector index for the indicator that is dominant for each latent variable. Range of values should be the number of indicators, which is nrow(W0). Length should be ncol(W0)
#' @param nbt Number of boostraps -- though this should really be removed
#' @param testEquivalence TRUE/FALSE for whether comparison with Matlab should be made
#' @param swap_step List the stages where the R computation should be swapped for Matlab's to better identify the root of the divergences
#' @param itmax Maximum number of ALS iterations
#' @param ceps Minimum amount of absolute change in the estimates between ALS iterations before ending the optimization.
#' @author Michael S. Truong
#' 
#' @return
#' @export
#' @importFrom MASS ginv 
#' @importFrom R.matlab readMat
#' @importFrom here here
#' @examples
#' require(here)
#' require(readxl)
#' 
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' NumberofJobInterviews <~ Interview1 + Interview2
#' NumberofJobOffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' Honesty_Humility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' Openness_to_Experience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ Honesty_Humility + Emotionality + Extraversion + Agreeableness + Conscientiousness + Openness_to_Experience
#' NumberofJobInterviews ~ NetworkingBehavior
#' NumberofJobOffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' 
#' igsca_sim_in <- extract_parseModel(model = tutorial_igsca_model,
#'                                    data = dat,
#'                                    ind_domi_as_first = TRUE)
#'
#'
#' (igsca_sim_out <- igsca_sim(Z0 = igsca_sim_in$Z0, W0 = igsca_sim_in$W0, C0 = igsca_sim_in$C0,
#'                             B0 = igsca_sim_in$B0, lv_type = igsca_sim_in$lv_type,
#'                             ov_type = igsca_sim_in$ov_type, ind_domi = igsca_sim_in$ind_domi,
#'                             nbt = 0,
#'                             devmode = TRUE,
#'                             swap_step = "noswap")
#'                             )
igsca_sim <- function(Z0, W0, C0, B0, lv_type, ov_type, ind_domi, nbt,
                      swap_step = c("noswap", "prepare_for_ALS", "first_iteration_update",
                                    "first_factor_update", "first_composite_update",
                                    "flip_signs_ind_domi"),
                      itmax = 100, ceps = 0.001, devmode =  FALSE,
                      devdir = list("dev", "Notes", "data")) {
  
# Safety Checks ------------------------------------------------------------
  
  swap_step <- match.arg(swap_step)
  # TODO: This should be extended and completed to make igsca into a stand-alone function. However, cSEM already has its own safety checks that make this essentially unnecessary
  stopifnot("The number of bootstraps should be a non-negative integer" = nbt >= 0)

# Auxiliary Variables -----------------------------------------------------
  ncase <- nrow(Z0)
  nvar <- ncol(Z0)
  nlv <- ncol(W0)
  ntv <- nvar+ nlv
  
  windex <- which(c(W0) == 1) 
  cindex <- which(c(C0) == 1)
  bindex <- which(c(B0) == 1)
  

# Bootstrap Cycle ---------------------------------------------------------
for(nb in seq_len(nbt+1)) { 

## Initial Estimates and Preparation -------------------------------------
  prepared_for_ALS <- prepare_for_ALS(
    Z0 = Z0,
    W0 = W0,
    B0 = B0,
    nvar = nvar,
    ncase = ncase,
    nlv = nlv,
    ov_type = ov_type
  )
  ### Updates W, C, B, V, Z, D, U, Gamma
  list2env(prepared_for_ALS, envir = environment())
  
  if(isTRUE(devmode)) {
    # Takes entire environment and saves it as rds for later comparison
    as.list.environment(x = environment(), all.names = TRUE) |>
      saveRDS(file = here::here(devdir, "R_out", swap_step, "prepare_forALS.rds"))
    
    if (swap_step == "prepare_for_ALS") {
      swapdir <- c(devdir, "matlab_out", "prepare_for_ALS.MAT")
      R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
        list2env(envir = environment())
    } 
  }
  
## Alternating Least Squares Algorithm -------------------------------------

### While Counters and Initial Estimates ----------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (any(bindex)) {
      est <- B[bindex] 
    } else {
      est <- C[cindex]
    }    
    est0 <- est + 1 
    it <- 0
    
    if (isTRUE(devmode)) {
      composite_counter = 0
      factor_counter = 0
    }
### Alternate Between Updating Weights and Loadings -----------------------
    while ((sum(abs(est0 - est)) > ceps) && (it <= itmax)) {
      # TODO: Review better way to do this while condition

      # Counter Things
      it <- it + 1
      est0 <- est
      
#### Update Weights ----------------------------------------------------------
      # FIXME: warning("X and WW probably aren't weights, the naming for this section is a misnomer")
      updated_X_weights <- update_X_weights(
        Z = Z,
        U = U,
        D = D,
        C = C,
        nlv = nlv,
        B = B
      )
      # Creates X (for updating Composites) and WW (for updating Factors)
      list2env(updated_X_weights, envir = environment())
      
#### Iterative Update of LVs -------------------------------------------------
      A <- cbind(C, B) 
      
      for (j in seq_len(nlv)) {
      # After each cycle, the Gamma, W and V matrices are updated
        tot <- nvar + j
        windex_j <- (W0[, j] == 1)
        Xj <- X[, windex_j]
        
        
        if (lv_type[j] == 0) {
          # Update Composite LV
          theta <-
            update_composite_LV(
              ntv = ntv,
              tot = tot, # Changes per lv_update iteration
              nlv = nlv,
              j = j, # Changes per lv_update iteration
              W = W, # Changes per lv_update iteration
              A = A,
              V = V, # Changes per lv_update iteration
              X = X,
              windex_j = windex_j # Changes per lv_update iteration
            )
          
          
          if(exists("composite_counter") && isTRUE(devmode)) {
            
            if (composite_counter == 0) {
              
              as.list.environment(x = environment(), all.names = TRUE) |>
                saveRDS(file = here::here(devdir, "R_out", swap_step, "first_composite_update.rds"))

              composite_counter = composite_counter + 1
              
              if (swap_step == "first_composite_update") {
                swapdir <- c(devdir, "matlab_out", "first_composite_update.MAT")
                R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
                  list2env(envir = environment())
              } 
            }
          }
        } else if (lv_type[j] == 1) {
          # Update Factor LV
          theta <-
             update_factor_LV(
               WW = WW, 
               windex_j = windex_j, # Changes per iteration
               j = j # Changes per iteration
               )
          
          if(exists("factor_counter") && isTRUE(devmode)) {
            
            if (factor_counter == 0) {
              
              as.list.environment(x = environment(), all.names = TRUE) |>
                saveRDS(file = here::here(devdir, "R_out", swap_step, "first_factor_update.rds"))
              
              factor_counter = factor_counter + 1
              
              if (swap_step == "first_factor_update") {
                swapdir <- c(devdir, "matlab_out", "first_factor_update.MAT")
                R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
                  list2env(envir = environment())
                
              } 
            }
          }
        } else {
          stop("lv_type should either be 1 for factors or 0 for composites")
        }
        # This is where the 'actual' updating occurs, in terms of the Gamma matrix, Weights and V(?)
        # FIXME: warning("Is Matlab doing a 2-norm here?")
        theta <- theta / norm(Xj %*% theta, "2")
        
        Gamma[, j] <- Xj %*% theta
        W[windex_j, j] <- theta
        V[windex_j, tot] <- theta
        
        if(it == 1 && isTRUE(devmode)) {
          
          as.list.environment(x = environment(), all.names = TRUE) |>
            saveRDS(file = here::here(devdir, "R_out", swap_step, "first_iteration_update_theta_Gamma_W_V.rds"))
          
          if (swap_step == "first_iteration_update") {
            swapdir <- c(devdir, "matlab_out", "first_iteration_update_theta_Gamma_W_V.MAT")
            R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
              list2env(envir = environment())
          }
        }
    }
      
#### Update Loadings, Path Coefficients and Disturbance Terms ----------

      updated_C_B_D <- update_C_B_D(
        X = X,
        nvar = nvar,
        Gamma = Gamma,
        cindex = cindex,
        C = C,
        nlv = nlv,
        bindex = bindex,
        B = B,
        ncase = ncase,
        D = D,
        Z = Z,
        ov_type = ov_type
      )
      # Updates C, B, D, uniqueD, est, and U
      list2env(updated_C_B_D, envir = environment())
    }
    
    if(it == 1 && isTRUE(devmode)) {
      
      as.list.environment(x = environment(), all.names = TRUE) |>
        saveRDS(file = here::here(devdir, "R_out", swap_step, "first_iteration_update_C_B_D_uniqueD_est.rds"))
      
      if (swap_step == "first_iteration_update") {
        swapdir <- c(devdir, "matlab_out", "first_iteration_update_C_B_D_uniqueD_est.MAT")
        R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
          list2env(envir = environment())
      }
      }
    }

## Flip Signs for Factors and cOmposites Based on Dominant Indicators --------

    flipped_signs <-
      flip_signs_ind_domi(
        nlv = nlv,
        Z = Z,
        ind_domi = ind_domi,
        j = j,
        Gamma = Gamma,
        C = C,
        B = B
      )
    # Updates Gamma, C and B
    list2env(flipped_signs, envir = environment())
    
    if (isTRUE(devmode)) {
      
      as.list.environment(x = environment(), all.names = TRUE) |>
        saveRDS(file = here::here(devdir, "R_out", swap_step, "flip_signs_ind_domi.rds"))
      
      if (swap_step == "flip_signs_ind_domi") {
        swapdir <- c(devdir, "matlab_out", "flip_signs_ind_domi.MAT")
        R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
          list2env(envir = environment())
      }
    }
    
    
    mW <- W[windex] 
    # TODO: Commented out in matlab version, why?
    # mC <- c[cindex] 
    mC <- t(C)[t(C0) == 1] 
    mB <- B[bindex] 
    mD <- uniqueD
  
  
  if (isTRUE(devmode)) {
    
    as.list.environment(x = environment(), all.names = TRUE) |>
      saveRDS(file = here::here(devdir, "R_out", swap_step, "end_results.rds"))
    
    # Heuristic check for when a function is returning more arguments than what is needed
    # Although in-order to maintain testability against matlab, this may have to be kept...
    # TODO: This could be improved later
    # devcheck <- print_important_objects()
    # 
    # extra_arguments <- lapply(list("prepared_for_ALS" = prepared_for_ALS,
    #                                "updated_X_weights" = updated_X_weights,
    #                                "updated_C_B_D" = updated_C_B_D),
    #                           FUN = \(x)  names(x)[!(names(x) %in% devcheck)])
    # extra_arguments 
    
  }
    return(
      list(
        "Weights (mW)" = mW,
        "Loadings (mC)" =  mC,
        "Path Coefficients (mB)" =  mB,
        "Uniqueness Terms (mD)" = mD
      )
    )

}


#' A parseModel extractor function for the purposes of running I-GSCA code example
#'
#' @param model Specified Model in lavaan style
#' @param data Dataframe 
#' @param ind_domi_as_first Boolean for whether the first indicator for each latent factor should be chosen as the dominant indicator
#'
#' @return Z0, W0, B0, C0, lv_type, ov_type, ind_domi
#' @export
#' @importFrom csem parseModel
#' @examples
#' 
#' require(here)
#' require(readxl)
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' NumberofJobInterviews <~ Interview1 + Interview2
#' NumberofJobOffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' Honesty_Humility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' Openness_to_Experience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ Honesty_Humility + Emotionality + Extraversion + Agreeableness + Conscientiousness + Openness_to_Experience
#' NumberofJobInterviews ~ NetworkingBehavior
#' NumberofJobOffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' extract_parseModel(model = tutorial_igsca_model, data = dat, ind_domi_as_first = TRUE)
extract_parseModel <- function(model, data, ind_domi_as_first = TRUE) {
  
  # TODO: Probably shouldn't use cSEM:: within the cSEM package
  csemify <- cSEM::parseModel(.model = model)
  
  Z0 <- data[, csemify$indicators]
  B0 <- csemify$structural
  W0 <- t(csemify$measurement)
  lv_type <- csemify$construct_type == "Common factor"
  ov_type <- csemify$construct_type == "Composite"
  
  # TODO: More parsimonious expression to get to C0
  W0_to_C0 <- W0
  W0_to_C0[, ov_type] <- 0
  C0 <- W0_to_C0
  
  stopifnot(
    "Data matrix (Z0) does not correctly correspond with Weights matrix (W0)" = identical(colnames(Z0), rownames(W0)),
    "Weights matrix (W0) does not correctly correspond with Structural matrix (B0)" = identical(colnames(W0), colnames(B0)),
    "Construct indicator does not correctly correspond with Weights Matrix (W0)" = identical(names(csemify$construct_type), colnames(W0))
  )
  
  if(isTRUE(ind_domi_as_first)) {
  # Row Indices of W0 that correspond to the dominant indicator for each factor/composite
  # TODO: Change to W0 to W0[, lv_type] if it turns out that the correction should only apply to latent factors
  ind_domi <- apply(W0, 2, as.logical) |> 
    apply(2, which, TRUE) |>
    lapply(FUN = \(x) x[[1]]) |>
    unlist()
  } else {
      ind_domi <- NA
    }
  
  return(list(
    "Z0" = Z0,
    "B0" = B0,
    "W0" = W0,
    "lv_type" = lv_type,
    "ov_type" = ov_type,
    "C0" = C0,
    "ind_domi" = ind_domi
  ))
}

#' Writes the extracted matrices to run with igsca_sim_test.m
#'
#' @param extracted_matrices Object returned by extract_parseModel
#'
#' @return
#' @export
#' @importFrom here here
#' @examples
#' 
#' require(here)
#' require(readxl)
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' NumberofJobInterviews <~ Interview1 + Interview2
#' NumberofJobOffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' Honesty_Humility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' Openness_to_Experience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ Honesty_Humility + Emotionality + Extraversion + Agreeableness + Conscientiousness + Openness_to_Experience
#' NumberofJobInterviews ~ NetworkingBehavior
#' NumberofJobOffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' write_for_matlab(extract_parseModel(model = tutorial_igsca_model, data = dat, ind_domi_as_first = TRUE))
write_for_matlab <- function(extracted_matrices) {
  
  indir <- list("dev", "Notes", "data", "matlab_in")
  extracted_matrices$lv_type <- as.numeric(extracted_matrices$lv_type)
  extracted_matrices$ov_type <- as.numeric(extracted_matrices$ov_type)
  
  mapply(
    write.csv,
    x = extracted_matrices,
    file = paste0(here::here(indir), "/", names(extracted_matrices), ".csv"),
    row.names = FALSE
  )
  
  invisible()  
}


#' One-to-One R Translation of gsca_inione.m from Heungsun Hwang
#' 
#' Initializes the values for I-GSCA
#' @param Z0 matrix...
#' @param W0 Weights
#' @param B0 Path Coefficients
#' 
#' @author Michael S. Truong
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' require(here)
#' require(readxl)
#' 
#' tutorial_gsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' NumberofJobInterviews <~ Interview1 + Interview2
#' NumberofJobOffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model Forced into Composite
#' Honesty_Humility <~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality <~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion <~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness <~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness <~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' Openness_to_Experience <~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ Honesty_Humility + Emotionality + Extraversion + Agreeableness + Conscientiousness + Openness_to_Experience
#' NumberofJobInterviews ~ NetworkingBehavior
#' NumberofJobOffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' 
#' gsca_sim_in <- extract_parseModel(model = tutorial_gsca_model,
#'                                    data = dat,
#'                                    ind_domi_as_first = TRUE)
#'                          
#' (gsca_inione_test <- gsca_inione(Z0 = gsca_sim_in$Z0, W0 = gsca_sim_in$W0, B0 = gsca_sim_in$B0))
gsca_inione <- function(Z0, W0, B0) {
  N <- nrow(Z0)
  J <- nrow(W0)
  P <- ncol(W0)
  TRep <- J + P 
  C0 <- t(W0)
  A0 <- cbind(C0, B0)
  
  windex <- which(W0 != 0)
  aindex <- which(A0 != 0)
  
  Z <- scale(Z0, center = TRUE, scale = TRUE) * sqrt(N) / sqrt(N-1)
  # Random Values to W and A
  W <- W0
  A <- A0
  
  W[windex] <- rep(1, length(windex)) 
  A[aindex] <- runif(length(aindex), min = 0, max = 1)
  
  Gamma <- Z %*% W
  W <- W / (t(sqrt(diag(t(Gamma) %*% Gamma))) |> rep(each = nrow(W))) 
  V <- cbind(diag(J), W) 
  Gamma <- Z %*% W
  
  Psi <- Z %*% V
  vecPsi <- c(Psi) 
  it <- 0
  f0 <- 100000000
  imp <- 100000000
  while ((it <= 300) && (imp > 0.0001)) {
   it <- it + 1
   Phi = kronecker(diag(TRep), Gamma)
   Phi <- Phi[, aindex]
   A[aindex] <- solve(t(Phi) %*% Phi, t(Phi)) %*% vecPsi
   
   for (p in seq_len(P)) {
     t_lil <- J + p
     windex_p <- which(W0[, p] != 0)
     m <- matrix(0, nrow = 1, ncol = TRep)
     m[t_lil] <- 1
     a <- A[p, ]
     beta <- m - a
     H1 <- diag(P)
     H2 <- diag(TRep)
     H1[p, p] <- 0
     H2[t_lil, t_lil] <- 0
     Delta <- (W %*% H1 %*% A) - (V %*% H2)
     vecZDelta <- c(Z %*% Delta)
     
     XI <- kronecker(t(beta), Z)
     XI <- XI[, windex_p]
     theta = MASS::ginv(t(XI) %*% XI) %*% t(XI) %*% vecZDelta
     zw = Z[,windex_p]%*%theta;
     
     theta <- sqrt(N) * theta / norm(zw, "2") 
     W[windex_p, p] <- theta
     V[windex_p, t_lil] <- theta
   }
   Gamma <- Z %*% W
   Psi <- Z %*% V
   dif <- Psi - Gamma %*% A
   f <- sum(diag(t(dif) %*% dif))
   imp <- f0 - f
   f0 <- f
   vecPsi <- c(Psi)
  }
  C <- A[, 1:J]
  B <- A[, (J+1) : ncol(A)]
  return(list("W" = W, "C" = C, "B" = B, "it" = it))
}

#' Flip signs of Gamma, Loadings and Path-Coefficients Cells Based on Dominant Indicator
#'
#' @param nlv 
#' @param Z 
#' @param ind_domi 
#' @param j 
#' @param Gamma 
#' @param C 
#' @param B 
#'
#' @return Gamma 
#' @export
#'
#' @examples
flip_signs_ind_domi <- function(nlv, Z, ind_domi, j, Gamma, C, B) {

  for (j in seq_len(nlv)) {
    if ((t(Z[, ind_domi[j]]) %*% Gamma[, j]) < 0) {
      Gamma[, j] <- (-1 * Gamma[, j])
      C[j,] <- (-1 * C[j,])
      B[B,] <- (-1 * B[j,])
      B[, j] <- (-1 * B[, j])
    }
  }
  return(list("Gamma" = Gamma, "C" = C, "B" = B))
}


#' Prepare for ALS Algorithm
#'
#' Internal I-GSCA function
#'
#' @param Z0 
#' @param W0 
#' @param B0 
#' @param nvar 
#' @param ncase 
#' @param nlv 
#' @param ov_type 
#'
#' @return
#' @export
#'
#' @examples
prepare_for_ALS <- function(Z0, W0, B0, nvar, ncase, nlv, ov_type) {
  
  bZ0 <- Z0
  
  
  initial_est <-
    gsca_inione(Z0 = bZ0,
                W0 = apply(W0 != 0, 2, as.numeric), 
                B0 = apply(B0 != 0, 2, as.numeric)
    ) 
  W <- initial_est$W
  C <- initial_est$C
  B <- initial_est$B
  
  
  V <- cbind(diag(nvar), W) 
  Z <- scale(bZ0, center = TRUE, scale = TRUE) / sqrt(ncase - 1) 
  Gamma <- Z %*% W 
  D <- diag(nvar)
  
  
  # TODO: Does the LAPACK argument for qr() still matter? Why does `complete` matter? 
  # Solution for Q is copied from estimators_weights.R
  Q <- qr.Q(qr(Gamma), complete =  TRUE) 
  F_o <- Q[, (nlv+1):ncase, drop = FALSE] 
  # FIXME: Should compare the svd of matlab and R, unsure if it matters that they don't match
  # FIXME: Come back to this stack overflow thing for proper citation
  # Ahmed Fasih https://stackoverflow.com/a/41972818
  svd_out <- (D %*% t(Z) %*% F_o) |>
    {\(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))}()
  u <- svd_out$u
  v <- svd_out$v
  
  # TODO: Utilde deviates from Matlab because of the SVD
  Utilde <- v[, 1:nvar] %*% t(u) 
  U <- F_o %*% Utilde 
  D <- diag(diag(t(U) %*% Z)) 
  D[ov_type == 0, ov_type == 0] <- 0
  
  
  return(
    list(
      "W" = W,
      "C" = C,
      "B" = B,
      "V" = V,
      "Z" = Z,
      "D" = D,
      "U" = U,
      "Gamma" = Gamma
      # "Utilde" = Utilde,
      # "v" = v,
      # "u" = u,
      # "F_o" = F_o,
      # "Q" = Q,
    )
  )
}

#' Update Loadings, Path-Coefficients and Uniqueness Terms After Updating Latent Variables
#'
#' @param X 
#' @param nvar 
#' @param Gamma 
#' @param cindex 
#' @param C 
#' @param nlv 
#' @param bindex 
#' @param B 
#' @param ncase 
#' @param D 
#' @param Z 
#' @param ov_type 
#'
#' @return
#' @export
#'
#' @examples
update_C_B_D <- function(X, nvar, Gamma, cindex, C, nlv, bindex, B, ncase, D, Z, ov_type) {
  
  t1 <- c(X)
  M1 <- kronecker(diag(nvar), Gamma) 
  M1 <- M1[, cindex]
  C[cindex] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)   
  
  t2 <- c(Gamma)
  M2 <- kronecker(diag(nlv), Gamma)
  M2 <- M2[, bindex]
  B[bindex] <- MASS::ginv(t(M2) %*% M2) %*% (t(M2) %*% t2)
  
  # Solution for Q is copied from estimators_weights.R
  Q <- qr.Q(qr(Gamma), complete =  TRUE) 
  F_o <- Q[, (nlv+1):ncase]
  # FIXME: warning("It's unclear to me whether this SVD is safe. The Matlab code seems to either do a 'normal' svd or a economy svd depending on the dimensionality of the input matrix. Should look into this more.")
  # https://www.mathworks.com/help/matlab/ref/double.svd.html?searchHighlight=svd&s_tid=srchtitle_support_results_1_svd#d126e1597915
  svd_out2 <- (D %*% t(Z) %*% F_o) |>
    {\(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))}()
  u <- svd_out2$u
  v <- svd_out2$v
  Utilde <- v[, 1:nvar] %*% t(u)
  U <- F_o %*% Utilde
  D <- diag(diag(t(U) %*% Z))
  D[ov_type != 1, ov_type != 1] <- 0
  
  
  if (any(bindex)) {
    est <- B[bindex] 
  } else {
    est <- C[cindex]
  }
  
  uniqueD <- diag(D) ^ 2
  
  return(
    list(
      "D" = D,
      "uniqueD" = uniqueD,
      "est" = est,
      "U" = U,
      "Utilde" = Utilde,
      "v" = v,
      "u" = u,
      "F_o" = F_o,
      "Q" = Q,
      "B" = B,
      "C" = C
    )
  )
}

#' Update Weights and X (?)
#'
#' @param Z 
#' @param U 
#' @param D 
#' @param C 
#' @param nlv 
#' @param B 
#'
#' @return
#' @export
#'
#' @examples
update_X_weights <- function(Z, U, D, C, nlv, B) {
  # X deviates from Matlab because it is an offspring of svd_out
  X <- Z - U %*% D
  # FIXME: warning("I don't quite understand why this is the equivalent of the Matlab expression, revisit page 44 of Hiebeler 2015 R and Matlab")
  # TODO: Should find alternative to solve() to avoid matrix inversions
  WW <-
    t(C) %*% solve((C %*% t(C) + diag(nlv) - 2 * B + (B %*% t(B))))
  return(list("X" = X, "WW" = WW))
}

#' Update Composite Latent Variables
#'
#' @param ntv 
#' @param tot 
#' @param nlv 
#' @param j 
#' @param W 
#' @param A 
#' @param V 
#' @param X 
#' @param windex_j 
#'
#' @return
#' @export
#'
#' @examples
update_composite_LV <- function(ntv, tot, nlv, j, W, A, V, X, windex_j) {
  
  
  # This updates the composite
  e <- matrix(0, nrow = 1, ncol = ntv)
  e[tot] <- 1 
  H1 <- diag(ntv)
  H2 <- diag(nlv)
  H1[tot, tot] <- 0
  H2[j, j] <- 0
  Delta <- (W %*% H2 %*% A) - (V %*% H1)
  
  
  vecZDelta <- c(X %*% Delta)
  beta <- e - A[j,]
  XI <- kronecker(t(beta), X)
  XI <- XI[, windex_j]
  
  theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta
  return(theta)
}


#' Update Factor Latent Variable
#'
#' @param WW 
#' @param windex_j 
#' @param j 
#'
#' @return
#' @export
#'
#' @examples
update_factor_LV <- function(WW, windex_j, j) {
  
  theta <- WW[windex_j, j]
  return(theta)
}


#' List arguments and objects that are important for fitting
#' 
#' Internal development function 
#' 
#' Anything returned in this list should always be returned by other functions so that they are up-to-date
#' 
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' devcheck <- print_important_objects()
print_important_objects <- function() {
  
  objects_used_by_functions <- lapply(list(extract_parseModel,
              flip_signs_ind_domi,
              gsca_inione,
              igsca_sim,
              prepare_for_ALS,
              update_C_B_D,
              update_composite_LV,
              update_factor_LV,
              update_X_weights,
              write_for_matlab), formals) |>
    unlist() |>
    names() |>
    unique()
  return(
    c(
      objects_used_by_functions,
      "uniqueD",
      "D",
      "est",
      "U",
      "B",
      "W",
      "Gamma",
      "C",
      "V",
      "X",
      "theta"
    ) |> unique()
  )
}


#' Read csv File and Convert to Matrix
#' 
#' Small helper function to read csv files and convert them to matrices immediately. This is to ease interoperability between matlab and R.
#' 
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
read.csv.to.matrix <- function(file) {
  
  read.csv(file) |>
    as.matrix()
}

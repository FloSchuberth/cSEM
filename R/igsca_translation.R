#' One-to-One R Translation of igsca_sim.m by Heungsun Hwang
#' 
#' This is a one-to-one R translation of the matlab implementation of I-GSCA by Heungsun Hwang. This implementation was also re-factored for ease of interpretation. 
#' 
#' The code here generally assumes that every indicator corresponds to a single latent variable.
#' 
#' @param z0 Input data-frame/matrix. Should contain named columns.
#' @param W0 Indicator matrix of weights: indicators (rows) and their corresponding latent variable (columns).
#' @param C0 Indicator matrix of loadings: indicators (rows) and their corresponding latent variable (columns).
#' @param B0 Square indicator matrix of path coefficients: from-latent-variable (rows) and to-latent-variable (columns). The order of latent variables should match the order in C0 and W0.
#' @param lv_type A logical vector that indices whether a latent variable (columns in W0 and C0) is a factor (TRUE) or composite (FALSE). Its length should be equal to the number of columns of W0 and C0. 
#' 
#' @param ov_type An indicator vector that indices whether an indicator (rows of W0 and C0) corresponds to a factor latent variable (1) or a composite latent variable (0). This vector is important for computing the uniqueness terms (D) because it zeros the entries for composite indicators. 
#' 
#' @param ind_domi A numeric vector that indices the dominant indicator for each latent variable. *It is to be clarified whether this should only apply to factor latent variables or also composite latent variables.* This is important for ensuring that the signs of the path-coefficients and loadings are consistent. It is sometimes the case in composite-based structural equation modelling methods that loadings/path-coefficients may have the opposite sign. The length of this vector should be equal to the number of latent variables and each value should represent the row number of the dominant indicator for that latent variable. 
#' 
#' @param nbt Non-functional argument for the number of bootstrap-iterations, intended for computation of standard-errors, 95% confidence intervals and p-values. 
#' @param devmode Logical argument (TRUE/FALSE) for whether devmode functionality should be enabled. When enabled, the function saves the function's environment as a .Rds at various stages of the algorithm. This was intended to facilitate the comparison between the R and Matlab implementations by narrowing down sources of discrepancy.
#' @param swap_step This is an argument for different key stages of the algorithm. When selected, it swaps the overlapping objects between the R and Matlab environment at that stage of the algorithm. Requires devmode to be TRUE to be functional. 
#' *'noswap'* means that no swapping will be made.
#' 
#' @param itmax Maximum number of iterations in the Alternating Least Squares (ALS) algorithm.
#' 
#' @param ceps Minimum amount of absolute change in the estimates of the path-coefficients (if B0 is non-zero) or the loadings (if B0 is all zero, meaning ther are no path-coefficients) between ALS iterations before ending the optimization.
#' 
#' @param devdir A list of folders to be provided to here::here() to store the saved environments from `swap_step`.
#' @param devmode_checkobj A development argument for checking whether the list of overlapping and non-overlapping objects between the R and Matlab implementations of I-GSCA are the same as before. 
#' @author Michael S. Truong
#' 
#' @returns List of 4 matrices that make up a fitted I-GSCA Model: (1) Weights, (2) Loadings, (3) Uniqueness Terms D^2, and (4) Path Coefficients.
#' @export
#' @importFrom MASS ginv 
#' @importFrom R.matlab readMat
#' @importFrom here here
#' @examples
#' require(here)
#' require(readxl)
#' 
#' # Specify the model according to GSCA Pro's example
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' # Load the data
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' # Pre-generate the required matrices for the algorithm and assume that the first indicator for each latent variable is the dominant indicator.
#' igsca_sim_in <- extract_parseModel(model = tutorial_igsca_model,
#'                                    data = dat,
#'                                    ind_domi_as_first = TRUE)
#'
#' # Fit the I-GSCA model
#' (igsca_sim_out <- igsca_sim(z0 = igsca_sim_in$z0, W0 = igsca_sim_in$W0, C0 = igsca_sim_in$C0,
#'                             B0 = igsca_sim_in$B0, lv_type = igsca_sim_in$lv_type,
#'                             ov_type = igsca_sim_in$ov_type, ind_domi = igsca_sim_in$ind_domi,
#'                             nbt = 0,
#'                             devmode = FALSE,
#'                             swap_step = "noswap")
#'                             )
igsca_sim <-
  function(z0,
           W0,
           C0,
           B0,
           lv_type,
           ov_type,
           ind_domi,
           nbt,
           swap_step = c(
             "noswap",
             "prepare_for_ALS",
             "first_iteration_update_theta_Gamma_W_V",
             "first_iteration_update_C_B_D_uniqueD_est",
             "first_factor_update",
             "first_composite_update",
             "flip_signs_ind_domi"
           ),
           itmax = 100,
           ceps = 0.001,
           devmode =  FALSE,
           devdir = list("tests", "comparisons", "igsca_translation"),
           devmode_checkobj = FALSE) {
    
  
# Safety Checks ------------------------------------------------------------
  
  swap_step <- match.arg(swap_step)
  # TODO: This could be extended and completed to make igsca into a stand-alone function. However, cSEM already has its own safety checks that make this essentially unnecessary
  stopifnot("The number of bootstraps should be a non-negative integer" = nbt >= 0)
  
# Auxiliary Variables -----------------------------------------------------
  ncase <- nrow(z0)
  nvar <- ncol(z0)
  nlv <- ncol(W0)
  ntv <- nvar + nlv
  
  windex <- which(c(W0) == 1) 
  cindex <- which(c(C0) == 1)
  bindex <- which(c(B0) == 1)
  
  # column and row names are kept aside because swapping with matlab's environment removes the column and row names. Arrays in matlabs don't have row or column names.
  special_names <- list()
  special_names$weights$row <- rownames(W0)
  special_names$weights$col <- colnames(W0)
  special_names$loadings$row <- rownames(C0)
  special_names$loadings$col <- colnames(C0)
  special_names$path$row <- rownames(B0)
  special_names$path$col <- colnames(B0)
  

# Bootstrap Cycle ---------------------------------------------------------
for(nb in seq_len(nbt+1)) { 

## Initial Estimates and Preparation -------------------------------------
  prepared_for_ALS <- prepare_for_ALS(
    z0 = z0,
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
      
      # Matlab vectors are read-back as vectors, which is undesirable for R
      # So here we record vectorness before the swap and re-apply it post-swap
      # using `convert_matlab2R()`
      vectorness <-
        lapply(as.list.environment(x = environment(), all.names = TRUE),
               is.vector)
      
      # Only take the objects that have overlapping names between R and Matlab
      mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
      mat_env <- mat_env[names(mat_env) %in% ls()] |>
        convert_matlab2R(vectorness = vectorness) 
      list2env(x = mat_env, envir = environment())
      
    } 
  }
  
## Alternating Least Squares Algorithm -------------------------------------

### While Counters and Initial Estimates ----------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (length(bindex) > 0) {
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
                
                vectorness <-
                  lapply(as.list.environment(x = environment(), all.names = TRUE),
                         is.vector)
                
                # Only take the objects that have overlapping names between R and Matlab
                mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
                mat_env <- mat_env[names(mat_env) %in% ls()] |>
                  convert_matlab2R(vectorness = vectorness) 
                list2env(x = mat_env, envir = environment())
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
                
                vectorness <-
                  lapply(as.list.environment(x = environment(), all.names = TRUE),
                         is.vector)
                
                # Only take the objects that have overlapping names between R and Matlab
                mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
                mat_env <- mat_env[names(mat_env) %in% ls()] |>
                  convert_matlab2R(vectorness = vectorness) 
                list2env(x = mat_env, envir = environment())
                
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
          
          if (swap_step == "first_iteration_update_theta_Gamma_W_V") {
            swapdir <- c(devdir, "matlab_out", "first_iteration_update_theta_Gamma_W_V.MAT")
            
            vectorness <-
              lapply(as.list.environment(x = environment(), all.names = TRUE),
                     is.vector)
            
            # Only take the objects that have overlapping names between R and Matlab
            mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
            mat_env <- mat_env[names(mat_env) %in% ls()] |>
              convert_matlab2R(vectorness = vectorness) 
            list2env(x = mat_env, envir = environment())
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
      
      if (swap_step == "first_iteration_update_C_B_D_uniqueD_est") {
        swapdir <- c(devdir, "matlab_out", "first_iteration_update_C_B_D_uniqueD_est.MAT")
        
        vectorness <-
          lapply(as.list.environment(x = environment(), all.names = TRUE),
                 is.vector)
        
        # Only take the objects that have overlapping names between R and Matlab
        mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
        mat_env <- mat_env[names(mat_env) %in% ls()] |>
          convert_matlab2R(vectorness = vectorness) 
        list2env(x = mat_env, envir = environment())
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
        
        vectorness <-
          lapply(as.list.environment(x = environment(), all.names = TRUE),
                 is.vector)
        
        # Only take the objects that have overlapping names between R and Matlab
        # Note: Previous iteration of this put ls() in an anonymous function, which does not work. 
        # For example, this {\(x) ls()}() only shows x 
        
        mat_env <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) 
        mat_env <- mat_env[names(mat_env) %in% ls()] |>
          convert_matlab2R(vectorness = vectorness) 
        list2env(x = mat_env, envir = environment())
      }
    }
    
    
    mW <- W[windex] 
    # This was commented out in the Matlab version because the loadings dimensions
    # are transposed relative to the weights matrix
    # mC <- c[cindex] 
    mC <- t(C)[t(C0) == 1] 
    mB <- B[bindex] 
    mD <- uniqueD
  
  
  # Format Output
  Weights <- W
  rownames(Weights) <- special_names$weights$row
  colnames(Weights) <- special_names$weights$col
  
  Loadings <- t(C)
  rownames(Loadings) <- special_names$loadings$row
  colnames(Loadings) <- special_names$loadings$col
  
  PathCoefficients <- B
  rownames(PathCoefficients) <- special_names$path$row
  colnames(PathCoefficients) <- special_names$path$col
  
  # browser()
  
  names(uniqueD) <- rownames(Loadings)
  
  
  if (isTRUE(devmode)) {
    as.list.environment(x = environment(), all.names = TRUE) |>
      saveRDS(file = here::here(devdir, "R_out", swap_step, "end_results.rds"))
    
    if (isTRUE(devmode_checkobj)) {
    ### Checking that objects in R and Matlab are named similarly
    R_names <- ls()
    swapdir <- c(devdir, "matlab_out", "end_results.MAT")
    mat_names <- R.matlab::readMat(here::here(swapdir), fixNames = FALSE) |>
      names()
    
    # browser()
    
    non_overlapping_names <- table(c(R_names, mat_names)) |>
      {\(x) x[x < 2]}() |>
      names()
    # Objects in R that don't appear in Matlab
    RObj_NotIn_Matlab <- R_names[R_names %in% non_overlapping_names]
    # Objects in Matlab that don't appear in R
    MatObj_NotIn_R <- mat_names[mat_names %in% non_overlapping_names]
    
    # For Maintenance
    # "As manually checked on December 17/2023, I'm fairly sure that the non-overlapping names are OK.
    # Most of the Matlab 'unique' names actually show-up in the R sub-functions.
    # A small minority simply don't appear for parsimony (C_t, crindex) and some
    # probably aren't legal/good-idea R names (t). The R-unique object names are
    # intermediaries from functions"
    # For Maintenance 
    
    return(list("RObj_NotIn_Matlab" = RObj_NotIn_Matlab, "MatObj_NotIn_R" =  MatObj_NotIn_R))
   } 
  }
    return(
      list(
        "Weights" = Weights,
        "Loadings" =  Loadings,
        "Path Coefficients" =  PathCoefficients,
        "Uniqueness Terms" = uniqueD
      )
    )

}


#' A parseModel extractor function for the purposes of running I-GSCA code example
#'
#' @param model Specified Model in lavaan style
#' @param data Dataframe 
#' @param ind_domi_as_first Boolean for whether the first indicator for each latent factor should be chosen as the dominant indicator
#'
#' @return z0, W0, B0, C0, lv_type, ov_type, ind_domi
#' @export
#' @examples
#' 
#' require(here)
#' require(readxl)
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' extract_parseModel(model = tutorial_igsca_model, data = dat, ind_domi_as_first = TRUE)
extract_parseModel <-
  function(model, data, ind_domi_as_first = TRUE) {
    # Note: parseModel is from cSEM internal
    csemify <- parseModel(.model = model)
    
    z0 <- data[, csemify$indicators]
    # browser()
    # B0 <- csemify$structural
    B0 <- t(csemify$structural)
    W0 <- t(csemify$measurement)
    
    lv_type <- csemify$construct_type == "Common factor"
    
    # Constructing ov_type
    ov_type <- vector(length = ncol(csemify$measurement))
    names(ov_type) <- colnames(csemify$measurement)
      
    
    for (lv in rownames(csemify$measurement)) {
      for (indicator in colnames(csemify$measurement)) {
        if (csemify$construct_type[lv] == "Common factor") {
        
          if (csemify$measurement[lv, indicator] == 1) {
            ov_type[indicator] <- TRUE
          }
          
        }
      }
    }
    
    C0 <- W0
    
    
    stopifnot(
      "Data matrix (z0) does not correctly correspond with Weights matrix (W0)" = identical(colnames(z0), rownames(W0)),
      "Weights matrix (W0) does not correctly correspond with Structural matrix (B0)" = identical(colnames(W0), colnames(B0)),
      "Construct indicator does not correctly correspond with Weights Matrix (W0)" = identical(names(csemify$construct_type), colnames(W0))
    )
    
    if (isTRUE(ind_domi_as_first)) {
      # Row Indices of W0 that correspond to the dominant indicator for each factor/composite
      # TODO: Change to W0 to W0[, lv_type] if it turns out that the correction should only apply to latent factors
      ind_domi <- apply(W0, 2, as.logical) |>
        apply(2, which, TRUE) |>
        lapply(FUN = \(x) x[[1]]) |>
        unlist()
    } else {
      ind_domi <- NA
    }
    
    return(
      list(
        "z0" = z0,
        "B0" = B0,
        "W0" = W0,
        "lv_type" = lv_type,
        "ov_type" = ov_type,
        "C0" = C0,
        "ind_domi" = ind_domi
      )
    )
  }

#' Writes the extracted matrices to run with igsca_sim_test.m
#'
#' @param extracted_matrices Object returned by extract_parseModel
#'
#' @return
#' @importFrom here here
#' @examples
#' 
#' require(here)
#' require(readxl)
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' write_for_matlab(extract_parseModel(model = tutorial_igsca_model, data = dat, ind_domi_as_first = TRUE))
write_for_matlab <- function(extracted_matrices) {
  indir <- list("tests", "comparisons", "igsca_translation", "matlab_in")
  extracted_matrices$lv_type <-
    as.numeric(extracted_matrices$lv_type)
  extracted_matrices$ov_type <-
    as.numeric(extracted_matrices$ov_type)
  
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
#' @param z0 matrix...
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
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model Forced into Composite
#' HonestyHumility <~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality <~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion <~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness <~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness <~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience <~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx")) 
#' 
#' 
#' gsca_sim_in <- extract_parseModel(model = tutorial_gsca_model,
#'                                    data = dat,
#'                                    ind_domi_as_first = TRUE)
#'                          
#' (gsca_inione_test <- gsca_inione(z0 = gsca_sim_in$z0, W0 = gsca_sim_in$W0, B0 = gsca_sim_in$B0))
gsca_inione <- function(z0, W0, B0) {
  N <- nrow(z0)
  J <- nrow(W0)
  P <- ncol(W0)
  TRep <- J + P
  C0 <- t(W0)
  A0 <- cbind(C0, B0)
  
  windex <- which(W0 != 0)
  aindex <- which(A0 != 0)
  
  Z <- scale(z0, center = TRUE, scale = TRUE) * sqrt(N) / sqrt(N - 1)
  # Random Values to W and A
  W <- W0
  A <- A0
  
  W[windex] <- rep(1, length(windex))
  A[aindex] <- runif(length(aindex), min = 0, max = 1)
  
  Gamma <- Z %*% W
  W <-
    W / (t(sqrt(diag(t(
      Gamma
    ) %*% Gamma))) |> rep(each = nrow(W)))
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
      a <- A[p,]
      beta <- m - a
      H1 <- diag(P)
      H2 <- diag(TRep)
      H1[p, p] <- 0
      H2[t_lil, t_lil] <- 0
      Delta <- (W %*% H1 %*% A) - (V %*% H2)
      vecZDelta <- c(Z %*% Delta)
      
      XI <- kronecker(t(beta), Z)
      XI <- XI[, windex_p]
      theta <- MASS::ginv(t(XI) %*% XI) %*% t(XI) %*% vecZDelta
      zw <- Z[, windex_p] %*% theta
      
      
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
  B <- A[, (J + 1):ncol(A)]
  return(list(
    "W" = W,
    "C" = C,
    "B" = B,
    "it" = it
  ))
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
#'
#' @examples
flip_signs_ind_domi <- function(nlv, Z, ind_domi, j, Gamma, C, B) {
  for (j in seq_len(nlv)) {
    if ((t(Z[, ind_domi[j]]) %*% Gamma[, j]) < 0) {
      Gamma[, j] <- (-1 * Gamma[, j])
      C[j, ] <- (-1 * C[j, ])
      B[B, ] <- (-1 * B[j, ])
      B[, j] <- (-1 * B[, j])
    }
  }
  return(list(
    "Gamma" = Gamma,
    "C" = C,
    "B" = B
  ))
}


#' Prepare for ALS Algorithm
#'
#' Internal I-GSCA function
#'
#' @param z0 
#' @param W0 
#' @param B0 
#' @param nvar 
#' @param ncase 
#' @param nlv 
#' @param ov_type 
#'
#' @return
#'
#' @examples
prepare_for_ALS <- function(z0, W0, B0, nvar, ncase, nlv, ov_type) {
  bz0 <- z0
  
  
  initial_est <-
    gsca_inione(
      z0 = bz0,
      W0 = apply(W0 != 0, 2, as.numeric),
      B0 = apply(B0 != 0, 2, as.numeric)
    )
  W <- initial_est$W
  C <- initial_est$C
  B <- initial_est$B
  
  
  V <- cbind(diag(nvar), W)
  Z <- scale(bz0, center = TRUE, scale = TRUE) / sqrt(ncase - 1)
  Gamma <- Z %*% W
  D <- diag(nvar)
  
  
  # TODO: Does the LAPACK argument for qr() still matter? Why does `complete` matter?
  # Solution for Q is copied from estimators_weights.R
  Q <- qr.Q(qr(Gamma), complete =  TRUE)
  F_o <- Q[, (nlv + 1):ncase, drop = FALSE]
  # FIXME: Should compare the svd of matlab and R, unsure if it matters that they don't match
  # FIXME: Come back to this stack overflow thing for proper citation
  # Ahmed Fasih https://stackoverflow.com/a/41972818
  svd_out <- (D %*% t(Z) %*% F_o) |>
    {
      \(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))
    }()
  u <- svd_out$u
  v <- svd_out$v
  
  # TODO: Utilde deviates from Matlab because of the SVD
  Utilde <- v[, 1:nvar] %*% t(u)
  U <- F_o %*% Utilde
  D <- diag(diag(t(U) %*% Z))
  D[ov_type == 0, ov_type == 0] <- 0
  
  
  return(list(
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
  ))
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
#'
#' @examples
update_C_B_D <-
  function(X,
           nvar,
           Gamma,
           cindex,
           C,
           nlv,
           bindex,
           B,
           ncase,
           D,
           Z,
           ov_type) {
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
    F_o <- Q[, (nlv + 1):ncase]
    # FIXME: warning("It's unclear to me whether this SVD is safe. The Matlab code seems to either do a 'normal' svd or a economy svd depending on the dimensionality of the input matrix. Should look into this more.")
    # https://www.mathworks.com/help/matlab/ref/double.svd.html?searchHighlight=svd&s_tid=srchtitle_support_results_1_svd#d126e1597915
    svd_out2 <- (D %*% t(Z) %*% F_o) |>
      {
        \(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))
      }()
    u <- svd_out2$u
    v <- svd_out2$v
    Utilde <- v[, 1:nvar] %*% t(u)
    U <- F_o %*% Utilde
    D <- diag(diag(t(U) %*% Z))
    D[ov_type != 1, ov_type != 1] <- 0
    
    
    if (length(bindex) > 0) {
      est <- B[bindex]
    } else {
      est <- C[cindex]
    }
    # browser()
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
update_composite_LV <-
  function(ntv, tot, nlv, j, W, A, V, X, windex_j) {
    # This updates the composite
    e <- matrix(0, nrow = 1, ncol = ntv)
    e[tot] <- 1
    H1 <- diag(ntv)
    H2 <- diag(nlv)
    H1[tot, tot] <- 0
    H2[j, j] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)
    
    
    vecZDelta <- c(X %*% Delta)
    beta <- e - A[j, ]
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

#' Convert Matlab to R
#' 
#' Matlab vectors are matrices when read into R, so here we perform the conversions to make them more R-friendly. This is also a placeholder for any future needed adjustments to make the Matlab environment swappable with R
#'
#' @param mat_env The (piped-in) matlab lab environment from readMat()
#' @param vectorness Vector object that lists the different objects in the pre-swap environment that were vectors
#'
#' @return Corrected Matlab
#' @export
#'
#' @examples
convert_matlab2R <- function(mat_env, vectorness) {
  # Want both objects that were vectors and that exist in the matlab environment
  # Some objects are R-specific
  matrix2vectorize <-
    names(vectorness[vectorness == TRUE &
                       (names(vectorness) %in% names(mat_env))])
  
  for (name in matrix2vectorize) {
    # TODO: I wonder if there's a way to parallelize this somehow...
    # I couldn't get lapply to work, presumably because of some environment issue
    mat_env[[name]] <- as.vector(mat_env[[name]])
  }
  # mat_env[[as.list(matrix2vectorize)]] |>
  #   as.list() |>
  #   lapply(\(matrix2vectorize) mat_env[[matrix2vectorize]] <- as.vector(mat_env[[matrix2vectorize]]))
  return(mat_env)
}

#' Helper Function for parsing GSCA Pro V1.2.1 Full Results Files 
#' 
#' The commented headers require more work to get working
#' 
#' @return
#' @export
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom here here
#' 
#' @examples
parse_GSCAPro_FullResults <-
  function(file_path = list(
    "tests",
    "comparisons",
    "igsca_translation",
    "GSCAPro_1_2_1Output",
    "Tutorial_IGSCA_model1_full_result_ver0.csv"
  ),
  headers = c(# "Model fit measures",
    "Weights",
    "Loadings",
    "Path coefficients",
    # "Component/Factor correlations",
    # "HTMT",
    # "Construct quality measures",
    # "Fornell-Larcker criterion values",
    # "R squared values of indicators in measurement model",
    # "R squared values of component/correlations in structural model",
    # "VIFs",
    # "F squared values",
    # "Sample correlations (lower diagonal) & Residual correlations (upper diagonal)",
    # "Correlations between indicators and component/correlations",
    "Unique D^2")) {
    file <- here::here(file_path)
    
    to_parse <-
      data.table::fread(file = file, fill = TRUE, sep = ",")
    
    col_1_strings <- unlist(to_parse[, 1])
    idx_of_headers <- which((col_1_strings %in% headers) == TRUE)
    
    idx_of_blank_rows <- which((col_1_strings %in% "") == TRUE)
    
    # Gets the idx of idx_of_blank_rows of the first blank that occurs after the header
    list_of_blanks <-
      lapply(idx_of_headers, FUN = \(x) which.max(x < idx_of_blank_rows))
    idx_first_blank_after_header <-
      lapply(list_of_blanks, FUN = \(x) idx_of_blank_rows[x]) |>
      unlist()
    
    # Add one to headers because the first row with the column names is one row after the header name
    nrows_of_tables <-
      idx_first_blank_after_header - (idx_of_headers + 1)
    
    if (length(idx_of_headers) != length(nrows_of_tables)) {
      stop("The number of rows and the the number of indices of headers does not match")
    }
    
    tables <- mapply(
      FUN = data.table::fread,
      file = file,
      fill = TRUE,
      skip = idx_of_headers + 1,
      # Add one because in Version 1.2.1 of GSCA Pro the Header was one row ahead of the column labels
      nrows = nrows_of_tables,
      # Minus one for adjusting for the idx_of_headers and another for the empty-gaps between tables
      sep = ","
    )
    names(tables) <- headers
    
    # Write Tables to Directory
    mapply(
      FUN = data.table::fwrite,
      x = tables,
      file = here::here(file_path[1:(length(file_path) - 1)], paste0(names(tables), ".csv"))
    )
    
    
    
    return(tables)
    
  }

#' Converts Output of igsca functions into a table to facilitate comparisons
#' 
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#'
#' @param model 
#' @param weights 
#' @param loadings 
#' @param uniqueD 
#' @param paths 
#' @importFrom lavaan lavaanify
#' @return
#' @export
#'
#' @examples
get_lavaan_table_igsca_matrix <- function(model, weights, loadings, uniqueD, paths) {
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"),]
  # Pre-allocate Columns
  table <-
    cbind(table, list(
      "weights" = 0,
      "loadings" = 0,
      "uniqueD" = 0,
      "paths" = 0
    ))
  
  # Slide in weights
  for (indicator in rownames(weights)) {
    for (lv in colnames(weights)) {
      table[((table$lhs == lv &
              table$rhs == indicator) &
              table$op %in% c("<~", "=~")), "weights"] <- weights[indicator, lv]
    }
  }
  
  # Slide in loadings
  for (indicator in rownames(loadings)) {
    for (lv in colnames(loadings)) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "loadings"] <- loadings[indicator, lv]
    }
  }
  
  # Slide in uniqueD
  for (indicator in names(uniqueD)) {
    # This assumes that every indicator only loads onto one factor
    # Cross-factor loadings will not work with this
      table[((table$rhs == indicator) &
               (table$op == "=~")), "uniqueD"] <- uniqueD[indicator]
  }
  
  # Slide in Paths
  for (lv_from in rownames(paths)) {
    for (lv_to in colnames(paths)) {
      table[((table$rhs == lv_from &
                table$lhs == lv_to) &
               table$op == "~"), "paths"] <- paths[lv_from, lv_to]
    }
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  
  return(table)
  
}


#' Takes GSCAPro input and Creates a Lavaan-style Table
#'
#' Assumes that every indicator loads onto only one latent variable (composite/factor)
#' 
#' @param gscapro_in 
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
get_lavaan_table_igsca_gscapro <- function(gscapro_in, model) {
  
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"),]
  # Pre-allocate Columns
  table <-
    cbind(table, list(
      "weights" = 0,
      "loadings" = 0,
      "uniqueD" = 0,
      "paths" = 0
    ))
  
  
  # Correct names of latent variables
  gscapro_in$Weights$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  gscapro_in$Weights$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  
  
  lv_idx <- list()
  span <- list() 
  
  lv_idx$weights <- which(is.na(gscapro_in$Weights$Estimate))
  # Number of rows after the row of the LV that correspond to that LV's indicators
  span$weights <- diff(lv_idx$weights)
  gscapro_in$Weights$lv <- ""
  gscapro_in$Weights$lv <-
    c(with(gscapro_in$Weights, rep(V1[is.na(Estimate) &
                                      (nchar(V1) > 0)], times = span$weights)), "")
  
  
  # Retrieving indicators from gscapro
  gscapro_indicators <-
    with(gscapro_in$Weights, V1[!(V1 %in% lv) & (nchar(V1) > 0)])
  gscapro_lv <-
    with(gscapro_in$Weights, unique(lv)[nchar(unique(lv)) > 0])
  
  rownames(gscapro_in$Weights)<-gscapro_in$Weights$V1
  
  for (indicator in gscapro_indicators) {
    for (lv in gscapro_lv) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "weights"] <- gscapro_in$Weights[indicator, "Estimate"]
    }
  }
  
  
  # Loadings Parser Shortcut
  gscapro_in$Loadings$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  gscapro_in$Loadings$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  if (identical(gscapro_in$Weights$V1, gscapro_in$Loadings$V1)) {
    rownames(gscapro_in$Loadings) <- gscapro_in$Loadings$V1
    
    for (indicator in gscapro_indicators) {
      for (lv in gscapro_lv) {
        table[((table$lhs == lv &
                  table$rhs == indicator) &
                 table$op %in% c("<~", "=~")), "loadings"] <-
          gscapro_in$Loadings[indicator, "Estimate"]
      }
    }
  } else {
    stop("Cannot take the shortcut of weight indicators for loadings")
  }
  
  # Paths need their own parser 
  
  gscapro_in$`Path coefficients` <-
    gscapro_in$`Path coefficients`[!is.na(gscapro_in$`Path coefficients`$Estimate), ]
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "->",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = " ")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  # ephpostfacto on November 6, 2009; editted by Jilber Urbina on Jan 1/2014 https://stackoverflow.com/a/1690753
  paths <-strsplit(gscapro_in$`Path coefficients`$V1, " ")
  gscapro_in$`Path coefficients`$lvfrom <-sapply(paths, FUN = \(.) .[1])
  gscapro_in$`Path coefficients`$lvto <-sapply(paths, FUN = \(.) .[2])
  
  for (lv_to_iter in gscapro_in$`Path coefficients`$lvto) {
    for (lv_from_iter in gscapro_in$`Path coefficients`$lvfrom) {
      table[((table$rhs == lv_from_iter &
                table$lhs == lv_to_iter) &
               table$op == "~"), "paths"] <- with(gscapro_in$`Path coefficients`, Estimate[(lvfrom == lv_from_iter) & (lvto == lv_to_iter)])
                 
    }
  }
  
  # Parsing for Uniqueness Terms
  for (indicator in names(gscapro_in$`Unique D^2`)) {
    table[((table$rhs == indicator) &
             (table$op == "=~")), "uniqueD"] <- gscapro_in$`Unique D^2`[indicator]
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  return(table)
}
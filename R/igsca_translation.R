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
#'r
#' (igsca_sim_out <- igsca_sim(Z0 = igsca_sim_in$Z0, W0 = igsca_sim_in$W0, C0 = igsca_sim_in$C0,
#'                             B0 = igsca_sim_in$B0, lv_type = igsca_sim_in$lv_type,
#'                             ov_type = igsca_sim_in$ov_type, ind_domi = igsca_sim_in$ind_domi,
#'                             nbt = 0,
#'                             testEquivalence = TRUE)
#'                             )
igsca_sim <- function(Z0, W0, C0, B0, lv_type, ov_type, ind_domi, nbt, testEquivalence = FALSE,
                      swap_step = c("data", "qr1", "gsca_initialization", "svd1", "qr2", "svd2"),
                      itmax = 100, ceps = 0.001) {
  

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

## Initial Estimates and Preparation ---------------------------------------
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
    

## Alternating Least Squares Algorithm -------------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (any(bindex)) {
      est <- B[bindex] 
    } else {
      est <- C[cindex]
    }
    est0 <- est + 1 # It's set early because of the upcoming condition
    it <- 0
    
    
    # TODO: Review better way to do this while condition
    while ((sum(abs(est0 - est)) > ceps) && (it <= itmax)) {
      it <- it + 1
      est0 <- est
      # TODO: X deviates from Matlab because it is an offspring of svd_out
      X <- Z - U %*% D
      
      
      warning("I don't quite understand why this is the equivalent of the Matlab expression, revisit page 44 of Hiebeler 2015 R and Matlab")
      # TODO: Should find alternative to solve() to avoid matrix inversions
      WW <-
        t(C) %*% solve((C %*% t(C) + diag(nlv) - 2 * B + (B %*% t(B))))
      A <- cbind(C, B) 

### Update Each Latent Variable One-by-One------------------------------------
      for (j in seq_len(nlv)) {
        tot <- nvar + j
        windex_j <- (W0[, j] == 1)
        Xj <- X[, windex_j]

#### Composite-Specific Update -----------------------------------------------
        if (lv_type[j] == 0) {
          warning("There may be a bug in using lv_type because if it's numeric then it might only ever call the first column/row perhaps boolean would be better?")
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
          
        } else if (lv_type[j] == 1) {
#### Factor-Specific Updates -------------------------------------------------
          theta <- WW[windex_j, j]
        } else {
          stop("lv_type should either be 1 for factor or 0 for composites")
        }
        warning("Is Matlab doing a 2-norm here?")
        theta <- theta / norm(Xj %*% theta, "2")
        Gamma[, j] <- Xj %*% theta
        W[windex_j, j] <- theta
        V[windex_j, tot] <- theta
    }

### Update Loadings, Structural Coefficients and Disturbance Terms ----------

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
      warning("It's unclear to me whether this SVD is safe. The Matlab code seems to either do a 'normal' svd or a economy svd depending on the dimensionality of the input matrix. Should look into this more.")
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
    }
    

# Flip Signs for Factors and Composites Based on Dominant Indicators --------

    for (j in seq_len(nlv)) {
      if ((t(Z[, ind_domi[j]]) %*% Gamma[, j]) < 0) {
        Gamma[, j] <- (-1 * Gamma[, j])
        C[j,] <- (-1 * C[j,])
        B[B,] <- (-1 * B[j,])
        B[, j] <- (-1 * B[, j])
      }
    }
    
    mW <- W[windex] 
    # mC <- c[cindex] # Commented out in matlab
    mB <- B[bindex] 
    mD <- uniqueD
    
    crindex <- (t(C0) == 1)
    C_t <- t(C)
    mC <- C_t[crindex] 
  }
  
  return(list("mW" = mW, "mC" =  mC, "mB" =  mB, "mD" = mD))
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
  browser()
  
  indir <- list("dev", "Notes", "data", "matlab")
  # TODO: Continue filling this out so that it iterates through all the objects using lapply and automatically substitutes names
  lapply(extracted_matrices, write.csv, )
  write.csv(extracted_matrices$Z0, file = here::here(indir, "Z0.csv"))
  write.csv()
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
#' z_ini <- read.csv("~/Documents/RStudio/cSEM/dev/Notes/data/rick_mg.csv") |>
#'            subset(select=-gender)
#'            
#' w_ini <- matrix(data = c(1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          1, 0, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 1, 0, 0,
#'                          0, 0, 1, 0,
#'                          0, 0, 1, 0,
#'                          0, 0, 1, 0,
#'                          0, 0, 1, 0,
#'                          0, 0, 0, 1,
#'                          0, 0, 0, 1,
#'                          0, 0, 0, 1),
#'                 byrow = TRUE,
#'                 nrow = ncol(z_ini)
#'                 )
#'                 
#' b_ini <- matrix(data = c(0, 1, 0, 0,
#'                          0, 0, 1, 1,
#'                          0, 0, 0, 0,
#'                          0, 0, 0, 0),
#'                 byrow = TRUE,
#'                 nrow = ncol(w_ini)
#'                 )
#'                          
#' (gsca_inione_test <- gsca_inione(Z0 = z_ini, W0 = w_ini, B0 = b_ini))
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
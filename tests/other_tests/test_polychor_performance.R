devtools::load_all()


# Old version (for comparison)
calculateIndicatorCorOld <- function(
  .X_cleaned           = NULL,
  .approach_cor_robust = "none"
){
  
  is_numeric_indicator <- lapply(.X_cleaned, is.numeric)
  
  only_numeric_cols <- all(unlist(is_numeric_indicator))
  
  if(.approach_cor_robust != "none" && !only_numeric_cols) {
    stop2("Setting `.approach_cor_robust = ", .approach_cor_robust, "` requires all",
          " columns of .data to be numeric.")
  }
  
  ## polycor::hetcor() is relatively slow. If all columns are numeric use cor
  ## directly
  switch (.approach_cor_robust,
          "none" = {
            if(only_numeric_cols) {
              S <- cor(.X_cleaned)
              cor_type <- "Pearson" 
              thres_est = NULL
            } else {
              
              # Indicator's correlation matrix
              S <- matrix(0, ncol = ncol(.X_cleaned), nrow = ncol(.X_cleaned),
                             dimnames = list(colnames(.X_cleaned), colnames(.X_cleaned)))
              # matrix containing the type of correlation 
              cor_type <- S
              
              # list for the thresholds
              thres_est <- NULL
              
              # temp is used to only calculate the correlations between two 
              # indicators once (upper triangular matrix)
              temp <- colnames(.X_cleaned)
              for(i in colnames(.X_cleaned)){
                temp <- temp[temp!=i]
                for(j in temp){
                  # If both indicators are not continous, the polychoric 
                  # correlation is calculated
                  if (is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == FALSE){
                    # The polycor package gives a list with the polychoric correlation and
                    # the thresholds estimates
                    cor_temp <- polycor::polychor(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$row.cuts
                    thres_est[[j]] <- cor_temp$col.cuts
                    
                    # If one indicator is continous, the polyserial correlation 
                    # is calculated.Note: polyserial needs the continous 
                    # indicator as the first argument.
                  }else if(is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == TRUE){
                    # The polycor package gives the polyserial correlation and the thresholds
                    cor_temp <- polycor::polyserial(.X_cleaned[,j], .X_cleaned[,i], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$cuts
                    thres_est[[j]] <- NA
                  }else if(is_numeric_indicator[[i]] == TRUE & is_numeric_indicator[[j]] == FALSE){
                    cor_temp <- polycor::polyserial(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[j]] <- cor_temp$cuts
                    thres_est[[i]] <- NA
                    
                    # If both indicators are continous, the Pearson correlation
                    # is calculated.
                  }else{
                    S[i,j] <- cor(.X_cleaned[,i], .X_cleaned[,j])
                    cor_type[i,j] <- "Pearson"
                    thres_est[[i]] <- NA
                    thres_est[[j]] <- NA
                  }
                }
              }
              S <- S + t(S)
              diag(S) <- 1
              cor_type <- unique(c(cor_type))
              cor_type <- cor_type[cor_type != "0"]
              
            
              
              # The lavCor function does no smoothing in case of empty cells, which creates problems during bootstrap
              # # Use lavCor function from the lavaan package for the calculation of the polychoric and polyserial correlation 
              # # No smoothing is conducted to ensure positive definiteness of the correlation matrix
              # S <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "cor")
              # 
              # # Estimate thresholds
              # thres_est <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "th")
              # 
              # # Define type of correlation, that can either be polyserial or polychoric
              # type_var=unlist(sapply(.X_cleaned,class))
              # 
              # # if at least one numeric variable is included the polyserial correlation is applied
              # if('numeric' %in% type_var){
              #   cor_type = "Polyserial"
              # } else { #only if all variables are categorical the type of correlation is set to polychoric
              #   cor_type = "Polychoric"
              # }

            }
          },

          "mcd" = {
            S <- MASS::cov.rob(.X_cleaned, cor = TRUE, method = "mcd")$cor
            S[upper.tri(S) == TRUE] = t(S)[upper.tri(S) == TRUE]

            cor_type <-  "Robust (MCD)"
            
            thres_est = NULL
          },
          "spearman" = {
            S <- cor(.X_cleaned, method = "spearman")

            cor_type <-  "Robust (Spearman)"
            
            thres_est = NULL
          }
  )
  # (TODO) not sure how to name the "type" yet and what to do with it. Theoretically,
  # a polycoric correlation could also be used with GSCA or some other non-PLS-PM method.
  list(S = S, cor_type = cor_type, thres_est = thres_est)
}


vars <- paste0("x", seq_len(9))
S <- matrix(
  c(1.000, 0.827, 0.848, 0.454, 0.444, 0.451, 0.181, 0.179, 0.175,
    0.827, 1.000, 0.815, 0.441, 0.435, 0.438, 0.166, 0.160, 0.152,
    0.848, 0.815, 1.000, 0.442, 0.425, 0.438, 0.168, 0.164, 0.157,
    0.454, 0.441, 0.442, 1.000, 0.926, 0.932, 0.414, 0.413, 0.415,
    0.444, 0.435, 0.425, 0.926, 1.000, 0.917, 0.419, 0.416, 0.421,
    0.451, 0.438, 0.438, 0.932, 0.917, 1.000, 0.410, 0.411, 0.416,
    0.181, 0.166, 0.168, 0.414, 0.419, 0.410, 1.000, 0.833, 0.846,
    0.179, 0.160, 0.164, 0.413, 0.416, 0.411, 0.833, 1.000, 0.820,
    0.175, 0.152, 0.157, 0.415, 0.421, 0.416, 0.846, 0.820, 1.000),
  nrow = 9, ncol = 9
)


thr5 <- c(-Inf, -2.2, 0.8, 0.2, 0.7, 1.9, Inf)
thr1 <- c(-Inf, 0.3, Inf)

n <- 2000
X <- as.data.frame(mvtnorm::rmvnorm(n, mean = rep(0, 9), sigma = S))
colnames(X) <- vars

Z5 <- as.data.frame(lapply(X, FUN = cut, breaks = thr5, ordered_result = TRUE))
Z1 <- as.data.frame(lapply(X, FUN = cut, breaks = thr1, ordered_results = TRUE))

colnames(Z5) <- paste0("z5", vars)
colnames(z1) <- paste0("z1", vars)

# ------------------------------------------------------------------------------
# Compare Numerical Results -- Polychoric Correlations
# ------------------------------------------------------------------------------

# Compare results to polycor::polychor()
round(abs(cSEM:::calculateIndicatorCor(Z5)$S - calculateIndicatorCorOld(Z5)$S), 3)
#>      z5x1 z5x2 z5x3  z5x4  z5x5  z5x6 z5x7 z5x8 z5x9
#> z5x1    0    0    0 0.000 0.000 0.000    0    0    0
#> z5x2    0    0    0 0.000 0.000 0.000    0    0    0
#> z5x3    0    0    0 0.000 0.000 0.000    0    0    0
#> z5x4    0    0    0 0.000 0.044 0.052    0    0    0
#> z5x5    0    0    0 0.044 0.000 0.030    0    0    0
#> z5x6    0    0    0 0.052 0.030 0.000    0    0    0
#> z5x7    0    0    0 0.000 0.000 0.000    0    0    0
#> z5x8    0    0    0 0.000 0.000 0.000    0    0    0
#> z5x9    0    0    0 0.000 0.000 0.000    0    0    0

# Compare results to lavaan::lavCor()
round(abs(cSEM:::calculateIndicatorCor(Z5)$S - lavaan::lavCor(Z5, ordered = colnames(Z5))), 3)
#>      z5x1 z5x2 z5x3 z5x4 z5x5 z5x6 z5x7 z5x8 z5x9
#> z5x1    0                                        
#> z5x2    0    0                                   
#> z5x3    0    0    0                              
#> z5x4    0    0    0    0                         
#> z5x5    0    0    0    0    0                    
#> z5x6    0    0    0    0    0    0               
#> z5x7    0    0    0    0    0    0    0          
#> z5x8    0    0    0    0    0    0    0    0     
#> z5x9    0    0    0    0    0    0    0    0    0

round(abs(cSEM:::calculateIndicatorCor(Z1)$S - lavaan::lavCor(Z1, ordered = colnames( Z1))), 3)
#>    x1 x2 x3 x4 x5 x6 x7 x8 x9
#> x1  0                        
#> x2  0  0                     
#> x3  0  0  0                  
#> x4  0  0  0  0               
#> x5  0  0  0  0  0            
#> x6  0  0  0  0  0  0         
#> x7  0  0  0  0  0  0  0      
#> x8  0  0  0  0  0  0  0  0   
#> x9  0  0  0  0  0  0  0  0  0

# ------------------------------------------------------------------------------
# Compare Numerical Results -- Polyserial Correlations
# ------------------------------------------------------------------------------
M <- cbind(Z5, X)

# Compare results to polycor::polychor()
round(abs(calculateIndicatorCor(M)$S - calculateIndicatorCorOld(M)$S), 3)
#>       z5x1  z5x2  z5x3  z5x4  z5x5  z5x6  z5x7  z5x8  z5x9    x1    x2    x3    x4    x5    x6    x7    x8    x9
#> z5x1 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.011 0.002 0.003 0.001 0.001 0.002 0.008 0.008 0.008
#> z5x2 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.012 0.001 0.002 0.003 0.008 0.004 0.006 0.005 0.007
#> z5x3 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.013 0.002 0.001 0.001 0.006 0.003 0.002 0.004 0.003
#> z5x4 0.000 0.000 0.000 0.000 0.041 0.051 0.000 0.000 0.000 0.004 0.005 0.002 0.001 0.005 0.002 0.005 0.006 0.008
#> z5x5 0.000 0.000 0.000 0.041 0.000 0.044 0.000 0.000 0.000 0.001 0.002 0.006 0.000 0.001 0.003 0.008 0.007 0.004
#> z5x6 0.000 0.000 0.000 0.051 0.044 0.000 0.000 0.000 0.000 0.003 0.004 0.002 0.003 0.003 0.001 0.002 0.000 0.002
#> z5x7 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.007 0.003 0.002 0.000 0.003 0.002 0.001 0.001 0.004
#> z5x8 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.006 0.004 0.006 0.004 0.006 0.010 0.002 0.001 0.004
#> z5x9 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.002 0.000 0.002 0.002 0.003 0.001 0.001 0.001 0.001
#> x1   0.011 0.012 0.013 0.004 0.001 0.003 0.007 0.006 0.002 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x2   0.002 0.001 0.002 0.005 0.002 0.004 0.003 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x3   0.003 0.002 0.001 0.002 0.006 0.002 0.002 0.006 0.002 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x4   0.001 0.003 0.001 0.001 0.000 0.003 0.000 0.004 0.002 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x5   0.001 0.008 0.006 0.005 0.001 0.003 0.003 0.006 0.003 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x6   0.002 0.004 0.003 0.002 0.003 0.001 0.002 0.010 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x7   0.008 0.006 0.002 0.005 0.008 0.002 0.001 0.002 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x8   0.008 0.005 0.004 0.006 0.007 0.000 0.001 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> x9   0.008 0.007 0.003 0.008 0.004 0.002 0.004 0.004 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

# Compare results to lavaan::lavCor()
#   Here we get a few .004 differences. This is simply due to lavaan using a lower/upper bound of +/-0.995
#   whilst we use a lower/upper bound of +/- .999
round(abs(calculateIndicatorCor(M)$S - lavaan::lavCor(M, ordered = colnames(Z5))), 3)
#>       z5x1  z5x2  z5x3  z5x4  z5x5  z5x6  z5x7  z5x8  z5x9    x1    x2    x3    x4    x5    x6    x7    x8    x9
#> z5x1 0.000                                                                                                      
#> z5x2 0.000 0.000                                                                                                
#> z5x3 0.000 0.000 0.000                                                                                          
#> z5x4 0.000 0.000 0.000 0.000                                                                                    
#> z5x5 0.000 0.000 0.000 0.000 0.000                                                                              
#> z5x6 0.000 0.000 0.000 0.000 0.000 0.000                                                                        
#> z5x7 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                                                  
#> z5x8 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                                            
#> z5x9 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                                      
#> x1   0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                                
#> x2   0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                          
#> x3   0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                                    
#> x4   0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                              
#> x5   0.000 0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                        
#> x6   0.000 0.000 0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000                  
#> x7   0.000 0.000 0.000 0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000            
#> x8   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000      
#> x9   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.004 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000


# ------------------------------------------------------------------------------
# Becnmark
# ------------------------------------------------------------------------------
# polychoric correlations 5 (finite) thresholds
rbenchmark::benchmark(
  old = calculateIndicatorCorOld(Z5),
  new = calculateIndicatorCor(Z5),
  lav = lavaan::lavCor(Z5, ordered = colnames(Z5))
)
#>   test replications elapsed relative user.self sys.self user.child sys.child
#> 3  lav          100   7.114    3.462     7.080    0.008          0         0
#> 2  new          100   2.055    1.000     2.042    0.002          0         0
#> 1  old          100  79.230   38.555    78.909    0.011          0         0

# polychoric correlations 2 (finite) thresholds
rbenchmark::benchmark(
  old = calculateIndicatorCorOld(Z1),
  new = calculateIndicatorCor(Z1),
  lav = lavaan::lavCor(Z5, ordered = colnames(Z1))
)
#>   test replications elapsed relative user.self sys.self user.child sys.child
#> 3  lav          100   5.569    2.573     5.545    0.007          0         0
#> 2  new          100   2.164    1.000     2.152    0.006          0         0
#> 1  old          100   9.579    4.427     9.543    0.003          0         0

# mixed
rbenchmark::benchmark(
  old = calculateIndicatorCorOld(cbind(X, Z1, Z5)),
  new = calculateIndicatorCor(cbind(X, Z1, Z5)),
  lav = lavaan::lavCor(cbind(X, Z1, Z5), ordered = c(colnames(Z1), colnames(Z5))),
)
#>   test replications elapsed relative user.self sys.self user.child sys.child
#> 3  lav          100  38.570    1.090    38.421    0.017          0         0
#> 2  new          100  35.393    1.000    35.038    0.233          0         0
#> 1  old          100 157.684    4.455   157.164    0.042          0         0

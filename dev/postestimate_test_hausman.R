#' Hausman test
#' 
#' Calculates a bootstrap-based Hausman test that can be used to compare 
#' OLS to 2SLS estimates 
#' or 2SLS to 3SLS estimates \insertCite{Hausman1978,Wong1996}{cSEM} 
#' (TODO: Needs to be checked whether the implemention is correct, I doubt).
#' 
#' @usage testHausman(
#'  .object         = NULL,
#'  .seed           = args_default()$.seed,
#'  .alpha          = args_default()$.alpha,
#'  .R              = args_default()$.R,
#'  .R2             = args_default()$.R2,
#'  .vcv_asymptotic = args_default()$.vcv_asymptotic
#'  )
#' 
#' @param .object A 2SLS or 3SLS object 
#' 
#' @inheritParams  csem_arguments
#' 
#' @inherit csem_test return
#' 
#' @seealso [csem()], [foreman()], [cSEMResults]
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export

testHausman <- function(
  .object         = NULL,
  .seed           = args_default()$.seed,
  .alpha          = args_default()$.alpha,
  .R              = args_default()$.R,
  .R2             = args_default()$.R2,
  .vcv_asymptotic = args_default()$.vcv_asymptotic
){
  # Implementation is based on:
  # Wong (1996) - Bootstrapping Hausman's exogeneity test
  
  UseMethod("testHausman")
}

#' @describeIn testHausman (TODO)
#' @export

testHausman.cSEMResults_default <- function(
  .object         = NULL,
  .seed           = args_default()$.seed,
  .alpha          = args_default()$.alpha,
  .R              = args_default()$.R,
  .R2             = args_default()$.R2,
  .vcv_asymptotic = args_default()$.vcv_asymptotic
  ){
  
  # Check whether either 2SLS or 3SLS was used 
  if(!(.object$Information$Arguments$.approach_paths %in% c("2SLS", "3SLS"))){
    stop2("In order to conduct a Hausman test, the structural model must be either estimated by 2SLS or 3SLS.")
  }
  
  if(.object$Information$Arguments$.approach_paths == "3SLS"){
    stop2("Comparison between 3SLS and 2SLS estimates is not implemented yet.")
  }
  
  
  # If structural model was estimated by 2SLS, the estimates should be compared to OLS
  # OLS estimator is efficient and consistent under H0
  # 2SLS is only consistent under H0
  if(.object$Information$Arguments$.approach_paths == "2SLS"){
    # Estimate model with OLS
    arguments_efficient <- .object$Information$Arguments
    
    # Remove instruments and set estimator path to OLS
    arguments_efficient$.approach_paths <- 'OLS'
    
    # FOR FUTURE: Overwrite cSEM_model$.instruments!!!!! 
    arguments_efficient$.instruments <- NULL # set to NULL to avoid warning message
    
    # Reestimation by OLS
    res_efficient=do.call(csem,arguments_efficient)
    
    # Under H0 the 2SLS is consistent but not asymptotically efficient 
    
    # If it is already a resample object select the default estimation.
    res_consistent <- .object
  }
  
  # If structural model was estimated by 3SLS, these estimates should be compared to 2SLS:
  # 3SLS estimator is efficient under H0
  # 2SLS is only consistent but not efficient under H0
  if(.object$Information$Arguments$.approach_paths == "3SLS"){
    # Estimate model with 2SLS
    arguments_consistent <- .object$Information$Arguments
    
    # Set estimator path to 2SLS
    arguments_consistent$.approach_paths <- '2SLS'
    
    
    # Reestimation by 2SLS
    res_consistent=do.call(csem,arguments_consistent)
    
    # Under H0 the 3SLS is asymptotically efficient and consistent
    res_efficient <- .object
  }
  
  
  
  # Bootstrap results of efficient and consistent estimator
  # I deliaberatly ignore inadmissible solution to ensure that both habe the same number of bootstrap.
  # For the future that should be allowed
  
  if(is.null(.seed)) {
    .seed <- sample(.Random.seed, 1)
  }

  boot_efficient <- resamplecSEMResults(res_efficient,.seed = .seed,.handle_inadmissibles = 'ignore',.R = .R2)
  coef_efficient <- boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Original
  
  # Bootstrap 2SLS estimates with same seed as the OLS estimate
  boot_consistent <- resamplecSEMResults(res_consistent,.seed = .seed,.handle_inadmissibles = 'ignore',.R = .R2)
  coef_consistent <- boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Original
  
  
  # dependent variables of the equations in which instruments are used
  # One could also think about investigating all euqations however, 
  # this currently leads to problems in bootstrap because of no variation 
  m <- .object$Information$Model$structural
  dep_vars <- names(.object$Information$Model$instruments)
  
  indep_vars <- lapply(dep_vars, function(x){
    colnames(m)[which(m[x, ,drop = FALSE] == 1, arr.ind = TRUE)[, "col"]]})
  names(indep_vars) <- dep_vars
  
  endo_vars <- lapply(dep_vars, function(x){setdiff(indep_vars[[x]],colnames(.object$Information$Model$instruments[[x]]))})
  names(endo_vars) = dep_vars
  
  belongs <- lapply(dep_vars,function(x){paste(x,endo_vars, sep = ' ~ ')})
  names(belongs) <- dep_vars
  
  # # Ols approach 
  # test=strsplit(names(boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Original) , split = ' ~ ')
  # 
  # test1 = sapply(test, function(x){
  #   x[1]
  #     })
  # 
  # belongs <- lapply(dep_vars, function(x){
  #   test1 == x
  # })
  # names(belongs) <- dep_vars
  # # 
  # # indep_vars =  lapply(belongs, function(x){
  # #        sapply(test,function(x){x[2]})[x]
  # #    })
  
  # calculation of the test statistic perhaps this can also be done with the resample function
  # Needs to be discussed 
  teststat <- sapply(dep_vars, function(x){
    
    # calculate the VCV of the, under H0, asymptotically efficient and consistent estimator,
    # but inconsistent under H1
    VCV_efficient <-   cov(boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]],
                                                                                                           drop = FALSE])
    
    # calculate the VCV of the under H0 and H1 consistent estimator
    VCV_consistent <- cov(boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], 
                                                                                                           drop = FALSE])
    
    # calculate the test statistic
    para_diff <- as.matrix(coef_consistent[belongs[[x]]]-coef_efficient[belongs[[x]]])
    
    # There are two ways to calculate the VCV of the difference:
    # Either using the asymptotic VCV, i.e., VCV(beta_consistent) - VCV(beta_efficient) or
    # as VCV(beta_2SLS - beta_OLS)
    # Problem with the asymptotic VCV is that you can get negative variances
    if(.vcv_asymptotic == TRUE){
      VCV_diff <- VCV_consistent - VCV_efficient
      if(!matrixcalc::is.positive.semi.definite(VCV_diff)){
        stop2("The difference of two variance-covariance matrices is not positive semi-definite.")
      }
    }
    
    
    # If we remove inadmissible results from the bootstrap, we muss ensure that the two matrices have the same dimension
    if(.vcv_asymptotic == FALSE){
      VCV_diff <- cov(boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE]-
                        boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE])
    }
    
    # Calculation of the test statistic
    c(t(para_diff)%*%solve(VCV_diff)%*%para_diff)
    
  })# end sapply
  
  # names(teststat) <- dep_vars
  
  ## Calculate the reference distribution of the test statistic
  
  # Extract the construct scores
  scores <- .object$Estimates$Construct_scores
  
  # Calculate the predicted values and the residuals
  uandyhat <- lapply(dep_vars, function(x){
    # Calculate the predicted values of the dependent variable
    pred <- scores%*%t(m[x,,drop = FALSE])
    # Calculate the residuals
    u <- scores[,x,drop=FALSE] - pred
    list(u=u,pred=pred)
  })
  
  names(uandyhat) <- dep_vars
  
  # Needs to be done equation per equation as it is not valid to replace all 
  # the scores of the dependent variables by their predicted values at once  
  ref_dist <- lapply(dep_vars, function(dep_var){
    
    # collect the instruments for that equation
    instr <- res_consistent$Information$Model$instruments[dep_var]
    
    # Adjust the structrual model as only the equation of the considered dependent variable should be estimated
    str_model=res_efficient$Information$Model$structural[dep_var,,drop = FALSE]
    
    # Eventually replace by parseModel function
    model_star=list(structural = str_model,
                    model_type = res_consistent$Information$Model$model_type,
                    instruments = instr)
    
    # list to store the reference distribution 
    refdist <- list()
    
    # A progress bar needs to be added, see testOMF
    
    for(bb in 1:.R) {
      # draw with replacement from u
      u_star=sample(uandyhat[[dep_var]][['u']],size = nrow(uandyhat[[dep_var]][['u']]),replace = T)
      
      # u_star <- u[sample(1:nrow(u),size=nrow(u),replace=T),]
      
      # calculate new scores of dependent variable
      dep_star=uandyhat[[dep_var]][['pred']]+u_star
      
      # create new matrix where the scores of the dependent variable are replaced by the new predicted scores
      scores_star = scores
      scores_star[,dep_var] <-dep_star
      colnames(scores_star)
      
      P_star <- calculateConstructVCV(.C = cor(scores_star), #cor ensures that the predicted scores are standardized
                                             .Q = res_consistent$Estimates$Reliabilities) #From where the reliabilities
      
      
      # In the next step these scores are used to obtain the OLS and 2SLS estimates
      # in case of PLSc this is tricky as we need the reliabilities
      # As an outcome, we obtain the OLS and the 2SLS estimate of the considered equation
      if(.object$Information$Arguments$.approach_paths == "2SLS"){
        
        efficient_star <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                              .csem_model = model_star,
                                              .approach_paths = 'OLS',
                                              .H = scores_star,
                                              .normality = res_efficient$Information$Arguments$.normality,
                                              .P = P_star,
                                              .Q = res_efficient$Estimates$Reliabilities)
        
        consistent_star <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                               .csem_model = model_star,
                                               .approach_paths = '2SLS',
                                               .H = scores_star,
                                               .normality = res_efficient$Information$Arguments$.normality,
                                               .P = P_star,
                                               .Q = res_efficient$Estimates$Reliabilities,
                                               .instruments = instr)
      }
      
      
      if(.object$Information$Arguments$.approach_paths == "3SLS"){
        consistent_star <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                               .csem_model = model_star,
                                               .approach_paths = '2SLS',
                                               .H = scores_star,
                                               .normality = res_efficient$Information$Arguments$.normality,
                                               .P = P_star,
                                               .Q = res_efficient$Estimates$Reliabilities,
                                               .instruments = instr)
        
        efficient_star <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                              .csem_model = model_star,
                                              .approach_paths = '3SLS',
                                              .H = scores_star,
                                              .normality = res_efficient$Information$Arguments$.normality,
                                              .P = P_star,
                                              .Q = res_efficient$Estimates$Reliabilities,
                                              .instruments = instr)
      }
      
      
      
      # calculate the difference
      diff_star <- efficient_star$Path_estimates[,endo_vars[[dep_var]],drop = FALSE] -
        consistent_star$Path_estimates[,endo_vars[[dep_var]],drop = FALSE]
      
      
      
      
      # bootstrap sample to obtain the variance of the diff_star
      boot_star=lapply(1:.R2, function(x){
        scores_temp=dplyr::sample_n(as.data.frame(scores_star),size=nrow(scores_star),replace=T)
        # daten=as.data.frame(scale(temp))
        
        P_temp <- calculateConstructVCV(.C = cor(scores_temp), #cor ensures that the predicted scores are standardized
                                               .Q = res_consistent$Estimates$Reliabilities) # From where do the reliabilities come
        
        
        # calculate the difference
        # What about the reliabilities in case of PLSc?
        
        if(.object$Information$Arguments$.approach_paths == "2SLS"){
          efficient_temp <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                                .csem_model = model_star,
                                                .approach_paths = 'OLS',
                                                .H = scores_temp,
                                                .normality = res_efficient$Information$Arguments$.normality,
                                                .P = P_temp,
                                                .Q = res_efficient$Estimates$Reliabilities)
          
          consistent_temp <- estimatePath(.approach_nl = res_consistent$Information$Arguments$.approach_nl,
                                                 .csem_model = model_star,
                                                 .approach_paths = '2SLS',
                                                 .H = scores_temp,
                                                 .normality = res_consistent$Information$Arguments$.normality,
                                                 .P = P_temp,
                                                 .Q = res_consistent$Estimates$Reliabilities,
                                                 .instruments = instr)
        }
        
        
        if(.object$Information$Arguments$.approach_paths == "3SLS"){
          consistent_temp <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                                 .csem_model = model_star,
                                                 .approach_paths = '2SLS',
                                                 .H = scores_temp,
                                                 .normality = res_efficient$Information$Arguments$.normality,
                                                 .P = P_temp,
                                                 .Q = res_efficient$Estimates$Reliabilities,
                                                 .instruments = instr)
          
          efficient_temp <- estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                                .csem_model = model_star,
                                                .approach_paths = '3SLS',
                                                .H = scores_temp,
                                                .normality = res_efficient$Information$Arguments$.normality,
                                                .P = P_temp,
                                                .Q = res_efficient$Estimates$Reliabilities,
                                                .instruments = instr)
        }
        
        diff_temp <- consistent_temp$Path_estimates[,endo_vars[[dep_var]], drop =FALSE] - efficient_temp$Path_estimates[,endo_vars[[dep_var]], drop = FALSE]
        
        out <- list(diff = diff_temp, efficient = efficient_temp$Path_estimates[,endo_vars[[dep_var]]],
                    consistent = consistent_temp$Path_estimates[,endo_vars[[dep_var]], drop =FALSE])
        
        return(out)
      })
      
      
      boot_star <- purrr::transpose(boot_star)
      
      # calculate the VCV 
      
      if(.vcv_asymptotic == FALSE){
        diff_boot <- do.call(rbind,boot_star$diff)
        vcv_star = cov(diff_boot)
      }
      
      if(.vcv_asymptotic == TRUE){
        VCV_efficient_star <- cov(do.call(rbind,boot_star$efficient))
        VCV_consistent_star <- cov(do.call(rbind,boot_star$consistent))
        
        vcv_star = VCV_consistent_star - VCV_efficient_star 
      }
      
      # calculate the test statistic
      refdist[[bb]]=c(diff_star%*%solve(vcv_star)%*%t(diff_star))
    }#end for loop
    
    do.call(c,refdist)
    
  }) #end lapply
  
  
  names(ref_dist) <- dep_vars
  
  ref_dist_matrix <- do.call(rbind, ref_dist) 
  
  # Order the significance levels
  .alpha <- .alpha[order(.alpha)]
  
  critical_values <- matrixStats::rowQuantiles(ref_dist_matrix, 
                                               probs =  1-.alpha, drop = FALSE)
  
  ## Compare critical value and teststatistic
  decision <- teststat < critical_values
  # TRUE = no evidence against the H0 --> not reject
  # FALSE --> reject
  
  # Return output
  out <- list(
    "Test_statistic"     = teststat,
    "Critical_value"     = critical_values,
    "Decision"           = decision,
    "Information"        = list(
      # "Number_admissibles" = ncol(ref_dist_matrix),
      # "Total_runs"         = counter + n_inadmissibles,
      "Total_runs"         = ncol(ref_dist_matrix),
      "Bootstrap_values"   = ref_dist,
      "Consistent_estimator"        = "2SLS",
      "Efficient_estimator"  =  if(.object$Information$Arguments$.approach_paths == "2SLS") {"OLS"
      } else if(.object$Information$Arguments$.approach_paths == "3SLS"){"3SLS"}
    )
  )
  
  class(out) <- "cSEMTestHausman"
  return(out)
}


#' @describeIn testHausman (TODO)
#' @export

testHausman.cSEMResults_multi <- function(
  .object         = NULL,
  .seed           = args_default()$.seed,
  .alpha          = args_default()$.alpha,
  .R              = args_default()$.R,
  .R2             = args_default()$.R2,
  .vcv_asymptotic = args_default()$.vcv_asymptotic
  ){
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, testHausman.cSEMResults_2ndorder,
           .seed           = .seed,
           .alpha          = .alpha,
           .R              = .R,
           .R2             = .R2,
           .vcv_asymptotic = .vcv_asymptotic
    )
  } else {
    lapply(.object, testHausman.cSEMResults_default,
           .seed           = .seed,
           .alpha          = .alpha,
           .R              = .R,
           .R2             = .R2,
           .vcv_asymptotic = .vcv_asymptotic
           )
  }
}

#' @describeIn testHausman (TODO)
#' @export

testHausman.cSEMResults_2ndorder <- function(
  .object         = NULL,
  .seed           = args_default()$.seed,
  .alpha          = args_default()$.alpha,
  .R              = args_default()$.R,
  .R2             = args_default()$.R2,
  .vcv_asymptotic = args_default()$.vcv_asymptotic
  ){
  
  stop2("Hausman test is not yet implemented for second-order models.")
}

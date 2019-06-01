### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.

# Empirical Bayes correction for loadings suggested by Dijkstra (2018) can be 
# extended to other parameter as well .
# I recommend to extend it to construct correlations instead of the path 
# coefficients and restimate them as we now the bounds of such correlations
# However, this requires that we bootstrap the SEs for the construct correlations.
quasiEmpiricalBayesCorrection <- function(.object,.method=c('median','mean')){
  # Loadings=.object$Estimates$Loading_estimates
  # Loadings[1,2]=2
  # Indmatrix=which(abs(Loadings)>1,arr.ind = F)

  
  L=.object$Estimates$Loading_estimates[.object$Information$Model$measurement!=0]
  

    
    sig_hat=apply(.object$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled,2,sd)
    
    # adjustedLoading=c()
    if(.method=='median'){
    for(i in which(abs(L)>1)){
      A=pnorm((-1-L[i])/sig_hat[i])
      B=pnorm((1-L[i])/sig_hat[i])
      # Overwrite the old loadings
      L[i]=L[i] + sig_hat[i] * qnorm(mean(c(A,B)),0,1)
      }
    }

    if(.method=='mean'){
      for(i in which(abs(L)>1)){
        A=pnorm((-1-L[i])/sig_hat[i])
        B=pnorm((1-L[i])/sig_hat[i])
        # Overwrite the old loadings
        L[i]=L[i] + sig_hat[i]/(B-A) *( dnorm((-1-L[i])/sig_hat[i],0,1) - dnorm((1-L[i])/sig_hat[i],0,1))
      }
     }
    
    # Overwrite the old loadings
    .object$Estimates$Loading_estimates[.object$Information$Model$measurement!=0]=L
    
    .object
  }
 
predict=function(.object, testDataset){
  
  # Check whether the test dataset contains the same inidcators as the train dataset
  if(sum(!(colnames(.object$Information$Model$measurement)%in% colnames(testDataset)))!=0){
    stop('The same indicators as in the original estimation need to be provided') #Perhaps we do not need to be tha strict as we only need the exogenous indicators
  }
  
  # Return error if model is non-linear. See danks et al how this can be addressed.
  if(.object$Information$Model$model_type != 'Linear'){
    stop2('Currenlty, predictPLS works only for linear models.')
  }
  
  if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
    stop2('Currently, predictPLS works only in combination with Pearson correlation.')
  }
  
  # Perform check of the provided dataset
  # Needs to be implemented
  
  # Order the provided dataset
  testData <- testDataset[,colnames(.object$Information$Model$measurement)]
  
  
  # save descriptives of the original unscaled train dataset
  trainData <- .object$Information$Arguments$.data
  mean_train <- colMeans(trainData)
  sd_train <- apply(trainData,2,sd)
  
  # Scale the test dataset with the descriptives of taining dataset
  testDatascale <- sapply(colnames(testData),function(x){
    (testData[,x]-mean_train[x])/sd_train[x]
  })
  
  W_train <- .object$Estimates$Weight_estimates
  Loadings_train <- .object$Estimates$Loading_estimates
  
  # Path coefficients estimates based on the trainig dataset
  pathtrain <- .object$Estimates$Path_estimates

  # Identifiy exogenous construct in the structural model
  Cons_exo <- rownames(pathtrain)[rowSums(pathtrain)==0]
  Cons_endo <- rownames(pathtrain)[rowSums(pathtrain)!=0]
  
  # Path coefficients of exogenous and endogenous constructs
  B_train      <- .object$Estimates$Path_estimates[Cons_endo, Cons_endo, drop = FALSE]
  Gamma_train  <- .object$Estimates$Path_estimates[Cons_endo, Cons_exo, drop = FALSE]
  
  # Predict scores for the exogenous constructs
  exogscores <- testDatascale%*%t(W_train[Cons_exo,,drop = FALSE])
  
  
  # calculate predictions of the endogenous constructs
  endoscores <- exogscores%*%t(Gamma_train) %*% solve(diag(nrow(B_train)) - t(B_train))
  
  
  xhat <- endoscores %*% Loadings_train[Cons_endo,,drop = FALSE]
  
  # Denormalize predictions
  xhatrescale= sapply(colnames(xhat),function(x){
    xhat[,x]*sd_train[x]+mean_train[x]
  } )
  
  # Return the predicted values of the indicators connected to endogenous constructs
  Ind_endo=colnames(.object$Information$Model$measurement)[colSums(.object$Information$Model$measurement[Cons_endo,])!=0]

  # Calculate the difference between original and predicted values
  residuals = testData[,Ind_endo] - xhatrescale[,Ind_endo]
  
  list(Pred_Indicators_endo = xhatrescale[,Ind_endo],
       Residuals_Indicators_endo = residuals, 
       Pred_Construct_endo = endoscores)
  
  
  }


 
#' Output of floodlight analysis
#'
#' Calculate the the effect of an independent variable depending on different values
#' of its moderator variable to perform a floodlight analysis \insertCite{Spiller2013}{cSEM}.
#' Moreover, the Johnson-Neyman point(s) is calculated. In doing so, it is considered 
#' whether there is a sign switch in the lower and upper confidence intervals. 
#' If so, the position of the switch is returned.
#' 
#'
#' @usage effect_moderator_two_way=function(.object=args_default()$.object,
#'                                          .steps = seq(-1.5,1.5,0.01),
#'                                          .alpha = args_default()$.alpha, 
#'                                          .dependent = NULL,
#'                                          .independent = NULL,
#'                                          .moderator = NULL )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
effect_moderator_two_way=function(.object=args_default()$.object,
                                  .steps = args_default()$.steps,
                                  .alpha = args_default()$.alpha,
                                  .dependent = args_default()$.dependent, 
                                  .independent = args_default()$.independent,
                                  .moderator = args_default()$.moderator ){
  
  # Check whether .object is of class cSEMResults_resampled
  if(!("cSEMResults_resampled" %in% class(.object))){
    stop2('.object needs to be of class cSEMResults_resampled')
  }
  
  # Check whether it is a non-linear model
  if("Nonlinear" != .object$Information$Model$model_type){
    stop2("Model type must be nonlinear.")
  }
  
  # Works only for one significance level, i.e., no vector of significances is allowed
  if(length(.alpha)!=1){
    stop2('Only one significant level is allowed, not a vector')
  }
  
  # Check whether dependent, independent, and moderator variable are provided
  if(is.null(.dependent) | is.null(.independent) | is.null(.moderator)){
    stop2("Not all variables (dependent, independent, and/or moderator) are supplied.")
  }
  
  structural = .object$Information$Model$structural
  
  # Character string containing the names of the dependent and independent variables
  dependent_vars=rownames(structural[rowSums(structural)!=0,,drop=FALSE])
  independent_vars = colnames(structural[,colSums(structural[.dependent,,drop = FALSE])!=0,drop = FALSE])
  
  # Check whether dependent variable is valid
  if(!(.dependent %in% dependent_vars)){
    stop2("Defined dependent variables is not a dependent variable in the specified model.")
  }
  
  
  # Check whether supplied variables are also used in the model
  if(!all(c(.independent, .moderator) %in% colnames(structural))){
    stop2("Independent and moderator variable do not occur together in the equation of the dependent variable.")
  }
  
  # Add stop if the variables are included in a higher-order moderation, e.g., cubic term 
  
  # Possible names of the interaction term. One could adjust the arguement that sth lik x.z needs to be provided to .moderator
  possible_names=c(paste(.independent,.moderator, sep= '.'),paste(.moderator,.independent, sep= '.'))
  
  # Name of the interaction term 
  name_interaction = possible_names[possible_names %in% independent_vars]
  
  if(length(name_interaction) != 1){
    stop2("The defined interaction term does not exist in the model or is not part of the equation of the dependent variable.")
  }
  
  
  # Effect names
  name_mod_effect = paste(.dependent, name_interaction, sep=' ~ ')
  name_single_effect = paste(.dependent, .independent, sep=' ~ ')
  

  dataplot_temp = lapply(.steps, function(x){
    
    # bootstrap effect of x on y at level step, i.e., z=.steps
    effect_boot=.object$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,name_single_effect]+
      .object$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,name_mod_effect]*x 
    
    # Value of the originally estimated effect at level .steps
    effect_at_steps=.object$Estimates$Estimates_resample$Estimates1$Path_estimates$Original[name_single_effect]+
      .object$Estimates$Estimates_resample$Estimates1$Path_estimates$Original[name_mod_effect]*x
    
    bounds = quantile(effect_boot,c(.alpha/2,1-.alpha/2))

    # Prepare and return output
    out = c(bounds[1],
      bounds[2],
      de = effect_at_steps,
      value = x)
    names(out) = c('lb', 'ub', 'direct_effect', 'steps')
    
    return(out)
  })
 
  out = do.call(rbind,dataplot_temp)
  
# Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  pos_ub = which(diff(sign(out[,'ub']))!=0)
  pos_lb = which(diff(sign(out[,'lb']))!=0) 
# Prepare and return output
    res=list(out = out, 
             Johnson_Neyman_points = list(JNlb = c(x = out[,'steps'][pos_lb], y = out[,'lb'][pos_lb]),
                                          JNub = c(x = out[,'steps'][pos_ub], y = out[,'ub'][pos_ub])),
             Information = c(alpha=.alpha, independent = .independent, moderator= .moderator, dependent = .dependent))
    class(res) = c("Two_Way_Effect", class(res))
    return(res)
} 



#' Plot the output of floodlight analysis
#'
#' It plots the direct effect of an independent variable depending on the levels of the 
#' moderator variable. Moreover, the confidence interval are plotted including the 
#' Johnson-Neyman points. 
#'
#' @usage plot(.TWobject)
#'
#'
#' @export
plot.Two_Way_Effect = function(.TWobject){
  
  require(ggplot2)
  plot1=ggplot(as.data.frame(.TWobject$out),aes(x=.TWobject$out[,'steps'],y=.TWobject$out[,'direct_effect']))+
    geom_line()+
    geom_ribbon(aes(ymin=.TWobject$out[,'lb'], ymax=.TWobject$out[,'ub']),alpha=0.2)+
    labs(x=paste('Level of ',.TWobject$Information['moderator']) , 
    y=paste('Effect of', .TWobject$Information['independent'], 'on \n', .TWobject$Information['dependent']))+
    theme_bw()+
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
    theme(panel.grid.minor = element_blank())
  
  # Plot Johnson-Neyman points, if they exist in the considered range
  if(length(.TWobject$Johnson_Neyman_points$JNlb)==2){
    JN=.TWobject$Johnson_Neyman_points$JNlb
    plot1 = plot1+
      geom_point(x=JN['x'],y=JN['y'],size=2)  
  }
  
  if(length(.TWobject$Johnson_Neyman_points$JNub)==2){
    JN=.TWobject$Johnson_Neyman_points$JNub
    plot1 = plot1+
      geom_point(x=JN['x'],y=JN['y'],size=2)  
  }
  
  # Plot
  plot1
}
  


#' Hausman test
#' 
#' Calculates a bootstrap-based Hausman test that can be used to compare OLS to 2SLS estimates \insertCite{Wong1996}{cSEM}
#' or 2SLS to 3SLS estimates (Needs to be checked whether this can be done, I think so).
#' 
#' @usage testHausman=function(.object,
#'  .seed=1234,
#'  .alpha = args_default()$.alpha,
#'  .R = args_default()$.R,
#'  .R2 = args_default()$.R2,
#'  .vcv_asymptotic = args_default()$.vcv_asymptotic)
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

testHausman=function(.object,
                     .seed=1234,
                     .alpha = args_default()$.alpha,
                     .R = args_default()$.R,
                     .R2 = args_default()$.R2,
                     .vcv_asymptotic = args_default()$.vcv_asymptotic) {
  # Implementation is based on:
  # Wong (1996) - Bootstrapping Hausman's exogeneity test
  
  UseMethod("testHausman")
}

#' @describeIn testHausman (TODO)
#' @export

testHausman.cSEMResults_default=function(.object,
                     .seed=1234,
                     .alpha = args_default()$.alpha,
                     .R = args_default()$.R,
                     .R2 = args_default()$.R2,
                     .vcv_asymptotic = args_default()$.vcv_asymptotic){

  # Check whether either 2SLS or 3SLS was used 
  if(!(.object$Information$Arguments$.approach_paths %in% c("2SLS", "3SLS"))){
    stop2("In order to conduct a Hausman test, the structural model must be either estimated by 2SLS or 3SLS.")
  }
  
  
  # If structural model was estimated by 2SLS, the estimates should be compared to OLS
  if(.object$Information$Arguments$.approach_paths == "2SLS"){
  # Estimate model with OLS
  arguments_efficient <- .object$Information$Arguments
  
  # Remove instruments and set estimator path to OLS
  arguments_efficient$.approach_paths <- 'OLS'
  arguments_efficient$.instruments <- NULL #Why do I have to overwrite the instruments? Shouldn't it be enough to set the estimator to OLS
  }
  
  # If structural model was estimated by 3SLS, the estimates should be compared to 2SLS
  if(.object$Information$Arguments$.approach_paths == "3SLS"){
    # Estimate model with 2SLS
    arguments_efficient <- .object$Information$Arguments
    
    # Set estimator path to 2SLS
    arguments_efficient$.approach_paths <- '2SLS'
    }
  
  # Reestimation by OLS or 2SLS
  res_efficient=do.call(csem,arguments_efficient)
  
  # Bootstrap OLS/2SLS estimates
  # I deliaberatly ignore inadmissible solution to ensure that both habe the same number of bootstrap.
  # For the future that should be allowed
  boot_efficient <- resamplecSEMResults(res_efficient,.seed = .seed,.handle_inadmissibles = 'ignore',.R = .R2)
  
  coef_efficient <- boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Original
  
  # Bootstrap 2SLS estimates with same seed as the OLS estimate
  boot_consistent <- resamplecSEMResults(.object,.seed = .seed,.handle_inadmissibles = 'ignore',.R = .R2)
  
  coef_consistent <- boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Original
  
  
  # dependent variables of the equations in which instruments have been used
  # One could also think about investigating all euqations however, 
  # this currently leads to problems in bootstrap because of no variation 
  m <- .object$Information$Model$structural
  dep_vars <- names(.object$Information$Model$instruments)
  
  ## Consider all equations
  # dep_vars <- rownames(m)[rowSums(m)!=0]
  
  test=strsplit(names(boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Original) , split = ' ~ ')
  
  test1 = sapply(test, function(x){
    x[1]
      })
  
  belongs <- lapply(dep_vars, function(x){
    test1 == x
  })
  names(belongs) <- dep_vars
  
  indep_var =  lapply(belongs, function(x){
         sapply(test,function(x){x[2]})[x]
     })
  
  # calculation of the test statistic
  teststat <- sapply(dep_vars, function(x){
  
    # calculate the VCV of the OLS estimates
    VCV_efficient <-   cov(boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE])
  
    # calculate the VCV of the 2SLS estimates
    VCV_consistent <- cov(boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE])
  
    # calculate the test statistic
    para_diff <- as.matrix(coef_consistent[belongs[[x]]]-coef_efficient[belongs[[x]]])
    
    # There are two ways to calculate the VCV of the difference:
    # Either using the asymptotic VCV, i.e., VCV(beta_OLS) - VCV(beta_2SLS) or
    # as VCV(beta_2SLS - beta_OLS)
    # Problem with the asymptotic VCV is that you can get negative variances
    if(.vcv_asymptotic == TRUE){
      VCV_diff <- VCV_consistent - VCV_efficient
    }
    
    
    # If we remove inadmissible results from the bootstrap, we muss ensure that the two matrices have the same dimension
    if(.vcv_asymptotic == FALSE){
      VCV_diff <- cov(boot_consistent$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE]-
                      boot_efficient$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,belongs[[x]], drop = FALSE])
    }
    
    # Calculation of the test statistic
    nrow(.object$Information$Data)*t(para_diff)%*%VCV_diff%*%para_diff
    
  })
  
  
  ## Calculate the reference distribution of the test statistic
  
  # Extract the construct scores
  scores <- .object$Estimates$Construct_scores
  
  uandyhat <- lapply(dep_vars, function(x){
    # Calculate the predicted values of the dependent variable
    pred <- scores%*%t(m[x,,drop = FALSE])
    # Calculate the residuals
    u <- scores[,x,drop=FALSE] - pred
    list(u=u,pred=pred)
  })
  
  names(uandyhat) <- dep_vars
  
  # ref_dist <- list()
  
  # uandyhattrans = purrr::transpose(uandyhat)
  # 
  # u <- do.call(cbind,uandyhattrans$u)
  # yhat <- do.call(cbind,uandyhattrans$pred)
  # 
  # Needs to be done equation per equation as it is not valid to replace all 
  # the scores of the dependent variables by their predicted values at once  
  ref_dist <- lapply(dep_vars, function(dep_var){
    
    # collect the instruments for that equation
    instr <- res_efficient$Information$Model$instruments[dep_var]
    
    # Adjust the structrual model as only the equation of the considered dependent variable should be estimated
    str_model=res_efficient$Information$Model$structural[dep_var,,drop = FALSE]
    model_star=list(structural = str_model,
                   model_type = res_efficient$Information$Model$model_type,
                   instruments = instr)
  
    # indep_vars = colnames(str_model)[str_model!=0]
    
    # list to store the reference distribution 
    refdist <- list()
    
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
      
      P_star <- cSEM:::calculateConstructVCV(.C = cor(scores_star), #cor ensures that the predicted scores are standardized
                                            .Q = res_efficient$Estimates$Reliabilities,
                                            .csem_model = res_efficient$Information$Model)
      
      # In the next step these scores are used to obtain the OLS and 2SLS estimates
      # in case of PLSc this is tricky as we need the reliabilities
      # As an outcome, we obtain the OLS and the 2SLS estimate of the considered equation
      if(.object$Information$Arguments$.approach_paths == "2SLS"){
      efficient_star <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                     .csem_model = model_star,
                                     .approach_paths = 'OLS',
                                     .H = scores_star,
                                     .normality = res_efficient$Information$Arguments$.normality,
                                     .P = P_star,
                                     .Q = res_efficient$Estimates$Reliabilities)
      
      consistent_star <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                     .csem_model = model_star,
                                     .approach_paths = '2SLS',
                                     .H = scores_star,
                                     .normality = res_efficient$Information$Arguments$.normality,
                                     .P = P_star,
                                     .Q = res_efficient$Estimates$Reliabilities,
                                     .instruments = instr)
      }
      
      
      if(.object$Information$Arguments$.approach_paths == "3SLS"){
        efficient_star <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                              .csem_model = model_star,
                                              .approach_paths = '2SLS',
                                              .H = scores_star,
                                              .normality = res_efficient$Information$Arguments$.normality,
                                              .P = P_star,
                                              .Q = res_efficient$Estimates$Reliabilities,
                                              .instruments = instr)
        
        consistent_star <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                               .csem_model = model_star,
                                               .approach_paths = '3SLS',
                                               .H = scores_star,
                                               .normality = res_efficient$Information$Arguments$.normality,
                                               .P = P_star,
                                               .Q = res_efficient$Estimates$Reliabilities,
                                               .instruments = instr)
        }
      
    
      
      # calculate the difference
      diff_star <- efficient_star$Path_estimates[,indep_var[[dep_var]],drop = FALSE] -
        consistent_star$Path_estimates[,indep_var[[dep_var]],drop = FALSE]
     
      
      
       
      # bootstrap sample to obtain the variance of the diff_star
      boot_star=lapply(1:.R2, function(x){
        scores_temp=dplyr::sample_n(as.data.frame(scores_star),size=nrow(scores_star),replace=T)
        # daten=as.data.frame(scale(temp))
        
        P_temp <- cSEM:::calculateConstructVCV(.C = cor(scores_temp), #cor ensures that the predicted scores are standardized
                                              .Q = res_efficient$Estimates$Reliabilities,
                                              .csem_model = res_efficient$Information$Model)
        
        # calculate the difference
        
        if(.object$Information$Arguments$.approach_paths == "2SLS"){
        efficient_temp <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                       .csem_model = model_star,
                                       .approach_paths = 'OLS',
                                       .H = scores_temp,
                                       .normality = res_efficient$Information$Arguments$.normality,
                                       .P = P_temp,
                                       .Q = res_efficient$Estimates$Reliabilities)
        
        consistent_temp <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                        .csem_model = model_star,
                                        .approach_paths = '2SLS',
                                        .H = scores_temp,
                                        .normality = res_efficient$Information$Arguments$.normality,
                                        .P = P_temp,
                                        .Q = res_efficient$Estimates$Reliabilities,
                                        .instruments = instr)
        }
        
        
        if(.object$Information$Arguments$.approach_paths == "3SLS"){
          efficient_temp <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                         .csem_model = model_star,
                                         .approach_paths = '2SLS',
                                         .H = scores_temp,
                                         .normality = res_efficient$Information$Arguments$.normality,
                                         .P = P_temp,
                                         .Q = res_efficient$Estimates$Reliabilities,
                                         .instruments = instr)
          
          consistent_temp <- cSEM:::estimatePath(.approach_nl = res_efficient$Information$Arguments$.approach_nl,
                                          .csem_model = model_star,
                                          .approach_paths = '3SLS',
                                          .H = scores_temp,
                                          .normality = res_efficient$Information$Arguments$.normality,
                                          .P = P_temp,
                                          .Q = res_efficient$Estimates$Reliabilities,
                                          .instruments = instr)
        }
        
        diff_temp <- efficient_temp$Path_estimates[,indep_var[[dep_var]]] - consistent_temp$Path_estimates[,indep_var[[dep_var]]]
        
        out <- list(diff = diff_temp, OLS = efficient_temp$Path_estimates[,indep_var[[dep_var]]],
                    TSLS = consistent_temp$Path_estimates[,indep_var[[dep_var]]])
        
        return(out)
        })
      
      
      boot_star <- purrr::transpose(boot_star)
      
      # calculate the VCV 
      
      if(.vcv_asymptotic == FALSE){
        diff_boot <- do.call(rbind,boot_star$diff)
        vcv_star = cov(diff_boot)
      }
      
      if(.vcv_asymptotic == TRUE){
        VCV_efficient_star <- cov(do.call(rbind,boot_star$OLS))
        VCV_consistent_star <- cov(do.call(rbind,boot_star$TSLS))

        vcv_star = VCV_consistent_star - VCV_efficient_star 
      }
      
      # calculate the test statistic
      refdist[[bb]]=nrow(scores_star)*diff_star%*%solve(vcv_star)%*%t(diff_star)
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
      "Number_admissibles" = ncol(ref_dist_matrix),
      # "Total_runs"         = counter + n_inadmissibles,
      "Bootstrap_values"   = ref_dist,
      "Consistent_estimator"        = .object$Information$Arguments$.approach_paths
    )
  )
  
  class(out) <- "cSEMTestHausman"
  return(out)
}


#' @describeIn testHausman (TODO)
#' @export

testHausman.cSEMResults_multi <- function(.object,
                                          .seed=1234,
                                          .alpha = args_default()$.alpha,
                                          .R = args_default()$.R,
                                          .R2 = args_default()$.R2,
                                          .vcv_asymptotic = args_default()$.vcv_asymptotic){
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, testHausman.cSEMResults_2ndorder,
           .seed = 1234,
           .alpha = args_default()$.alpha,
           .R = args_default()$.R,
           .R2 = args_default()$.R2,
           .vcv_asymptotic = args_default()$.vcv_asymptotic)
  } else {
    lapply(.object, testHausman.cSEMResults_default,
           .seed = 1234,
           .alpha = args_default()$.alpha,
           .R = args_default()$.R,
           .R2 = args_default()$.R2,
           .vcv_asymptotic = args_default()$.vcv_asymptotic)
    }
}

#' @describeIn testHausman (TODO)
#' @export

 testHausman.cSEMResults_2ndorder <- function(.object,
       .seed = 1234,
       .alpha = args_default()$.alpha,
       .R = args_default()$.R,
       .R2 = args_default()$.R2,
       .vcv_asymptotic = args_default()$.vcv_asymptotic){
   
   stop2("Hausman test is not yet implemented for second-order models.")
 }

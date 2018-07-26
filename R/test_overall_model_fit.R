# Test for overall model fit

testOverallModelFit <- function(
  .object,
  .dropInadmissibles=TRUE,
  .alpha=0.05,
  .runs=499,
  ...
  ) {
  
  # Extract required infromation 
  X=.object$Informatio$Dataset
  nObs=.object$Information$Number_of_observations
  S=.object$Estimates$Indicator_VCV
  Sigma_hat=fitted(.object)
  
  
  # Calculate test statistic
  teststat=c(dG=dG(.object), SRMR=SRMR(.object), dL=dL(.object))
  
  S_half=solve(expm::sqrtm(S))
  Sig_half=expm::sqrtm(Sigma_hat)
  
  # Transform dataset, see Beran & Srivastava (1985)
  
  
  X_trans = X%*%S_half%*%Sig_half
  colnames(X_trans)=colnames(X)
  
  # Calculate reference distribution
  ref_dist=lapply(replicate(.runs, X_trans[sample(1:nObs,replace=TRUE),], simplify = FALSE),function(x){
    Est_temp=csem(.data =x,
                  .model = .object$Information$Model,
                  .approach_cf =  .object$Information$Approach_cf,
                  .approach_paths =.object$Information$Path_approach,
                  .approach_weights= .object$Information$Weight_approach,
                  .disattenuate = .object$Information$Disattenuate,
                  .dominant_indicators =.object$Information$Dominant_indicators, 
                  .estimate_structural =.object$Information$Estimate_structural,
                  .ignore_structural_model =.object$Information$Ignore_structural_model,
                  .iter_max = .object$Information$Number_max_iterations,
                  .PLS_mode =.object$Information$PLS_modes,
                  .PLS_weight_scheme_inner = .object$Information$PLS_inner_weighting_scheme,
                  .tolerance =  .object$Information$Tolerance
    )
    
    status_code=status(Est_temp)
    
    # if it is controlled for inadmissible
    if(.dropInadmissibles){
      if(0 %in% status_code){
        c(dG = dG(Est_temp), SRMR = SRMR(Est_temp), dL=dL(Est_temp)) 
      }
      # else, i.e., dropInadmissible == FALSE
    }else{ 
      c(dG = dG(Est_temp), SRMR = SRMR(Est_temp), dL=dL(Est_temp))
    }
  })
  ref_dist_matrix=do.call(cbind,ref_dist)
  critical_value=apply(ref_dist_matrix,1,quantile,1-.alpha)
  
  list(Test_statistic=teststat,Critical_value=critical_value, Number_admissibles=nrow(ref_dist_matrix)) 
  
}
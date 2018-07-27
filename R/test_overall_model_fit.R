#' @export 

# Test for overall model fit

Test_for_overall_model_fit=function(.object=args_default()$.model,
                                    .dropInadmissibles=args_default()$.model,
                                    .alpha=args_default()$.model,
                                    .runs=args_default()$.model,
                                    ...){
  
  # Extract required infromation 
  X=.object$Information$Data
  nObs=nrow(X)
  S=.object$Estimates$Indicator_VCV
  Sigma_hat=fitted(.object)
  
  
  # Calculate test statistic
  teststat=c(dG=dG(.matrix1=S,.matrix2=Sigma_hat), 
             SRMR=SRMR(.object), 
             dL=dL(.matrix1=S,.matrix2=Sigma_hat))
  
  S_half=solve(expm::sqrtm(S))
  Sig_half=expm::sqrtm(Sigma_hat)
  
  # Transform dataset, see Beran & Srivastava (1985)
  
  
  X_trans = X%*%S_half%*%Sig_half
  colnames(X_trans)=colnames(X)
  
  # Collect arguments
  arguments=.object$Information$Arguments
  
  # Calculate reference distribution
  ref_dist=lapply(replicate(.runs, X_trans[sample(1:nObs,replace=TRUE),], simplify = FALSE),function(x){
    # Est_temp=csem(.data =x,
    #               .model = .object$Information$Model,
    #               .approach_cf =  .object$Information$Approach_cf,
    #               .approach_paths =.object$Information$Path_approach,
    #               .approach_weights= .object$Information$Weight_approach,
    #               .disattenuate = .object$Information$Disattenuate,
    #               .dominant_indicators =.object$Information$Dominant_indicators, 
    #               .estimate_structural =.object$Information$Estimate_structural,
    #               .ignore_structural_model =.object$Information$Ignore_structural_model,
    #               .iter_max = .object$Information$Number_max_iterations,
    #               .PLS_mode =.object$Information$PLS_modes,
    #               .PLS_weight_scheme_inner = .object$Information$PLS_inner_weighting_scheme,
    #               .tolerance =  .object$Information$Tolerance)
   
    arguments[[".data"]] <- x
    
    Est_temp=do.call(csem,arguments)
    
    status_code=status(Est_temp)
    
    # if it is controlled for inadmissible
    if(.dropInadmissibles){
      if(0 %in% status_code){
        c(dG = dG(.matrix1=Est_temp$Estimates$Indicator_VCV,.matrix2=fitted(Est_temp)),
          SRMR = SRMR(Est_temp),
          dL=dL(.matrix1=Est_temp$Estimates$Indicator_VCV,.matrix2=fitted(Est_temp))) 
      }
      # else, i.e., dropInadmissible == FALSE
    }else{ 
      c(dG = dG(.matrix1=Est_temp$Estimates$Indicator_VCV,.matrix2=fitted(Est_temp)),
        SRMR = SRMR(Est_temp),
        dL=dL(.matrix1=Est_temp$Estimates$Indicator_VCV,.matrix2=fitted(Est_temp)))
    }
  })
  ref_dist_matrix=do.call(cbind,ref_dist)
  critical_value=apply(ref_dist_matrix,1,quantile,1-.alpha)
  
  
  if(length(.alpha)>1){
    decision = apply(critical_value,1,function(x){teststat<x})
  }
  if(length(.alpha)==1){
    decision = teststat<critical_value
  }
  
  list(Test_statistic=teststat,Critical_value=critical_value, Decision=decision, Number_admissibles=ncol(ref_dist_matrix)) 
  
}
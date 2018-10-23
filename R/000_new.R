### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


#' AVE
#'
#' Computes the average variance extracted \insertCite{Raykov1997}{cSEM} based on the standardized loadings.
#'
#' @usage AVE(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @references 
#' \insertAllCited{}
#'
#'
#' @export
#'
AVE=function(.object=args_default()$.object,
             .only_common_factors=args_default()$.only_common_factors){
  construct_names=names(.object$Information$Model$construct_type)
  
  # Extract loadings
  L=.object$Estimates$Loading_estimates

  AVEs=sapply(construct_names, function(x){
    lam=c(L[x,L[x,]!=0])
    ave=sum(lam^2)/(sum(lam^2)+sum(1-lam^2))
    ave
    
  })
  
  names(AVEs)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    AVEs[co_names]=NULL
  }
  
  return(AVEs) 
}



#' CR
#'
#' Computes the composite reliability \insertCite{Raykov1997}{cSEM} based on standardized loading.
#'
#' @usage CR(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @references 
#' 
#' \insertAllCited{}
#'
#' @export
#'
CR=function(.object=args_default()$.object,
            .only_common_factors=args_default()$.only_common_factors){
  construct_names=names(.object$Information$Model$construct_type)
  
  # Extract loadings
  L=.object$Estimates$Loading_estimates
  
  
  # Calculate CR for all constructs
  CRs=sapply(construct_names, function(x){
    lam=c(L[x,L[x,]!=0])
    cr=sum(lam)^2/(sum(lam)^2+sum(1-lam^2))
    cr
    
  })
  
  names(CRs)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    CRs[co_names]=1
  }
  
  return(CRs)
}

#' Cronbach_alpha
#'
#' Computes Cronbach's alpha \insertCite{Cronbach1951}{cSEM} based on the correlation matrix 
#'
#' @usage Cronbach_alpha(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @references 
#' \insertAllCited{}
#'
#' @export
#'
Cronbach_alpha=function(.object=args_default()$.object,
            .only_common_factors=args_default()$.only_common_factors){
  construct_names=names(.object$Information$Model$construct_type)
  
  # calculate Cronbach's alpah for all constructs
  alphas=sapply(construct_names,function(x){
    relevant_indicators=colnames(.object$Information$Model$measurement[x,.object$Information$Model$measurement[x,]!=0,drop=FALSE])
    S_relevant=.object$Estimates$Indicator_VCV[relevant_indicators,relevant_indicators]
    alpha_temp=psych::alpha(S_relevant,delete=FALSE,na.rm=FALSE)
    alpha_temp$total$raw_alpha
    })
  
  names(alphas)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    alphas[co_names]=NULL
  }
  
  return(alphas)
}

#' Fornell_Larcker
#'
#' Computes the Fornell-Larcker criterion \insertCite{Fornell1981}{cSEM}.
#'
#' @usage Fornell_Larcker(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'@references
#'
#' \insertAllCited{}
#'
#' @export
#'
Fornell_Larcker=function(.object=args_default()$.object,
                                 .only_common_factors=args_default()$.only_common_factors){

  # Calculate the average variance extracted
  average_variance_extracted=AVE(.object,.only_common_factors=FALSE)

  .object

  construct_names=names(.object$Information$Model$construct_type)

  construct_vcv_sq=.object$Estimates$Construct_VCV^2
  diag(construct_vcv_sq)=average_variance_extracted[colnames(construct_vcv_sq)]
stop('Not implemented yet')

}


#' Estimates the distance between multiple groups. 
#'
#' 
#' The geodesic distance (dG) and the Euclidean distance (dE) is used
#' to assess group differences. Permutation is used to generate a 
#' reference distribtuion. 
#' 
#' @usage testMGDnew(
#'  .object             = args_default()$.object,
#'  .alpha              = args_default()$.alpha,
#'  .drop_inadmissibles = args_default()$.drop_inadmissibles,
#'  .runs               = args_default()$.runs,
#'  .saturated          = args_default()$.saturated,
#'  .type_vcv           = args_default()$.type_vcv
#'  ) 
#' 
#' @inheritParams csem_arguments
#' 
#' @inherit csem_test return
#'
#' @seealso [cSEMResults]
#'
#' @examples
#' \dontrun{
#' require(cSEM)
#' data(satisfaction)
#'
#' model <- "
#' # Structural model
#' QUAL ~ EXPE
#' EXPE ~ IMAG
#' SAT  ~ IMAG + EXPE + QUAL + VAL
#' LOY  ~ IMAG + SAT
#' VAL  ~ EXPE + QUAL
#'
#' # Measurement model
#'
#' EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
#' IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
#' LOY  =~ loy1  + loy2  + loy3  + loy4
#' QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
#' SAT  <~ sat1  + sat2  + sat3  + sat4
#' VAL  <~ val1  + val2  + val3  + val4
#' "
#' 
#' listData <- list(satisfaction[-3,], satisfaction[-5, ], satisfaction[-10, ])
#' out.cSEM <- csem(listData, model) 
#'
#' testMGDnew(.object = out.cSEM, .runs = 20, .parallel = TRUE, .type_vcv= 'construct')
#' }
#'
#' @export
#'
testMGDnew <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .runs                  = args_default()$.runs,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .show_progress         = args_default()$.show_progress
){
  
  ### Checks and errors ========================================================
  ## Check if cSEMResults object
  # if(class(.object) != "cSEMResults") {
  #   stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  # }
  
  # ## Check if .object contains estimates for at least two groups.
  # if(attr(.object, "single") == TRUE) {
  #   stop("At least two groups required.", call. = FALSE)
  # }
  
  ## Check if any of the group estimates are inadmissible
  if(!all(sapply(.object, function(x) sum(verify(x)) == 0))) {
    stop("Initial estimation results for at least one group are inadmissible.\n", 
         "See `lapply(.object, verify)` for details.",
         call. = FALSE)
  }
  
  # Verstehe ich nicht
  # Check if data for different groups is identical
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    stop("At least two groups are identical.", call. = FALSE)
  }
  
  ### Calculation===============================================================
  ## 1. Compute the test statistics
  teststat <- c(
    "dG" = calculateDistance(.matrices = fit(.object = .object,
                                             .saturated = .saturated,
                                             .type_vcv =.type_vcv),
                             .distance = "geodesic"),
    "dL" = calculateDistance(.matrices = fit(.object = .object,
                                             .saturated= .saturated,
                                             .type_vcv=.type_vcv),
                             .distance = "squared_euclidian"))
  
  ## 2. Permuation
  # Put data in a list
  org_data_list<- lapply(.object, function(x) x$Information$Data)
  
  # Collect initial arguments (from the first object)
  arguments <- .object[[1]]$Information$Arguments
  
  # 
  # list containing the results
  ref_dist = list()
  counter=1
  total_iterations=0
  
  repeat{
    # permutate data
    X_temp=permutateData(.matrices=org_data_list)
    
    # Replace the old dataset by the new one
    arguments[[".data"]] <- X_temp
    
    # Set .id etl kann dies auch ausserhalb von repeat
    arguments[[".id"]] <- "permID"
    
    # Estimate model
    Est_temp <- do.call(csem, arguments)   
    
    # Check status
    status_code <- verify(Est_temp)
    
    if(.handle_inadmissibles == 'drop' | .handle_inadmissibles == 'replace'){
      if(sum(unlist(status_code)) == 0){
        
        ref_dist[[counter]]= c('dG' = calculateDistance(.matrices = fit(.object = Est_temp, 
                                               .saturated  = .saturated,
                                               .type_vcv=.type_vcv), 
                                               .distance = "geodesic"),
                              'dL' = calculateDistance(.matrices = fit(.object = Est_temp, 
                                                                       .saturated  = .saturated,
                                                                       .type_vcv=.type_vcv),
                                                       .distance = "squared_euclidian"))
        
        counter=counter+1
        
       }else if(sum(unlist(status_code)) != 0 & .handle_inadmissibles == 'drop'){
        # if(.handle_inadmissibles == 'drop')#{
        ref_dist[[counter]]= NULL
        counter=counter+1
        }
      } else if(.handle_inadmissibles == 'ignore') { 
        
        ref_dist[[counter]]= c('dG' = calculateDistance(.matrices = fit(.object = Est_temp, 
                                                                       .saturated  = .saturated,
                                                                       .type_vcv=.type_vcv), 
                                                       .distance = "geodesic"),
                              'dL' = calculateDistance(.matrices = fit(.object = Est_temp, 
                                                                       .saturated  = .saturated,
                                                                       .type_vcv=.type_vcv),
                                                       .distance = "squared_euclidian"))
        counter=counter+1
      }
    
      total_iterations=total_iterations+1  
      
    # Break repeat loop if the necessary number of runs was succesful or 
      # 10'000 iterations have been done without sucess.
      if((counter - 1) == .runs) {
        break
      } else if(total_iterations == 10000) { 
        stop("Not enough admissible result.", call. = FALSE)
      }
    
    
  } #end repeat 
  

  ## Compute critical values 
  ref_dist_matrix <- do.call(cbind, ref_dist)
  critical_value  <- matrix(apply(ref_dist_matrix, 1, quantile, 1-.alpha), 
                            ncol = length(teststat),
                            dimnames = list(paste(.alpha*100, sep = "","%"), 
                                            names(teststat))
  )
  
  if (length(.alpha) > 1) {
    decision <- t(apply(critical_value, 1, function(x) {teststat < x}))
  }
  
  if (length(.alpha) == 1) {
    decision <- teststat < critical_value
  }
  
  # Return output
  out <- list(
    "Test_statistic"     = teststat,
    "Critical_value"     = critical_value, 
    "Decision"           = decision, 
    "Number_admissibles" = ncol(ref_dist_matrix),
    "Total_runs"         = total_iterations
  ) 
  
  class(out) <- "cSEMTestMGD"
  return(out)
}  
  
  



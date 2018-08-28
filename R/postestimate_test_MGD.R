#' Test for multi group differences (MGD)
#'
#' What it does (TODO).
#' 
#' More details here (TODO).
#' 
#' @usage testOH(
#'  .object             = args_default()$.model,
#'  .alpha              = args_default()$.alpha,
#'  .drop_inadmissibles = args_default()$.drop_inadmissibles,
#'  .parallel           = args_default()$.parallel,
#'  .runs               = args_default()$.runs,
#'  .show_progress      = args_default()$.show_progress
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
#' # still to implement
#' }
#'
#' @export

testMGD <- function(
  .object             = args_default()$.model,
  .alpha              = args_default()$.alpha,
  .drop_inadmissibles = args_default()$.drop_inadmissibles,
  .parallel           = args_default()$.parallel,
  .runs               = args_default()$.runs,
  .show_progress      = args_default()$.show_progress
  ){

  ### Checks and errors ========================================================
  ## Check if .object contains estimates for at least two groups.
  if(attr(.object, "single") == TRUE) {
    stop("At least two groups required.", call. = FALSE)
  }
  
  ## Check if any of the group estimates is inadmissible
  if(!all(sapply(.object, function(x) sum(verify(x)) == 0))) {
    stop("Results for at least one group are inadmissible.", call. = FALSE)
  }
  
  # Check if data for different groups is identical
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    stop("Identical data sets.", call. = FALSE)
  }
  
  ### Calculation===============================================================
  ## 1. Compute the test statistics
  teststat <- c(
    "dG" = calculateDistance(.matrices = lapply(.object, fit), 
                             .distance = "geodesic"),
    "dL" = calculateDistance(.matrices = lapply(.object, fit), 
                             .distance = "squared_euclidian")
    )
  
  ## 2. Permuation
  # Put data in a list
  listMatrices <- lapply(.object, function(x) x$Information$Data)
  
  # Collect initial arguments
  arguments <- .object[[1]]$Information$Arguments
  
  # Set .id
  arguments[[".id"]] <- "permID"
  
  # Results
  permEstimates <- list()
  
  ##  WITHOUT PARALLELIZATION --------------------------------------------------
  if(.parallel == FALSE) {
    if(.show_progress == TRUE) {
      # Progress bar
      pb <- txtProgressBar(min = 0, max = .runs, style = 3)
    }
    
    for(iPerm in 1:.runs){
      # Permutate data
      permData <- permutateData(listMatrices)
      # Replace .data 
      arguments[[".data"]] <- permData
      # Estimate using permutated data
      Est_tmp <- do.call(csem, arguments)
      # Check if all estimates produce admissible results
      status_code <- sapply(Est_tmp, verify)
      
      # Drop potential inadmissibles if required
      if(.drop_inadmissibles){
        ## If no inadmissibles exists continue as usual
        if(all(sapply(.object, function(x) sum(verify(x)) == 0))){
          
          permEstimates[[iPerm]] <- c(
            dG = calculateDistance(lapply(Est_tmp, fit), .distance = "geodesic"),
            dL = calculateDistance(lapply(Est_tmp, fit), .distance = "squared_euclidian")
          )
        } else {
          # Set to NULL
          permEstimates[[iPerm]] <- NULL
        }
        # else, i.e., dropInadmissible == FALSE
      } else {
        
        permEstimates[[iPerm]] <- c(
          dG = calculateDistance(lapply(Est_tmp, fit), .distance = "geodesic"),
          dL = calculateDistance(lapply(Est_tmp, fit), .distance = "squared_euclidian")
        )
      }
      if(.show_progress == TRUE){
        # Update Progress bar
        setTxtProgressBar(pb, iPerm)
      }
    }
    if(.show_progress==TRUE){
      # Close Progress bar
      close(pb)
    }
  }
  
  ## With PARALLEL VERSION -----------------------------------------------------
  if(.parallel == TRUE) {
    if(.show_progress==TRUE){
      print("Progress bar is not available for parallel computing.")
    }
    if(!requireNamespace("doSNOW")){
      stop("cSEM requires the doSNOW package")
    } 
    if(!requireNamespace("parallel")){
      stop("cSEM requires the parallel package")
    } 
    # Prepare parallelization
    core= detectCores()
    cl <- makeCluster(core)
    registerDoParallel(cl)
    
    permEstimates <- foreach(iPerm = 1:.runs) %dopar% {
                               # permutate data
                               permData <- permutateData(listMatrices)
                               # set Data
                               arguments[[".data"]] <- permData
                               # estimate 
                               Est_tmp <- do.call(csem, arguments)
                               # status codes
                               status_code=lapply(Est_tmp, verify)
                               
                               # if it is controlled for inadmissible
                               if(.drop_inadmissibles){
                                 if(!(FALSE %in% sapply(status_code, is.null))){
                                   c(dG = calculateDistance(
                                     .matrices =
                                       lapply(Est_tmp, fit),
                                     .distance = "geodesic"
                                   ),
                                   dL = calculateDistance(
                                     .matrices =
                                       lapply(Est_tmp, fit),
                                     .distance = "squared_euclidian"
                                   )
                                   )
                                 }else{
                                   # Set null if status code is false
                                   NULL
                                 }
                                 # else, i.e., dropInadmissible == FALSE
                               }else{
                                 c(
                                   dG = calculateDistance(
                                     .matrices = lapply(Est_tmp, fit),
                                     .distance = "geodesic"
                                   ),
                                   dL = calculateDistance(
                                     .matrices = lapply(Est_tmp, fit),
                                     .distance = "squared_euclidian"
                                   )
                                 )
                                 
                               }
                             }
    # Stop Cluster
    stopCluster(cl)
  }
  
  
  ref_dist       <- do.call(cbind, permEstimates)
  critical_value <- matrix(apply(ref_dist, 1, quantile, 1-.alpha), 
                           ncol = length(teststat), 
                           dimnames = list(paste(.alpha*100, sep = "","%"),
                                           names(teststat)))
  
  if(length(.alpha) > 1){
    decision <- t(apply(critical_value, 1, function(x){
      ifelse(teststat > x, TRUE, FALSE)
      }))
  }
  if(length(.alpha) == 1){
    decision <- ifelse(teststat > critical_value, TRUE, FALSE)
  }
  out <- list(
    "Hypothesis"         = paste0("H0: No significant difference between groups."),
    "Test_statistic"     = teststat, 
    "Critial_value"      = critical_value,
    "Decision"           = decision, 
    "Number_admissibles" = ncol(ref_dist))
  
  # define return class
  class(out) <- "cSEMTest"
  return(out)
}

permutateData <- function(.matrices = NULL){
  
  ### Checks and errors ========================================================
  ## Check if list and at least of length 2
  if(!is.list(.matrices) && length(.matrices) < 2) {
    stop("`.matrices` must be a list of at least length two.", call. = FALSE)
  }
  
  # ## Check if column names are identical
  # if (FALSE %in% sapply(.matrices,function(x) {
  #   identical(colnames(x), colnames(.matrices[[1]]))})) {
  #   stop("`.matrices` must have the same colnames.", call. = FALSE)
  # }
  
  ### Permutation ==============================================================
  
  # combine data
  combinedData <- do.call(rbind, .matrices)
  
  # create ID
  ID <- rep(1:length(.matrices),lengths(.matrices)/ncol(.matrices[[1]]))
  
  # add permID
  permData <- cbind(combinedData, permID = sample(ID))
  
  # If we want to return a list
  # l=lapply(1:length(.matrices),function(x){temp=permData[permData[,'permID']==x,]
  # temp[,-ncol(temp),drop=FALSE]})
  # return(l)
  
  return(permData)
}

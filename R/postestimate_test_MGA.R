#' Test for group differences
#'
#' What it does.
#' 
#' More details here.
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

testOverallMGA <- function(.object=args_default()$.model,
                           .drop_inadmissibles=args_default()$.drop_inadmissibles,
                           .alpha=args_default()$.alpha,
                           .runs=args_default()$.runs,
                           ...){
  
  # test if attribute "single" =  TRUE
  if(attr(.object, "single") == TRUE) {
    stop(".object cannot be single. 
         Overall test for group differences not possible.",   call. = FALSE)
  }
  
  # Check if .object is admissible
  if(FALSE %in% sapply(lapply(.object, verify), is.null)){
    stop("Initial estimation is inadmissible.", call. = FALSE)
  }
  
  # Check if there are identical data sets
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    stop("there are identical data sets.", call. = FALSE)
  }
  
  # 1: calculate test statistic
  teststat = c(dG = calculateDistance(.matrices = lapply(.object, fit), 
                                       .distance="geodesic"),
               dL = calculateDistance(.matrices = lapply(.object, fit), 
                                       .distance="squared_euclidian"))
  
  # 2: permutation procedure
  listMatrices <- list()
  for(iData in 1:length(.object)){
    listMatrices[[iData]] <- .object[[iData]]$Information$Data
  }

  # Collect initial arguments
  arguments=.object[[1]]$Information$Arguments
  arguments[[".id"]] <- "permID"
  
  # Results
  permEstimates <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = .runs, style = 3)
  
  for(iPerm in 1:.runs){
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
        permEstimates[[iPerm]] <-
          c(
            dG = calculateDistance(
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
        permEstimates[[iPerm]] <- NULL
      }
      # else, i.e., dropInadmissible == FALSE
    }else{
      permEstimates[[iPerm]] <-
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
    # Update Progress bar
    setTxtProgressBar(pb, iPerm)
  }
  # Close Progress bar
  close(pb)
  
  ref_dist <- do.call(cbind, permEstimates)
  critical_value <- matrix(apply(ref_dist, 1, quantile, 1-.alpha), 
                           ncol = length(teststat), 
                           dimnames = list(paste(.alpha*100, sep = "","%"),
                                           names(teststat)))
  
  if(length(.alpha)>1){
    decision = t(apply(critical_value,1,function(x){teststat<x}))
  }
  if(length(.alpha)==1){
    decision = teststat<critical_value
  }
  out <- list(Test_statistic = teststat, 
              Critial_value = critical_value,
              Decision = decision, 
              Number_admissibles = ncol(ref_dist))
  
  # define return class
  class(out) <- "cSEMTestResults"
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

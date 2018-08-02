#' @title Test for group differences
#'
#' @description x of this function.
#' 
#' @details Whaaaaaaaaaaaaaaaats up.  \deqn{sqrt{9} + b}
#' 
#' @usage testOverallMGA(.object=args_default()$.model,
#' .dropInadmissibles=args_default()$.dropInadmissibles,
#' .alpha=args_default()$.alpha,
#' .runs=args_default()$.runs,
#' ...)
#' 
#' @inheritParams csem_arguments
#' 
#' @inherit csem_testresults return
#' 
#' @references
#'   \insertAllCited{}
#'
#' @seealso [csem], [foreman]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export
testOverallMGA <- function(.object=args_default()$.model,
                           .dropInadmissibles=args_default()$.dropInadmissibles,
                           .alpha=args_default()$.alpha,
                           .runs=args_default()$.runs,
                           ...){
  # Check if .id variable is set
  if(!is.null(.object[[1]]$Information$Arguments$.id)){
    stop("No .id variable set. Overall test for group differences not possible.",
         call. = FALSE)
  }
  # Check if .object is admissible
  if(FALSE %in% sapply(lapply(.object, status), is.null)){
    stop("Initial estimation is inadmissible.", call. = FALSE)
  }
  
  # 1: calculate test statistic
  teststat = c(dG = calculateArithmDistance(.corMatrices = lapply(.object, fitted), 
                                       .distance="geodesic"),
               dL = calculateArithmDistance(.corMatrices = lapply(.object, fitted), 
                                       .distance="euclidean"))
  
  # 2: permutation procedure
  listMatrices <- list()
  for(iData in 1:length(.object)){
    listMatrices[[iData]] <- .object[[iData]]$Information$Data
  }
  # combine data
  totalData <- data.frame(do.call(rbind, listMatrices))
  

  # Collect initial arguments
  arguments=.object[[1]]$Information$Arguments
  arguments[[".id"]] <- "permID"
  
  # Results
  permEstimates <- list()
  
  for(iPerm in 1:.runs){
    # create ID
    ID <- rep(1:length(listMatrices),
              lengths(listMatrices)/ncol(listMatrices[[1]]))
    # random ID
    permData <- cbind(totalData, permID = as.character(sample(ID)))
    # set Data
    arguments[[".data"]] <- permData
    # estimate 
    Est_tmp <- do.call(csem, arguments)
    # status codes
    status_code=lapply(Est_tmp, status)
    
    # if it is controlled for inadmissible
    if(.dropInadmissibles){
      if(!(FALSE %in% sapply(status_code, is.null))){
        permEstimates[[iPerm]] <-
          c(
            dG = calculateArithmDistance(
              .corMatrices =
                lapply(Est_tmp, fitted),
              .distance = "geodesic"
            ),
            dL = calculateArithmDistance(
              .corMatrices =
                lapply(Est_tmp, fitted),
              .distance = "euclidean"
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
          dG = calculateArithmDistance(
            .corMatrices = lapply(Est_tmp, fitted),
            .distance = "geodesic"
          ),
          dL = calculateArithmDistance(
            .corMatrices = lapply(Est_tmp, fitted),
            .distance = "euclidean"
          )
        )
      
    }
  }  
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


#' Estimates the distance between multiple groups. 
#'
#' 
#' The geodesic distance (dG) and the Euclidean distance (dE) is used
#' to assess group differences. Permutation is used to generate a 
#' reference distribtuion. 
#' 
#' @usage testMGD(
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
#' testMGD(.object = out.cSEM, .runs = 20, .parallel = TRUE)
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
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  ## Check if .object contains estimates for at least two groups.
  if(attr(.object, "single") == TRUE) {
    stop("At least two groups required.", call. = FALSE)
  }
  
  ## Check if any of the group estimates are inadmissible
  if(!all(sapply(.object, function(x) sum(verify(x)) == 0))) {
    stop("Initial estimation results for at least one group are inadmissible.\n", 
         "See `lapply(.object, verify)` for details.",
         call. = FALSE)
  }
  
  # Check if data for different groups is identical
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    stop("At least two groups are identical.", call. = FALSE)
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
  
  ##  PREPARE PERMUTAITON --------------------------------------------------
  opts <- NULL
  # Progress bar
  if(.show_progress){
    # Progress bar
    pb <- txtProgressBar(min = 0, max = .runs, style = 3)
    if(.parallel){
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    }
  }
  
  # Parallel
  if(.parallel == TRUE){
    # Initiate Cluster
    nprocs = parallel::detectCores()
    cl <- parallel::makeCluster(nprocs)
    doSNOW::registerDoSNOW(cl)
    # DoParallel
    permEstimates <- foreach::foreach(iPerm = 1:.runs, .options.snow = opts) %dopar% {
      permutationProcedure(.object = .object,
                           .listMatrices = listMatrices, 
                           .arguments = arguments, 
                           .drop_inadmissibles = .drop_inadmissibles)
    }
    parallel::stopCluster(cl)
  }else{
    # Sequential
    permEstimates <- list()
    for(iPerm in 1:.runs){
      permEstimates[[iPerm]] <- permutationProcedure(.object = .object,
                                                     .listMatrices = listMatrices, 
                                                     .arguments = arguments,
                                                     .drop_inadmissibles = .drop_inadmissibles)
      # Update progress bar
      if(.show_progress){
        setTxtProgressBar(pb, iPerm)
      }
    }
  }

  # Close progress bar
  if(.show_progress){
    close(pb)
  }
  
  # Calculate Estimates
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
    "Test_statistic"     = teststat, 
    "Critical_value"     = critical_value,
    "Decision"           = decision, 
    "Number_admissibles" = ncol(ref_dist))
  
  # define return class
  class(out) <- "cSEMTestMGD"
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

permutationProcedure <- function(.object, .listMatrices, .arguments, .drop_inadmissibles){
  # Permutate data
  permData <- permutateData(.listMatrices)
  # Replace .data 
  .arguments[[".data"]] <- permData
  # Estimate using permutated data
  Est_tmp <- do.call(csem, .arguments)
  # Check if all estimates produce admissible results
  status_code <- sapply(Est_tmp, verify)
  # Drop potential inadmissibles if required
  if(.drop_inadmissibles){
    ## If no inadmissibles exists continue as usual
    if(all(sapply(.object, function(x) sum(verify(x)) == 0))){
      return(c(
        dG = calculateDistance(lapply(Est_tmp, fit), .distance = "geodesic"),
        dL = calculateDistance(lapply(Est_tmp, fit), .distance = "squared_euclidian")))
    } else {
      # return NULL
      return(NULL)
    }
    # else, i.e., dropInadmissible == FALSE
  } else {
    return(c(
      dG = calculateDistance(lapply(Est_tmp, fit), .distance = "geodesic"),
      dL = calculateDistance(lapply(Est_tmp, fit), .distance = "squared_euclidian")))
  }
}


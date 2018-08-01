#' Test Multigroup
#' @export
#'


# @Distance functions
calculateArithmDistance <- function(.corMatrices, .distance="geodesic"){
  # number of groups to compare
  groups <- length(.corMatrices)
  # result vector
  results <- c()
  # index
  index <- 1
  for (g in 2:(groups)) {
    for(h in 1:(g - 1)){
      
      switch(.distance,
             # Geodesic distance
             geodesic = {
               # Eigen values
               Eigen <- eigen(solve(.corMatrices[[h]]) %*% 
                                       .corMatrices[[g]])
               # calc distance
               results[index]  <-  sum((log(Eigen$values))^2)
               
               # results[index] <- 2 * dG(.corMatrices[[h]], .corMatrices[[g]])
               # update index
               index <- index + 1 
             },
             # Euclidean Distance
             euclidean = {
               # Estimate euclidean Distance
               results[index]  <- sqrt((sum((.corMatrices[[h]] -
                                               .corMatrices[[g]])^2))/
                                         (ncol(.corMatrices[[h]])*
                                            (ncol(.corMatrices[[h]]) + 1))) 
               index <- index + 1 
             },
             # Add other distance measures here
             
             # Distance does not exist
             {
               stop("This distance measure is not defined.")
             }
      )
    }
  }
  res <- (1 / (groups*(groups-1))) * sum(results)
  return(res)
}

#' Test statistic to compare group differences. 
#'
#' @param .objc a cSEM object
#' @param .runs the number 
#' @return The test statistic, the critical value, a logical value with 
#' the decision, and the number of admissible results.
#' @examples
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
  
  return(list(Test_statistic = teststat, 
              Critial_value = critical_value,
              Decision = decision, 
              Number_admissibles = ncol(ref_dist)))
}


# Test ---------------------------

model <- "
# Structural model
EXPE ~ IMAG
QUAL ~ EXPE
VAL  ~ EXPE + QUAL
SAT  ~ IMAG + EXPE + QUAL + VAL

LOY  ~ IMAG + SAT

# Measurement model

IMAG <~ imag1 + imag2 + imag3
EXPE <~ expe1 + expe2 + expe3
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
VAL  <~ val1  + val2  + val3
SAT  <~ sat1  + sat2  + sat3  + sat4
LOY  <~ loy1  + loy2  + loy3  + loy4
"

sat.plspm <- read.csv2("C:/Users/Michael Klesel/Dropbox/R Projekte/R cSEM/satisfaction.csv")

require(cSEM)
data(satisfaction)
a <- csem(.data = sat.plspm, .model = model, .id = "gender")
lapply(a, status)


testOverallMGA(.object = a, .runs = 20, .alpha = c(0.05, 0.01))



b <- csem(.data = satisfaction, .model = model, .PLS_weight_scheme_inner = "path")

Test_for_overall_model_fit(b)


require(matrixpls)
c <- matrixpls(S=cor(satisfaction), model = model)
sum(fitted.matpls(c)[rownames(fitted(b)),colnames(fitted(b))]-fitted(b))



# MGA
a_test <- testOverallMGA(a, .10)
a_test


testOverallMGA(a, 5)


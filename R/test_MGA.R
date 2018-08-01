#' Test Multigroup
#'
#'


# Calculate the arithmetic distance between groups
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

# Check if .object object has an id variable
checkMGA <- function(.object){
  if(is.null(.object[[1]]$Information$Arguments$.id)){
    return(FALSE)
  }
  return(TRUE)
}

testOverallMGA <- function(.object, .runs){
  if(!checkMGA(.object)){
    stop("No .id variable set. Overall test for group differences not possible.", call. = FALSE)
  }
  
  # 1: Overall distance
  geoDistance <- calculateArithmDistance(lapply(.object, fitted))
  
  # 2: permutation procedure
  # extract data
  listMatrices <- list()
  for(iData in 1:length(.object)){
    listMatrices[[iData]] <- .object[[iData]]$Information$Data
  }
  # combine data
  totalData <- data.frame(do.call(rbind, listMatrices))
  
  # Results
  permEstimates <- c()
  for(iPerm in 1:.runs){
    # create ID
    ID <- rep(1:length(listMatrices),lengths(listMatrices)/ncol(listMatrices[[1]]))
    # random ID
    permData <- cbind(totalData, permID = as.character(sample(ID)))
    
    # Collect arguments
    arguments=.object[[1]]$Information$Arguments
    arguments[[".data"]] <- permData
    arguments[[".id"]] <- "permID"
    
    permOut <- do.call(csem, arguments)
    
    permEstimates[iPerm] <- calculateArithmDistance(lapply(permOut, fitted))
  }
  return(c(geoDistance, quantile(permEstimates)))
}

# 
# 
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



testOverallMGA(a, 20)




b <- csem(.data = satisfaction, .model = model, .PLS_weight_scheme_inner = "path")

Test_for_overall_model_fit(b)


require(matrixpls)
c <- matrixpls(S=cor(satisfaction), model = model)
sum(fitted.matpls(c)[rownames(fitted(b)),colnames(fitted(b))]-fitted(b))



# MGA
a_test <- testOverallMGA(a, .10)
a_test


testOverallMGA(a, 5)


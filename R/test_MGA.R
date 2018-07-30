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

extractIndicator_VCV <- function(.cSEM){
  groups <- length(names(.cSEM))
  corMatrices <- list()
  for(iGroups in 1:groups){
    # TODO: implied
    corMatrices[[iGroups]] <- .cSEM[[iGroups]]$Estimates$Indicator_VCV
  }
  names(corMatrices) <- names(.cSEM)
  return(corMatrices)
}


# Check if .cSEM object has an id variable
checkMGA <- function(.cSEM){
  if(is.null(.cSEM[[1]]$Information$Arguments$.id)){
    return(FALSE)
  }
  return(TRUE)
}

testOverallMGA <- function(.cSEM, .permutations){
  if(!checkMGA(.cSEM)){
    stop("No .id variable set. Overall test for group differences not possible.")
  }
  # 1: Overall distance
  geoDistance <- calculateArithmDistance(extractIndicator_VCV(.cSEM))
  
  # extract data
  listMatrices <- list()
  for(iData in 1:length(.cSEM)){
    listMatrices[[iData]] <- .cSEM[[iData]]$Information$Data
  }
  # combine data
  totalData <- data.frame(do.call(rbind, listMatrices))
  
  # Results
  permEstimates <- c()
  for(iPerm in 1:.permutations){
    # create ID
    ID <- rep(1:length(listMatrices),lengths(listMatrices)/ncol(listMatrices[[1]]))
    # random ID
    permData <- cbind(totalData, permID = as.character(sample(ID)))
    
    # Collect arguments
    arguments=.cSEM[[1]]$Information$Arguments
    arguments[[".data"]] <- permData
    arguments[[".id"]] <- "permID"
    
    permOut <- do.call(csem, arguments)
    
    permEstimates[iPerm] <- calculateArithmDistance(extractIndicator_VCV(permOut))
  }
  return(c(geoDistance, permEstimates))
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
SAT  =~ sat1  + sat2  + sat3  + sat4
LOY  =~ loy1  + loy2  + loy3  + loy4
"


require(cSEM)
data("satisfaction")
b <- csem(.data = cSEM::satisfaction, .model = model)
b
listviewer::jsonedit(b, mode="view")

require(plspm)
data("satisfaction")
a <- csem(.data = satisfaction, .model = model, .id = "gender")

# MGA
a_test <- testOverallMGA(a, .10)
a_test


testOverallMGA(a, 5)




listviewer::jsonedit(a, mode="view")

cdata <- permutateData(a)
c <- csem(.data = list(cdata[[1]], cdata[[2]]), .model = model)





a$female$Information$Model

str(a)

calculateArithmDistance(a, "geodesic")
calculateArithmDistance(a, "euclidean")

a_perm <- permutateData(a)
calculateArithmDistance(a_perm, "geodesic")


listviewer::jsonedit(a, mode="view")

a$male$Information$Data
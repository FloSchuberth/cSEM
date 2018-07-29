#' Test Multigroup
#'
#'

# Calculate the arithmetic distance between groups
calculateArithmDistanceOld <- function(.cSEM, .distance="geodesic"){
  # number of groups to compare
  groups <- length(.cSEM)
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
               Eigen <- eigen(solve(.cSEM[[h]]$Estimates$Indicator_VCV) %*% 
                         .cSEM[[g]]$Estimates$Indicator_VCV)
               # calc distance
              results[index]  <-  sum((log(Eigen$values))^2) 
              # update index
              index <- index + 1 
             },
             # Euclidean Distance
             euclidean = {
               # Estimate euclidean Distance
               results[index]  <- sqrt((sum((.cSEM[[h]]$Estimates$Indicator_VCV -
                                               .cSEM[[g]]$Estimates$Indicator_VCV)^2))/
                                         (ncol(.cSEM[[h]]$Estimates$Indicator_VCV)*
                                            (ncol(.cSEM[[h]]$Estimates$Indicator_VCV)
                                             + 1))) 
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
    corMatrices[[iGroups]] <- .cSEM[[iGroups]]$Estimates$Indicator_VCV
  }
  names(corMatrices) <- names(.cSEM)
  return(corMatrices)
}

# permutate Data
permutateData <- function(.cSEM){
  if(!checkMGA(.cSEM)){
    stop("Only one group. Permutation not possible.")
  }
  # Store data in a List
  dataList <- list()
  for(iData in 1:length(.cSEM)){
    dataList[[iData]] <- .cSEM[[iData]]$Information$Data
  }

  dataFlat <- do.call("rbind", dataList)
  # Shuffle
  dataShuffled <- dataFlat[sample(1:nrow(dataFlat)), ]
  # Put together
  totalElements <- sum(lengths(dataList))
  sizeSubLists <- lengths(dataList)
  permData <- list()
  start <- 1
  end <- lengths(dataList)[[1]]
  for(i in 1:length(lengths(dataList))){
    permData[[i]] <-  matrix(t(dataShuffled)[start:end], 
                             ncol= ncol(dataShuffled), 
                             byrow = T, 
                             dimnames = list(c(),
                              colnames(dataFlat)))
    start <- start + lengths(dataList)[i]
    end <- end + lengths(dataList)[i]
  }
  names(permData) <- names(.cSEM)
  return(permData)
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
 
  .cSEM <- a
  
  # 1: Overall distance
  geoDistance <- calculateArithmDistance(extractIndicator_VCV(.cSEM))
  geoDistance
  
  # 2: Permutation
  .permutations <- 5
  
  # Results
  permEstimates <- list()
  
  #for(iPerm in 1:.permutations){
    # Permute data
    permData <- permutateData(.cSEM = .cSEM)
    # Get original arguments
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
    
    # New esimate for each element
    permOut <- cSEM::csem(satisfa, model)
    
    permData_VCV <- list()
    for(iGroup in 1:length(permData)){
      permData_VCV[[1]] <- cSEM::csem(permData[[1]], model)$Estimates$Indicator_VCV
      permData_VCV[[2]] <- cSEM::csem(permData[[2]], model)$Estimates$Indicator_VCV
    }
    
    permDistance <- calculateCorrectionFactors(permData_VCV)
    permDistance
     
  #}

   
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
a
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
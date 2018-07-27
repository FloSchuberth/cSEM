#' Test Multigroup
#'
#'

# Calculate the arithmetic distance between groups
calculateArithmDistance <- function(.cSEM, .distance="geodesic"){
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

# permutate Data
permutateData <- function(.cSEM){
  # TODO: check if group variable (ID) exist
  if(length(.cSEM)==1){
    stop("Only one group. Permutation not possible.")
  }
  # Store data in a List
  dataList <- list()
  for(iData in 1:length(.cSEM)){
    dataList[[iData]] <- .cSEM[[iData]]$Information$Data
  }
  # Put together and unlist
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
  return(permData)
}

testOverallMGA <- function(.cSEM, .permutations){
  if(length(.cSEM) = 1){
    stop("Only one group. Multigroup is not possible.")
  }
 
  # 1: Overall distance
  geoDistance <- calculateArithmDistance(.cSEM)
  geoDistance
  
  # 2: Permutation
  .cSEM <- a
  .permutations <- 5
  
  # Results
  permEstimates <- list()
  
  for(iPerm in 1:.permutations){
    args_needed[[".data"]] <- permutateData(.cSEM)
    
    permEstimates[iPerm] <- lappy(permData, .cSEM$female$Information$Arguments)
    
  }

   
}





# Test ---------------------------


require(cSEM)
require(plspm)
data(satisfaction)

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

a <- csem(.data = satisfaction, .model = model, .id = "gender")
a

a$female$Information$Model

str(a)

calculateArithmDistance(a, "geodesic")
calculateArithmDistance(a, "euclidean")

a_perm <- permutateData(a)
calculateArithmDistance(a_perm, "geodesic")


listviewer::jsonedit(a, mode="view")

a$male$Information$Data
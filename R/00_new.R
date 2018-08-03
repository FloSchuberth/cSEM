# Calculates the Heterotrait-Monotrait factor correlations, see Henseler et al. (2015)


HTMT = function(.object){
  S=.object$Estimates$Indicator_VCV
  
  # HTMT only for common factors
  cf_names=names(.object$Information$Construct_types[.object$Information$Construct_types=="Common factor"])
  
  stop('Not yet implemented') 
  
}

# inspired by matrixpls
convergenceCheck=function(.W_new=W,.W_old=W_iter,.kind='absolute',.tolerance=1e-05){
  
  switch (.kind,
          "absolute" = {max(abs(.W_old - .W_new))<.tolerance},
          "squared" = {max((.W_old-.W_new)^2)<.tolerance},
          "relative" ={max(abs((.W_old[.W_new != 0]-.W_new[.W_new != 0])/
                                 .W_new[.W_new != 0]))<.tolerance} 
  )
}

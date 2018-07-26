# Calculates the Heterotrait-Monotrait factor correlations, see Henseler et al. (2015)


HTMT = function(.object){
  S=.object$Estimates$Indicator_VCV
  
  # HTMT only for common factors
  cf_names=names(.object$Information$Construct_types[.object$Information$Construct_types=="Common factor"])
  
  stop('Not yet implemented') 
  
}
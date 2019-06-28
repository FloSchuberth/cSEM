#' Function that calculates the parameter differences across groups
#'
#' 
#' 
#' @usage parameter_difference(
#' .object=args_default()$.object,
#' .comparison=args_default()$.comparison)
#' 
#' @inheritParams csem_arguments
#' 
#'
#' @seealso [cSEMResults]
#'
parameter_difference=function(.object=args_default()$.object,
                              .comparison=args_default()$.comparison){
  

  # Parse model that indicates which parameters should be compared
  # if no model indicating the comparisons is provided, all parameters are compared
  if(is.null(.comparison)){
  model_comp = .object[[1]]$Information$Model  
  }else{
  model_comp=cSEM:::parseModel(.comparison,.check_errors = F)
  }
  
  # extract different types of constructs
  construct_type = .object[[1]]$Information$Model$construct_type
  
  
  # Create indication matrix for structural coefficients 
  path=.object[[1]]$Estimates$Path_estimates
  path_ind=path
  path_ind[]=0
  path_ind_temp=which(model_comp$structural==1,arr.ind = TRUE)

  if(!(dim(path_ind_temp)[1]==0)){  
    path_ind_temp=cbind(rownames(model_comp$structural)[path_ind_temp[, 'row']],
                        colnames(model_comp$structural)[path_ind_temp[, 'col']])
    
    path_ind[path_ind_temp]=1

    # check whether specified path coefficients occur in the original model
    path_ind_org = .object[[1]]$Information$Model$structural
    if(!all(path_ind_org[path_ind_temp]==1)){
      stop2("Path coefficients specified for comparison are not part of the original model.")
    }
    
  }
  
  # Create indication matrix for loadings
  load=.object[[1]]$Estimates$Loading_estimates
  load_ind=load
  load_ind[]=0
  
  cf_name= names(which(construct_type == 'Common factor'))
  cf_name_comp=intersect(rownames(model_comp$measurement),cf_name)
  
  load_ind_temp = which(model_comp$measurement[cf_name_comp,,drop=FALSE]==1,arr.ind = TRUE)
  
  if(!(dim(load_ind_temp)[1]==0)){
    load_ind_temp = cbind(rownames(model_comp$measurement[cf_name_comp,,drop = FALSE])[load_ind_temp[, 'row']],
                          colnames(model_comp$measurement[cf_name_comp,,drop = FALSE])[load_ind_temp[, 'col']])
    load_ind[load_ind_temp]=1
  
    # check whether specified loadings occur in the original model
    load_ind_org = .object[[1]]$Information$Model$measurement[cf_name,,drop=FALSE]
    if(!all(load_ind_org[load_ind_temp]==1)){
      stop2("Loadings specified for comparison are not part of the original model.")
    }
    }
  
  # Create indication matrix for weights
  weight=.object[[1]]$Estimates$Weight_estimates
  weight_ind=weight
  weight_ind[]=0
  
  co_name= names(which(construct_type == 'Composite'))
  co_name_comp=intersect(rownames(model_comp$measurement),co_name)
  
  weight_ind_temp = which(model_comp$measurement[co_name_comp,,drop=FALSE]==1,arr.ind = TRUE)
  
  if(!(dim(weight_ind_temp)[1]==0)){
    weight_ind_temp = cbind(rownames(model_comp$measurement[co_name_comp,,drop = FALSE])[weight_ind_temp[, 'row']],
                            colnames(model_comp$measurement[co_name_comp,,drop = FALSE])[weight_ind_temp[, 'col']])
    
    weight_ind[weight_ind_temp]=1
    
    # check whether specified weights occur in the original model
    weight_ind_org = .object[[1]]$Information$Model$measurement[co_name,,drop=FALSE]
    if(!all(weight_ind_org[weight_ind_temp]==1)){
      stop2("Weights specified for comparison are not part of the original model.")
    }
    
  }
  

  
  
  # Calculate differences
  # Path coefficients
  matrices_path=lapply(.object,function(x){x$Estimates$Path_estimates})
  
  temp <- utils::combn(matrices_path, 2, simplify = FALSE)
  diff_path=lapply(temp, function(x){
    x[[1]]-x[[2]]
  })
  
  # Loadings
  matrices_load=lapply(.object,function(x){x$Estimates$Loading_estimates})
  temp <- utils::combn(matrices_load, 2, simplify = FALSE)
  diff_load=lapply(temp, function(x){
    x[[1]]-x[[2]]
  })
  
  # Weights
  matrices_weight=lapply(.object,function(x){x$Estimates$Weight_estimates})
  temp <- utils::combn(matrices_weight, 2, simplify = FALSE)
  diff_weight=lapply(temp, function(x){
    x[[1]]-x[[2]]
  })
  
  names(diff_load) <- names(diff_weight) <- names(diff_path) <- sapply(temp, function(x) paste(names(x)[1], 'vs.', names(x)[2]))
  
  vec_and_name=function(.para_diff_matrix,
                        .ind_matrix,
                        .name_matrix,
                        .sep){
    if(!(dim(.name_matrix)[1]==0)){
      temp=c(.para_diff_matrix[.ind_matrix==1])
      # give names
      name=paste(.name_matrix[,1], .sep, .name_matrix[,2])
      names(temp)=name
      temp
    } else {
      return(NA)
    }
  }
  
  difference_path=lapply(diff_path, function(x) vec_and_name(.para_diff_matrix = x,
                               .ind_matrix = path_ind,
                               .name_matrix = path_ind_temp,
                               .sep = '~'))
  
  difference_load=lapply(diff_load, function(x) vec_and_name(.para_diff_matrix = x,
                                                             .ind_matrix = load_ind,
                                                             .name_matrix = load_ind_temp,
                                                             .sep = '=~'))
  
  difference_weight=lapply(diff_weight, function(x) vec_and_name(.para_diff_matrix = x,
                                                               .ind_matrix = weight_ind,
                                                               .name_matrix = weight_ind_temp,
                                                               .sep = '<~'))
  
  # merge list together
  out=mapply(function(x,y,z){c(x,y,z)},
             x=difference_path,y=difference_load,z=difference_weight,
             SIMPLIFY = FALSE)
  
  return(out)
}

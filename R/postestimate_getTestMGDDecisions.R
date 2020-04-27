#' Get decision from cSEMTestMGD object.
#'
#' This function summarizes the decisions of all tests passed in the .object. 
#' 
#' This function is particularily useful to get an overview of all decisions
#' made. 
#' 
#' @usage getTestMGDDecisions(.object)
#'
#' @inheritParams csem_arguments
#' 
#' @return A data.frame with the decisions.
#' 
#' @seealso [csem()], [testMGD()]
#' @export

getTestMGDDecisions <- function(.object){
  
  ## Install plotly if not already installed
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop2(
      "Package `dplyr` required. Use `install.packages(\"dplyr\")` and rerun.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop2(
      "Package `tidyr` required. Use `install.packages(\"tidyr\")` and rerun.")
  }
  
  
  # Chef is a testMGD object is used
  if(!inherits(.object, "cSEMTestMGD")){
    stop2("Object has not the expected class. Please ensure that .object
          is from class cSEMTestMGD")
  }
  
  # Help function --------------------------------------------------------------
  
  # Change CI alpha equivlanet (1-CI)
  changeCItoAlphaEquivalent <- function(.names){
    # .names <- c("99%","95%","90%")
    tmpNames <- readr::parse_number(.names)
    tmpNames <- 100-tmpNames
    tmpNames <- paste0(tmpNames, "%")
    return(tmpNames)
  }
  
  # Klesel ---------------------------------------------------------------------
  
  # Test only provides an overall decision 
  
  if("Klesel" %in% names(.object)){
    klesel <- .object$Klesel$Decision %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      # longer format
      tidyr::pivot_longer(cols = c(dG, dL), 
                          names_to = "distance", 
                          values_to = "decision") %>%
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Klesel", comparison = "overall")
  }
  
  # Chin  ----------------------------------------------------------------------
  
  if("Chin" %in% names(.object)){
    
    # overall decision
    chin1 <- .object$Chin$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::mutate(test = "Chin", comparison = "overall")
    
    chin2 <-  .object$Chin$Decision %>% 
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "correction") %>%
      # Check all comparisons
      dplyr::group_by(correction, alpha) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>% 
      # put results into format 
      tidyr::pivot_longer(cols = (-c(correction, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Chin")
  }
  
  # Sarstedt  ------------------------------------------------------------------
  
  if("Sarstedt" %in% names(.object)){
    
    # overall result
    sarstedt1 <- .object$Sarstedt$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::mutate(test = "Sarstedt", comparison = "overall")
    
    # decision based on single paths
    sarstedt2 <- .object$Sarstedt$Decision %>%
      purrr::modify_depth(2, dplyr::bind_rows) %>%
      purrr::map(dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "correction") %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(correction, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Sarstedt")
  }
  
  # Keil  ----------------------------------------------------------------------
  
  if("Keil" %in% names(.object)){
    
    # overall decision
    keil1 <- .object$Keil$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::mutate(test = "Keil", comparison = "overall")
    
    # decision based on single paths
    keil2 <- .object$Keil$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::group_by(correction, alpha) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(correction, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Keil")
    
  }
  
  # Nitzl  ---------------------------------------------------------------------
  
  if("Nitzl" %in% names(.object)){
    
    # overall decision
    nitzl1 <- .object$Nitzl$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::mutate(test = "Nitzl", comparison = "overall")
    
    # decision based on single paths
    nitzl2 <- .object$Nitzl$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::group_by(correction, alpha) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(correction, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Nitzl")
    
  }
  
  # Henseler  ------------------------------------------------------------------
  
  if("Henseler" %in% names(.object)){
    
    # overall decision
    henseler1 <- .object$Henseler$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::mutate(test = "Henseler", comparison = "overall")
    
    # decision based on single paths
    henseler2 <- .object$Henseler$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "correction") %>%
      dplyr::group_by(correction, alpha) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(correction, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      dplyr::mutate(test = "Henseler")
  }
  
  # CI para  -------------------------------------------------------------------
  
  if("CI_para" %in% names(.object)){
    
    # overall decision
    CIpara1 <- .object$CI_para$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_longer(cols = contains("CI"), 
                          names_to = "type_ci", 
                          values_to = "decision") %>%
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      # Change % as equivlanet to alpha
      dplyr::rename_at(dplyr::vars(contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(test = "CI_para", comparison = "overall")
    
    # decision based on single paths
    CIpara2 <-  .object$CI_para$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(3, ~ dplyr::select(., Name, Decision)) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "type_ci") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "comparison") %>%
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_wider(names_from = Name, values_from = Decision) %>%
      dplyr::group_by(alpha, type_ci) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(type_ci, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      # Change % as equivlanet to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(test = "CI_para")
    
  }
  # CI overlap  ----------------------------------------------------------------
  
  if("CI_overlap" %in% names(.object)){
    
    # overall decision
    CIoverlap1 <- .object$CI_overlap$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_longer(cols = dplyr::contains("CI"), 
                          names_to = "type_ci", 
                          values_to = "decision") %>%
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      # Change % as equivlanet to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(test = "CI_overlap", comparison = "overall")
    
    # decision based on single paths
    CIoverlap2 <-  .object$CI_overlap$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(3, ~ dplyr::select(., Name, Decision)) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "type_ci") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "comparison") %>%
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_wider(names_from = Name, values_from = Decision) %>%
      dplyr::group_by(alpha, type_ci) %>%
      dplyr::summarise_at(dplyr::vars(-comparison),any) %>%
      ungroup %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(type_ci, alpha)),
                          names_to = "comparison", 
                          values_to = "decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = alpha, values_from = decision) %>%
      # Change % as equivlanet to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(test = "CI_overlap")
  }
  
  # Summarize all results ------------------------------------------------------
  
  results <- dplyr::bind_rows(klesel, 
                              chin1, chin2,
                              sarstedt1, sarstedt2,
                              keil1, keil2,
                              nitzl1 , nitzl2,
                              henseler1, henseler2,
                              CIpara1, CIpara2,
                              CIoverlap1, CIoverlap2) %>%
    # Change order of columns
    dplyr::select(test, comparison, dplyr::contains("%"), 
                  correction, type_ci, distance)
  
  return(results)
  }

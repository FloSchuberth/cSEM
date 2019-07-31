#' Resample data 
#'
#' Resample data from a dataset using common resampling methods. 
#' For bootstrap or jackknife resampling, package users usually do not need to 
#' call this function but directly use [resamplecSEMResults()] instead.
#'
#' The function `resampleData()` is general purpose. It simply resamples data 
#' from a dataset according to the resampling method provided 
#' via the `.resample_method` argument and returns a list of resamples. 
#' Currently, `bootstrap`, `jackknife`, `permutation`, and  `cross-validation`
#' (both leave-one-out (LOOCV) and k-fold cross-validation) are implemented. 
#' 
#' The user may provide the dataset to resample from can be either explicitly via the `.data` 
#' argument or implicitly by providing a [cSEMResults] objects to `.object`
#' in which case the original data used in the call that created the 
#' [cSEMResults] object is used for resampling. 
#' If both, a [cSEMResults] object and a dataset via `.data` are provided 
#' the former is ignored. 
#' 
#' As [csem()] accepts a single dataset, a list of datasets as well as datasets
#' that contain a column name used to split the data into groups,
#' the [cSEMResults] object may contain multiple datasets.
#' In this case, resampling is done by dataset or group. Note that depending
#' on the number of datasets/groups provided this computation may be slower
#' as resampling will be repeated for each dataset/group. 
#' 
#' To split data provided via the `.data` argument into groups, the column name or 
#' the column index of the column containing the group levels to split the data 
#'  must be given to `.id`. If data that contains grouping is taken from 
#' a [cSEMResults] object, `.id` is taken from the object information. Hence, 
#' providing  `.id` is redundant in this case and therefore ignored.
#' 
#' The number of bootstrap or permutation runs as well as the number of 
#' cross-validation repetitions is given by `.R`. The default is
#' `499` but should be increased in real applications. See e.g.,
#' \insertCite{Hesterberg2015;textual}{cSEM}, p.380 for recommendations concerning
#' the bootstrap. For jackknife `.R` is ignored as it is based on the N leave-one-out datasets.
#' 
#' Choosing `resample_method = "permutation"` for ungrouped data causes an error
#' as permutation will simply reorder the observations which is usually not 
#' meaningful. If a list of data is provided 
#' each list element is assumed to represent the observations belonging to one
#' group. In this case, data is pooled and group adherance permutated.
#' 
#' For cross-validation the number of folds (`k`) defaults to `10`. It may be
#' changed via the `.cv_folds` argument. Setting `k = N` (where `N` is the 
#' number of observations) produces leave-one-out cross-validation samples.
#' Note: 1. `k` can not be larger than the `N`. 2. If `N/k` is not not an integer 
#' the last fold will have less observations.
#' 
#' Random number generation (RNG) uses the L'Ecuyer-CRMR RGN stream as implemented
#' in the \href{https://github.com/HenrikBengtsson/future.apply}{future.apply package} 
#' \insertCite{Bengtsson2018a}{cSEM}.
#' See [?future_lapply][future.apply::future_lapply] for details. By default
#' a random seed is choosen.
#' 
#' @usage resampleData(
#'  .object          = NULL,
#'  .data            = NULL,
#'  .resample_method = c("bootstrap", "jackknife", "permutation", 
#'                       "cross-validation"),
#'  .cv_folds        = 10,
#'  .id              = NULL,
#'  .R               = 499,
#'  .seed            = NULL
#' )
#'
#' @param .data A `data.frame`, a `matrix` or a `list` of data of either type. 
#'   Possible column types or classes of the data provided are: 
#'   "`logical`", "`numeric`" ("`double`" or "`integer`"), "`factor`" (ordered and unordered) 
#'   or a mix of several types. The data may also include
#'   **one** character column whose column name must be given to `.id`. 
#'   This column is assumed to contain group identifiers used to split the data into groups.
#'   If `.data` is provided, `.object` is ignored. Defaults to `NULL`.
#' @param .resample_method Character string. The resampling method to use. One of: 
#'  "*bootstrap*", "*jackknife*", "*permutation*", or "*cross-validation*". 
#'  Defaults to "*bootstrap*".
#' @param .R Integer. The number of bootstrap runs, permutation runs
#'   or cross-validation repetitions to use. Defaults to `499`.
#' @inheritParams csem_arguments
#' 
#' @return The structure of the output depends on the type of input and the 
#'   resampling method:
#' \describe{
#' \item{Bootstrap}{If a `matrix` or `data.frame` without grouping variable 
#'   is provided (i.e., `.id = NULL`), the result is a list of length `.R` 
#'   (default `499`). Each element of that list is a bootstrap (re)sample.
#'   If a grouping variable is specified or a list of data is provided 
#'   (where each list element is assumed to contain data for one group), 
#'   resampling is done by group. Hence, 
#'   the result is a list of length equal to the number of groups 
#'   with each list element containing `.R` bootstrap samples based on the 
#'   `N_g` observations of group `g`.}
#' \item{Jackknife}{If a `matrix` or `data.frame` without grouping variable 
#'   is provided (`.id = NULL`), the result is a list of length equal to the number
#'   of observations/rows (`N`) of the dataset provided. 
#'   Each element of that list is a jackknife (re)sample.
#'   If a grouping variable is specified or a list of data is provided 
#'   (where each list element is assumed to contain data for one group), 
#'   resampling is done by group. Hence, 
#'   the result is a list of length equal to the number of group levels 
#'   with each list element containing `N` jackknife samples based on the 
#'   `N_g` observations of group `g`.}
#' \item{Permutation}{If a `matrix` or `data.frame` without grouping variable 
#'   is provided an error is returned as permutation will simply reorder the observations.
#'   If a grouping variable is specified or a list of data is provided 
#'   (where each list element is assumed to contain data of one group), 
#'   group membership is permutated. Hence, the result is a list of length `.R`
#'   where each element of that list is a permutation (re)sample.}
#' \item{Cross-validation}{If a `matrix` or `data.frame` without grouping variable 
#'   is provided a list of length `.R` is returned. Each list element
#'   contains a list containing the `k` splits/folds subsequently
#'   used as test and training datasets.  
#'   If a grouping variable is specified or a list of data is provided 
#'   (where each list element is assumed to contain data for one group), 
#'   cross-validation is repeated `.R` times for each group. Hence, 
#'   the result is a list of length equal to the number of groups,
#'   each containing `.R` list elements (the repetitions) which in turn contain 
#'   the `k` splits/folds.
#'   }
#' }
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [cSEMResults], [resamplecSEMResults()]
#'
#' @examples
#' # ===========================================================================
#' # Using the raw data 
#' # ===========================================================================
#' ### Bootstrap (default) -----------------------------------------------------
#' 
#' res_boot1 <- resampleData(.data = satisfaction)
#' str(res_boot1, max.level = 3, list.len = 3)
#' 
#' ## To replicate a bootstrap draw use .seed:
#' res_boot1a <- resampleData(.data = satisfaction, .seed = 2364)
#' res_boot1b <- resampleData(.data = satisfaction, .seed = 2364)
#'                            
#' identical(res_boot1, res_boot1a) # TRUE
#' 
#' ### Jackknife ---------------------------------------------------------------
#' 
#' res_jack <- resampleData(.data = satisfaction, .resample_method = "jackknife")
#' str(res_jack, max.level = 3, list.len = 3)
#' 
#' ### Cross-validation --------------------------------------------------------
#' ## Create dataset for illustration:
#' dat <- data.frame(
#'   "x1" = rnorm(100),
#'   "x2" = rnorm(100),
#'   "group" = sample(c("male", "female"), size = 100, replace = TRUE),
#'   stringsAsFactors = FALSE)
#' 
#' ## 10-fold cross-validation (repeated 100 times)
#' cv_10a <- resampleData(.data = dat, .resample_method = "cross-validation", 
#'                       .R = 100)
#' str(cv_10a, max.level = 3, list.len = 3)
#' 
#' # Cross-validation can be done by group if a group identifyer is provided:
#' cv_10 <- resampleData(.data = dat, .resample_method = "cross-validation", 
#'                       .id = "group", .R = 100)
#' 
#' ## Leave-one-out-cross-validation (repeated 50 times)
#' cv_loocv  <- resampleData(.data = dat[, -3], 
#'                           .resample_method = "cross-validation", 
#'                           .cv_folds = nrow(dat),
#'                           .R = 50)
#' str(cv_loocv, max.level = 2, list.len = 3)
#' 
#' ### Permuation ---------------------------------------------------------------
#' 
#' res_perm <- resampleData(.data = dat, .resample_method = "permutation",
#'                          .id = "group")
#' str(res_perm, max.level = 2, list.len = 3)
#' 
#' # Forgetting to set .id causes an error
#' \dontrun{
#' res_perm <- resampleData(.data = dat, .resample_method = "permutation")
#' }
#' 
#' # ===========================================================================
#' # Using a cSEMResults object
#' # ===========================================================================
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
#' EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
#' IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
#' LOY  =~ loy1  + loy2  + loy3  + loy4
#' QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
#' SAT  =~ sat1  + sat2  + sat3  + sat4
#' VAL  =~ val1  + val2  + val3  + val4
#' "
#' a <- csem(satisfaction, model)
#' 
#' # Create bootstrap and jackknife samples
#' res_boot <- resampleData(a, .resample_method = "bootstrap", .R = 999)
#' res_jack <- resampleData(a, .resample_method = "jackknife")
#' 
#' # Since `satisfaction` is the dataset used the following approaches yield
#' # identical results.
#' res_boot_data   <- resampleData(.data = satisfaction, .seed = 2364)
#' res_boot_object <- resampleData(a, .seed = 2364)
#' 
#' identical(res_boot_data, res_boot_object) # TRUE
#' 
#' @export
#'
resampleData <- function(
  .object          = NULL,
  .data            = NULL,
  .resample_method = c("bootstrap", "jackknife", "permutation", 
                       "cross-validation"),
  .cv_folds        = 10,
  .id              = NULL,
  .R               = 499,
  .seed            = NULL
) {
  .resample_method <- match.arg(.resample_method, 
            c("bootstrap", "jackknife", "permutation", "cross-validation"))

  ## Set plan on how to resolve futures to "sequential" as it is virtually always
  ## slower to resample data using "multiprocess".; reset at the end
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan("sequential")
  
  ## Get data set
  if(is.null(.data)) {
    ## Get information according to class of object
    if(any(class(.object) %in% "cSEMResults_default")) {
      data <- as.data.frame(.object$Information$Arguments$.data)
      id   <- NULL
    } else if(any(class(.object) %in% "cSEMResults_multi")) {
      data        <- .object[[1]]$Information$Data_pooled
      data_split  <- lapply(.object, function(x) x$Information$Arguments$.data)
      id          <- ifelse(is.null(.object[[1]]$Information$Arguments$.id), 
                            "id", .object[[1]]$Information$Arguments$.id)
    } else if(any(class(.object) %in% "cSEMResults_2ndorder")) {
      data <- .object$First_stage$Information$Arguments$.data
      id   <- NULL
    } else {
      stop2("The following error occured in the `resampleData()` function:\n",
            "`object` must be of class cSEMResults."
            )
    }
  } else {
    
    ## Check data
    if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
      stop2(
        "The following error occured in the `resampleData()` function:\n",
        "Data must be provided as a `matrix`, a `data.frame` or a `list`. ", 
         ".data has class: ", class(.data)
        )
    }
    ## Select cases
    if(!is.null(.id) && !inherits(.data, "list")) {
      
      if(length(.id) != 1) {
        stop2(
          "The following error occured in the `resampleData()` function:\n",
          "`.id` must be a character string or an integer identifying one single column."
          )
      }
      
      if(is.matrix(.data)) {
        .data <- as.data.frame(.data)
      }
      
      data       <- .data 
      data_split <- split(.data, f = .data[, .id])
      id         <- .id

    } else if(any(class(.data) == "list")) {
      
      data <- do.call(rbind, .data)

      if(is.null(names(.data))) {
        data[, "id"] <- rep(paste0("Data_", 1:length(.data)), 
                            times = sapply(.data, nrow))
      } else {
        data[, "id"] <- rep(names(.data), times = sapply(.data, nrow))
      }
      data_split <- .data
      id <- "id"
      
      ##add names here
      
    } else {
      data <- .data
      id   <- NULL
    }
  }
  
  ## Create seed if none is given
  if(is.null(.seed)) {
    .seed <- sample(.Random.seed, 1)
  }
  
  ## Choose resampling method
  out <- switch (.resample_method,
    "jackknife"   = {
      if(exists("data_split")) {
        future.apply::future_lapply(data_split, function(y) future.apply::future_lapply(1:nrow(y), function(x) y[-x, ]))
      } else {
        future.apply::future_lapply(1:nrow(data), function(x) data[-x, ]) 
      }
    },
    "bootstrap"   = {
      if(exists("data_split")) {
        future.apply::future_lapply(data_split, function(y) {
          future.apply::future_lapply(1:.R, function(x) {
            y[sample(1:nrow(y), size = nrow(y), replace = TRUE), ]
          }, future.seed = .seed)
        })
      } else {
        future.apply::future_lapply(1:.R, function(x) {
          data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]}, 
          future.seed = .seed
        )
      }
    },
    "permutation" = {
      if(is.null(id)) {
        stop2(
          "The following error occured in the `resampleData()` function:\n",
          "No id column specified to permutate the data with."
        )
      } else {
        future.apply::future_lapply(1:.R, function(x) {
          cbind(data[,-which(colnames(data) == id)], "id" = sample(data[, id]))},
          future.seed = .seed)
      }
    },
    "cross-validation" = {
      # k-fold cross-validation (=draw k samples of equal size.).
      # Note the last sample may contain less observations if equal sized
      # samples are not possible
      if(exists("data_split")) {
        ## The number of folds cannot be larger than the minimum number of
        ## rows/observations for all data sets in data_split
        if(max(sapply(data_split, nrow)) < .cv_folds) {
          stop2(
            "The following error occured in the `resampleData()` function:\n",
            "The number of folds is larger than the number of observations", 
            " in at least one of the groups."
          )
        }
        future.apply::future_lapply(data_split, function(y) {
          future.apply::future_lapply(1:.R, function(x) {
            # shuffle data set
            y <- y[sample(1:nrow(y)), ]
            suppressWarnings(
              split(as.data.frame(y), rep(1:.cv_folds, 
                                          each = ceiling(nrow(y)/.cv_folds)))
            )
          }, future.seed = .seed)
        })
      } else {
        ## The number of folds cannot be larger than the minimum number of
        ## rows/observations for all data sets in data_split
        if(nrow(data) < .cv_folds) {
          stop2(
            "The following error occured in the `resampleData()` function:\n",
            "The number of folds is larger than the number of observations."
          )
        }
        # shuffle data
        future.apply::future_lapply(1:.R, function(x) {
          data <- data[sample(1:nrow(data)), ]
          suppressWarnings(
            split(as.data.frame(data), rep(1:.cv_folds, 
                                           each = ceiling(nrow(data)/.cv_folds)))
          )
        }, future.seed = .seed)
      }
    } # END cross-validation
  ) # END switch
  ## Return samples
  out
}

#' Resample cSEMResults 
#' 
#' Resample a [cSEMResults] object using bootstrap or jackknife resampling. 
#' The function is called by [csem()] if the user sets 
#' `csem(..., .resample_method = "bootstrap")` or 
#' `csem(..., .resample_method = "jackknife")` but may also be called directly.
#' 
#' Given `M` resamples (for bootstrap `M = .R` and for jackknife `M = N`, where
#' `N` is the number of observations) based on the data used to compute the
#' [cSEMResults] object provided via `.object`, `resamplecSEMResults()` essentially calls 
#' [csem()] on each resample using the arguments of the origianl call (ignoring any arguments
#' related to resampling) and returns estimates for each of a subset of 
#' practically useful resampled parameters/statistics computed by [csem()]. 
#' Currently, the following estimates are computed and returned by default based 
#' on each resample: Path estimates, Loading estimates, Weight estimates.
#' 
#' In practical application users may need to resample a specific statistic (e.g,
#' the heterotrait-monotrait ratio of correlations (HTMT) or differences between path 
#' coefficients such as beta_1 - beta_2).
#' Such statistics may be provided by a function `fun(.object, ...)` or a list of 
#' such functions via the `.user_funs` argument. The first argument of 
#' these functions must always be `.object`. 
#' Internally, the function will be applied on each  
#' resample to produce the desired statistic. Hence, arbitrary complicated statistics
#' may be resampled as long as the body of the function draws on elements contained
#' in the [cSEMResults] object only. Output of `fun(.object, ...)` should preferably 
#' be a (named) vector but matrices are also accepted. 
#' However, the output will be vectorized (columnwise) in this case. 
#' See the examples section for details.
#'
#' Both resampling the origianl [cSEMResults] object (call it "first resample") 
#' and resampling based on a resampled [cSEMResults] object (call it "second resample") 
#' are supported. Choices for the former 
#' are "*bootstrap*" and "*jackknife*". Resampling based on a resample is turned off
#' by default (`.resample_method2 = "none"`) as this significantly
#' increases computation time (there are now `M * M2` resamples to compute, where
#' `M2` is `.R2` or `N`).
#' Resamples of a resample are required, e.g., for the studentized confidence
#' interval computed by the [infer()] function. Typically, bootstrap resamples
#' are used in this case \insertCite{Davison1997}{cSEM}.
#' 
#' As [csem()] accepts a single dataset, a list of datasets as well as datasets
#' that contain a column name used to split the data into groups,
#' the [cSEMResults] object may contain multiple datasets. 
#' In this case, resampling is done by dataset or group. Note that depending
#' on the number of datasets/groups, the computation may be considerably
#' slower as resampling will be repeated for each dataset/group. However, apart
#' from speed considerations users don not need to worry about the type of
#' input used to compute the [cSEMResults] object as `resamplecSEMResults()`
#' is able to deal with each case.
#' 
#' The number of bootstrap runs for the first and second run are given by `.R` and `.R2`. 
#' The default is `499` for the first and `199` for the second run 
#' but should be increased in real applications. See e.g.,
#' \insertCite{Hesterberg2015;textual}{cSEM}, p.380, 
#' \insertCite{Davison1997;textual}{cSEM}, and
#' \insertCite{Efron2016}{cSEM} for recommendations.
#' For jackknife `.R` are `.R2` are ignored. 
#' 
#' Resampling may produce inadmissble results (as checked by [verify()]).
#' By default these results are dropped however users may choose to `"ignore"`
#' or `"replace"` inadmissble results in which resampling continious until
#' the necessary number of admissble results is reached.
#' 
#' The \pkg{cSEM} package supports (multi)processing via the \href{https://github.com/HenrikBengtsson/future}{future} 
#' framework \insertCite{Bengtsson2018}{cSEM}. Users may simply choose an evaluation plan
#' via `.eval_plan` and the package takes care of all the complicated backend 
#' issues. Currently, users may choose between standard single-core/single-session
#'  evaluation (`"sequential"`) and multiprocessing (`"multiprocess"`). The future package
#' provides other options (e.g., `"cluster"` or `"remote"`), however, they probably 
#' will not be needed in the context of the \pkg{cSEM} package as simulations usually
#' do not require high-performance clusters. Depeding on the operating system, the future
#' package will manage to distribute tasks to multiple R sessions (Windows)
#' or multiple cores. Note that multiprocessing is not necessary always faster
#' when only a "small" number of replications is required as the overhead of
#' initializing new sessions or distributing tasks to different cores 
#' will not immediatley be compensated by the avaiability of multiple sessions/cores.
#'
#' Random number generation (RNG) uses the L'Ecuyer-CRMR RGN stream as implemented in the
#' \href{https://github.com/HenrikBengtsson/future.apply}{future.apply package} \insertCite{Bengtsson2018a}{cSEM}.
#' It is independent of the evaluation plan. Hence, setting e.g., `.seed = 123` will
#' generate the same random number and replicates
#' for both `.eval_plan = "sequential"` and `.eval_plan = "multiprocess"`.
#' See [?future_lapply][future.apply::future_lapply] for details.
#' 
#' @usage resamplecSEMResults(
#'  .object                = NULL,
#'  .resample_method       = c("bootstrap", "jackknife"), 
#'  .resample_method2      = c("none", "bootstrap", "jackknife"), 
#'  .R                     = 499,
#'  .R2                    = 199,
#'  .handle_inadmissibles  = c("drop", "ignore", "replace"),
#'  .user_funs             = NULL,
#'  .eval_plan             = c("sequential", "multiprocess"),
#'  .seed                  = NULL,
#'  .sign_change_option    = c("none","individual","individual_reestimate",
#'                             "construct_reestimate"),
#'  ...
#' ) 
#'
#' @inheritParams csem_arguments
#' @param .resample_method Character string. The resampling method to use. One of: 
#'  "*bootstrap*" or "*jackknife*". Defaults to "*bootstrap*".
#' @param ... Further arguments passed to functions supplied to `.user_funs`.
#'   
#' @return The core structure is the same structure as that of `.object` with
#' the following elements added:
#' \itemize{
#' \item{ `$Estimates_resamples`: A list containing the `.R` resamples and
#' the original estimates for each of the resampled quantities (Path_estimates, 
#' Loading_estimates, Weight_estimates, user defined functions). 
#' Each list element is a list containing elements 
#' `$Resamples` and `$Original`. `$Resamples` is a `(.R x K)` matrix with each
#' row representing one resample for each of the `K` parameters/statistics.
#' `$Original` contains the original estimates (vectorized by column if the output of 
#' the user provided function is a matrix.}
#' \item {`$Information_resamples`: A list containing addtional information.}
#' }
#' Use `str(<.object>, list.len = 3)` on the resulting object for an overview.
#' 
#' @references
#'   \insertAllCited{} 
#'   
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \donttest{
#' model <- "
#' # Structural model
#' QUAL ~ EXPE
#' EXPE ~ IMAG
#' SAT  ~ IMAG + EXPE + QUAL + VAL
#' LOY  ~ IMAG + SAT
#' VAL  ~ EXPE + QUAL
#' 
#' # Measurement model
#' EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
#' IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
#' LOY  =~ loy1  + loy2  + loy3  + loy4
#' QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
#' SAT  =~ sat1  + sat2  + sat3  + sat4
#' VAL  =~ val1  + val2  + val3  + val4
#' "
#' 
#' a <- csem(satisfaction, model)
#' 
#' # Bootstrap and jackknife estimation
#' boot <- resamplecSEMResults(a)
#' jack <- resamplecSEMResults(a, .resample_method = "jackknife") 
#' 
#' # Bootstrap resampling with resampling each resample (totals .R*R2 runs).
#' boot1 <- resamplecSEMResults(a, .resample_method = "bootstrap", .R = 50,
#'                              .resample_method2 = "bootstrap", .R2 = 20,
#'                              .seed = 1303)
#'
#' # To get inferencial quanitites such as the estimated standard error or
#' # the percentile confidence intervall for each resampled quantity use infer()
#' inference <- infer(boot1)
#' inference$Path_estimates$Standard_deviation
#' inference$Path_estimates$CI_percentile
#' 
#' # As usual summarize() can be called directly
#' summarize(boot1)
#' 
#' # In the example above .R x .R2 = 50 x 20 = 1000. Multiprocessing will be
#' # faster here and is therefore recommended. Note that multiprocessing
#' # does not affect the random number generation
#' boot2 <- resamplecSEMResults(a, .resample_method = "bootstrap", .R = 50,
#'                              .resample_method2 = "bootstrap", .R2 = 20,
#'                              .eval_plan = "multiprocess", .seed = 1303)
#'
#' identical(boot1, boot2)
#' }                  
#' @export
#'
resamplecSEMResults <- function(
  .object                = NULL, 
  .resample_method       = c("bootstrap", "jackknife"),
  .resample_method2      = c("none", "bootstrap", "jackknife"),
  .R                     = 499,
  .R2                    = 199,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .user_funs             = NULL,
  .eval_plan             = c("sequential", "multiprocess"),
  .seed                  = NULL,
  .sign_change_option    = c("none","individual","individual_reestimate",
                             "construct_reestimate"),
  ...
  ) {
  
  ## Does .object already contain resamples
  if(any(class(.object) == "cSEMResults_resampled")) {
    stop2("The following issue was encountered in the `resamplecSEMResults()` function:\n",
          "The object provided already contains resamples.")
  }
  
  ## Check for the minimum number of necessary resamples
  if(.R < 3 | .R2 < 3) {
    stop2("The following error occured in the `resamplecSEMResults()` function:\n",
          "At least 3 resamples required.")
  }
  
  ## Has the object to use the data to resample from produced admissible results?
  if(sum(unlist(verify(.object))) != 0) {
    warning2(
      "The following issue was encountered in the `resamplecSEMResults()` function:\n",
      "Estimation based on the original data has produced inadmissible results.\n", 
      "This may be a sign that something is wrong.",
      " Resampling will continue but may not produce reliable results.")
  }
  
  ## Match arguments
  .resample_method  <- match.arg(.resample_method)
  .resample_method2 <- match.arg(.resample_method2)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .eval_plan <- match.arg(.eval_plan)
  .sign_change_option <- match.arg(.sign_change_option)
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, function(x) {
      resamplecSEMResults(
        .object               = x,
        .resample_method      = .resample_method,
        .resample_method2     = .resample_method2,
        .R                    = .R,
        .R2                   = .R2,
        .handle_inadmissibles = .handle_inadmissibles,
        .user_funs            = .user_funs,
        .eval_plan            = .eval_plan,
        .seed                 = .seed,
        .sign_change_option   = .sign_change_option,
        ...
      )
    })
    
    ## Add/ set class
    class(out) <- c(class(.object), "cSEMResults_resampled")
    return(out)
  }
  
  # Set plan on how to resolve futures; reset at the end
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(.eval_plan)
  
  ### Process original data ----------------------------------------------------
  ## Select relevant statistics/parameters/quantities and vectorize 
  Est_original <- selectAndVectorize(.object)
  
  ## Apply user defined function if specified
  user_funs <- if(!is.null(.user_funs)) {
    applyUserFuns(.object, .user_funs = .user_funs, ...)
  }
  
  ## Add output of the user functions to Est_original
  if(!is.null(.user_funs)) {
    Est_original <- c(Est_original, user_funs)
  }
  
  ### Resample and compute -----------------------------------------------------
  # Create seed if none is given
  if(is.null(.seed)) {
    .seed <- sample(.Random.seed, 1)
  }
  
  out <- resamplecSEMResultsCore(
    .object                = .object, 
    .resample_method       = .resample_method,
    .resample_method2      = .resample_method2,
    .R                     = .R,
    .R2                    = .R2,
    .handle_inadmissibles  = .handle_inadmissibles,
    .handle_inadmissibles2 = .handle_inadmissibles, 
    .user_funs             = .user_funs,
    .eval_plan             = .eval_plan,
    .seed                  = .seed,
    .sign_change_option    = .sign_change_option,
    ...
  )
  
  # Check if at least 3 admissible results were obtained
  n_admissibles <- length(out)
  if(n_admissibles < 3) {
    stop("The following error occured in the `resamplecSEMResults()` functions:\n",
         "Less than 2 admissible results produced.", 
         " Consider setting `.handle_inadmissibles == 'replace'` instead.",
         call. = FALSE)
  }
  
  # Turn list "inside out" and bind bootstrap samples to matrix 
  # - columns are variables
  # - rows are bootstrap runs
  out <- purrr::transpose(out) 
  
  if(.resample_method2 != "none") {
    out_2 <- out$Estimates2
    out   <- purrr::transpose(out$Estimates1)
  }
  
  out <- out %>% 
    lapply(function(x) do.call(rbind, x))
  
  # Add estimated quantities based on the the original sample/data
  out <- mapply(function(l1, l2) list("Resampled" = l1, "Original" = l2),
                l1 = out,
                l2 = Est_original,
                SIMPLIFY = FALSE)
  
  ## Return --------------------------------------------------------------------
  # The desired output is: a list of two:
  #  1. Resamples := a list of two
  #      1.1 Estimates1 := the .R resamples of the "outer" resampling run.
  #      1.2 Estimates2 := the .R2 resamples of the "inner" run or NA
  #  2. Information := a list of three
  #      2.1 Method  := the method used to obtain the "outer" resamples
  #      2.2 Method2 := the method used to obtain the "inner" resamples or NA
  #      2.3 Number_of_observations := the number of observations.
  #      2.4 Original_object        := the original cSEMResults object
  #
  # Since resamplecSEMResults is called recursivly within its body it is 
  # difficult to produce the desired output. There is now straightforward way 
  # to determine if a call is a recursive call or not
  # My workaround for now is to simply check if the argument provided for ".object"
  # is "Est_temp" since this is the argument for the recurive call. This
  # is of course not general but a dirty workaround. I guess it is save enough
  # though.
  is_recursive_call <- eval.parent(as.list(match.call()))$.object == "Est_temp"
  
  ## Return
  if(is_recursive_call) {
    out
  } else {
    if(.resample_method2 != "none") {
      out <- list("Estimates1" = out, "Estimates2" = out_2)
      
    } else {
      out <- list("Estimates1" = out, "Estimates2" = NA)
    }
    
    ## Add resamples and additional information 
    if(inherits(.object, "cSEMResults_2ndorder")) {
      resample_out <- c(.object$Second_stage$Information)
      resample_out[[length(resample_out) + 1]] <- list(
        "Estimates" = c(out),
        "Information_resample" =       list(
          "Method"                  = .resample_method,
          "Method2"                 = .resample_method2,
          "Number_of_admissibles"   = n_admissibles,
          "Number_of_observations"  = nrow(.object$Second_stage$Information$Data),
          "Number_of_runs"          = .R,
          "Number_of_runs2"         = .R2,
          "Seed"                    = .seed,
          "Handle_inadmissibles"    = .handle_inadmissibles,
          "Sign_change_option"      = .sign_change_option
        )
      )
      names(resample_out)[length(resample_out)] <- "Resamples"
      
      out <- list(
        "First_stage"   = .object$First_stage,
        "Second_stage"  = list(
          "Estimates" = .object$Second_stage$Estimates,
          "Information" = resample_out
        )
      )
      
      # Renew class for 2nd stage
      class(out$Second_stage) <- class(.object$Second_stage) 
      
    } else {
      # Estimates
      estim <- c(.object$Estimates)
      estim[[length(estim) + 1]] <- c(out)
      names(estim)[length(estim)] <- "Estimates_resample"
      
      # Information
      info <- c(.object$Information)
      info[[length(info) + 1]] <- list(
        "Method"                  = .resample_method,
        "Method2"                 = .resample_method2,
        "Number_of_admissibles"   = n_admissibles,
        "Number_of_observations"  = nrow(.object$Information$Data), # for infer()
        "Number_of_runs"          = .R,
        "Number_of_runs2"         = .R2,
        "Seed"                    = .seed,
        "Handle_inadmissibles"    = .handle_inadmissibles,
        "Sign_change_option"      = .sign_change_option
      )
      names(info)[length(info)] <- "Information_resample"
      
      out <- list(
        "Estimates"   = estim,
        "Information" = info
      )
    }
    
    ## Add/ set class
    class(out) <- c(class(.object), "cSEMResults_resampled")
    return(out)
  }
}

#' Core tasks of the resamplecSEMResults function
#' @noRd
#' 
resamplecSEMResultsCore <- function(
  .object                = args_default()$.object, 
  .resample_method       = args_default()$.resample_method,
  .resample_method2      = args_default()$.resample_method2,
  .R                     = args_default()$.R,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .handle_inadmissibles2 = NULL,
  .user_funs             = args_default()$.user_funs,
  .eval_plan             = args_default()$.eval_plan,
  .seed                  = args_default()$.seed,
  .sign_change_option    = args_default()$.sign_change_option,
  ...
) {
  
  ## Get arguments
  if(inherits(.object, "cSEMResults_2ndorder")) {
    info1 <- .object$First_stage$Information
    info2 <- .object$Second_stage$Information
    
    est1_normal <- .object$First_stage$Estimates
    est2_normal <- .object$Second_stage$Estimates
    
    summary_original <- summarize(.object)
    est1 <- summary_original$First_stage$Estimates
    est2 <- summary_original$Second_stage$Estimates
    
    args <- info2$Arguments_original

    ## Append the 2ndorder approach to args
    args$.approach_2ndorder <- info2$Approach_2ndorder
    
    ## Replace original data (which contains either an id column or is a list)
    ## by the data used for the group that is being resampled.
    if(!is.null(args$.id) | inherits(args$.data, "list")) {
      args$.data <- info1$Data
      args$.id <- NULL
    }
  } else {
    info1       <- .object$Information
    est_normal  <- .object$Estimates
    
    summary_original <- summarize(.object)
    est <- summary_original$Estimates
    
    args <- info1$Arguments 
  }
  
  ## Resample jackknife
  if(.resample_method == "jackknife") {
    resample_jack <- resampleData(.object, .resample_method = "jackknife") 
    .R <- length(resample_jack)
  } 
  
  ### Start resampling loop ====================================================
  Est_ls <- future.apply::future_lapply(1:.R, function(i) {
  # Est_ls <- lapply(1:.R, function(i) {
    # Replace the old dataset by a resampled data set (resampleData always returns
    # a list so for just one draw we need to pick the first list element)
    
    data_temp <- if(.resample_method == "jackknife") {
      resample_jack[[i]]
    } else {
      # We could use resampleData here but, bootstrap resampling directly is faster
      # (not surprising)
      # (compared both approaches using microbenchmark)
      data <- args[[".data"]]
      data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]
    }
    
    args[[".data"]] <- data_temp
    
    # Estimate model
    Est_temp <- if(inherits(.object, "cSEMResults_2ndorder")) {
      ## NOTE: using do.call(csem, args) would be more elegant but is much 
      # much much! slower (especially for larger data sets). 
      do.call(csem, args) 

    } else {
      # It is important to use foreman() here 
      # instead of csem() to allow for lapply(x, resamplecSEMResults_default) when x 
      # is of class cSEMResults_2ndorder.
      
      ## NOTE: 
      #  01.03.2019: Using do.call(foreman, args) would be more elegant but is much 
      #              much much! slower (especially for larger data sets).
      #
      #  15.05.2019: Apparently do.call(foreman, args) is not that bad
      #              after all. I did several comparisons using microbenchmark
      #              but there was no speed difference (anymore?!). When I compared and
      #              thought that do.call is slow I fixed other parts of 
      #              foreman as well...maybe they had been the real culprit.
      #              So we use do.call again, as it avoids redundancy
      do.call(foreman, args)
    }
    
    # Check status
    status_code <- sum(unlist(verify(Est_temp)))
    
    # Distinguish depending on how inadmissibles should be handled
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
      
      
      # ## Select relevant statistics/parameters/quantities
      x1 <- selectAndVectorize(Est_temp)
      
      ### Sign change correction for PLS-PM
      if(.sign_change_option != "none" && info1$Arguments$.approach_weights == "PLS-PM") {
        
        if(inherits(.object, "cSEMResults_2ndorder")) {
          
          est1_temp_normal <- Est_temp$First_stage$Estimates
          est2_temp_normal <- Est_temp$Second_stage$Estimates
        
          ### Individual_reestimate and construct_reestimate -------------------
          if(.sign_change_option == "individual_reestimate" | 
             .sign_change_option == "construct_reestimate") {
            
            ## First stage
            # Is there a sign difference between the first stage resample 
            # and original weight estimates. If so, which weights differ?
            sign_diff1 <- sign(est1_normal$Weight_estimates) != 
              sign(est1_temp_normal$Weight_estimates)
            
            # Only if at least one sign differs a correction is needed
            if(sum(sign_diff1) != 0) {
              
              W1_new <- est1_temp_normal$Weight_estimates
              
              # Individual_reestimate
              if(.sign_change_option == "individual_reestimate"){
                
                W1_new[sign_diff1] <- W1_new[sign_diff1] * (-1)
              } # END individual_reestimate
              
              if(.sign_change_option == "construct_reestimate"){
                
                # In line with Tenenhaus et al. (2005, p 177) loadings are used,
                # although they have not considered 2nd order models
                
                Lambda1 <- est1_normal$Loading_estimates
                Lambda1_temp <- est1_temp_normal$Loading_estimates
                
                Lambda_diff1 <- abs(rowSums(Lambda1 - Lambda1_temp))
                Lambda_sum1 <- abs(rowSums(Lambda1 + Lambda1_temp))
                
                W1_new[Lambda_diff1 > Lambda_sum1, ] <- W1_new[Lambda_diff1 > Lambda_sum1, ] * (-1)
              } # END construct_reestimate
              
              # Put corrected weights in a list to be able to supply them as
              # fixed weights to .PLS_modes
              W1_list <- lapply(1:nrow(W1_new), function(x) {
                temp <- W1_new[x, ]
                temp[temp !=0]
              })
              names(W1_list) <- rownames(W1_new)
              
              ## Replace modes by fixed sign-corrected weights
              args[[".PLS_modes"]] <- W1_list
              
              ## Reestimate
              # Note: calling csem() directly is faster, however, i guess
              # change option wont be used that often, so for now, I will
              # use the more elegant, while slower, do.call construction.
              Est_temp <- do.call(csem, args)
              
              ## Update
              est2_temp_normal <- Est_temp$Second_stage$Estimates
              
            } # END first stage
            
            ## Second stage
            # Is there a sign difference between the second stage (sign-corrected) 
            # resample weight estimates and the original weight estimates. 
            # If so, which weights differ?
            sign_diff2 <- sign(est2_normal$Weight_estimates) != 
              sign(est2_temp_normal$Weight_estimates)
  
            # Only if at least one sign differs a correction is needed
            if(sum(sign_diff2) != 0) {
              
              W2_new <- est2_temp_normal$Weight_estimates
            
              # Individual_reestimate
              if(.sign_change_option == "individual_reestimate"){
                
                W2_new[sign_diff2] <- W2_new[sign_diff2] * (-1)
              } # END individual_reestimate
              
              if(.sign_change_option == "construct_reestimate"){
                
                # In line with Tenenhaus et al. (2005, p 177) loadings are used,
                # although they have not considered 2nd order models
                
                Lambda2 <- est2_normal$Loading_estimates
                Lambda2_temp <- est2_temp_normal$Loading_estimates
                
                Lambda_diff2 <- abs(rowSums(Lambda2 - Lambda2_temp))
                Lambda_sum2 <- abs(rowSums(Lambda2 + Lambda2_temp))
                
                W2_new[Lambda_diff2 > Lambda_sum2, ] <- W2_new[Lambda_diff2 > Lambda_sum2, ] * (-1)
              } # END construct_reestimate
              
              # Put corrected weights in a list to be able to supply them as
              # fixed weights to .PLS_modes
              W2_list <- lapply(1:nrow(W2_new), function(x){
                temp <- W2_new[x, ]
                temp[temp != 0]
              })
              names(W2_list) <- rownames(W2_new)
              
              args[[".PLS_modes"]] = W2_list
              
              ## Reestimate
              # Note: calling csem() directly is faster, however, i guess
              # change option wont be used that often, so for now, I will
              # use the more elegant, while slower, do.call construction.
              Est_temp <- do.call(csem, args)
              
            } # END second stage
            
            # ## Update using estimates based on the sign-corrected weights
            x1 <- selectAndVectorize(Est_temp)
          } # END individual_reestimate, construct_reestimate
          
          if(.sign_change_option == "individual") {
            
            # Reverse the signs off ALL parameter estimates in a bootstrap run if 
            # their sign differs from the sign of the original estimation
            x1 <- reverseSign(.Est_temp = Est_temp, .summary_original = summary_original)
            
          }
      } else {
        est_temp_normal <- Est_temp$Estimates
        
        sign_diff <- sign(est_normal$Weight_estimates) != 
          sign(est_temp_normal$Weight_estimates)
        
        # Is there a difference in the signs of th weights? Otherwise no correction of the signs is done
        # Not sure whether this is a problem for the construct_reestimate approach which only compares the loadings
        # I think not.
        if(sum(sign_diff) !=0 ) {
          
          
          # Sign change correction individual_reestimate and construct_reestimate
          if(.sign_change_option == "individual_reestimate" | 
             .sign_change_option == "construct_reestimate") {
            
            W_new <- est_temp_normal$Weight_estimates
            
            # Individual_reestimate
            if(.sign_change_option == "individual_reestimate"){
              
              W_new[sign_diff] <- W_new[sign_diff] * (-1)
            } # END individual_reestimate
            
            if(.sign_change_option == "construct_reestimate"){
              
              Lambda <- est_normal$Loading_estimates
              Lambda_temp <- est_temp_normal$Loading_estimates
              
              Lambda_diff <- abs(rowSums(Lambda - Lambda_temp))
              Lambda_sum <- abs(rowSums(Lambda + Lambda_temp))
              
              W_new[Lambda_diff > Lambda_sum, ] <- W_new[Lambda_diff > Lambda_sum, ] * (-1)
            }
            
            # Put corrected weights in a list to be able to supply them as
            # fixed weights to .PLS_modes
            W_list <- lapply(1:nrow(W_new), function(x){
              temp <- W_new[x, ]
              temp[temp != 0]
            })
            names(W_list) <- rownames(W_new)
            
            args[[".PLS_modes"]] = W_list
            
            ## Reestimate
            # Note: calling csem() directly is faster, however, i guess
            # change option wont be used that often, so for now, I will
            # use the more elegant, while slower, do.call construction.
            Est_temp <- do.call(csem, args)
            
            ## Update using estimates based on the sign-corrected weights
            x1 <- selectAndVectorize(Est_temp)
            
          } # end if individual_reestimate, construct_reestimate
          
          if(.sign_change_option == 'individual') {

            # Reverse the signs off ALL parameter estimates in a bootstrap run if 
            # their sign differs from the sign of the original estimation
            x1 <- reverseSign(.Est_temp = Est_temp, .summary_original = summary_original)
            
          } # END .sign_change_option == "individual" 
        } # END sum(sign(.object$Estimates$Weight_estimates)!=sign(Est_temp$Estimates$Weight_estimates)
      } # END else (= not 2ndorder)
    } # END .sign_change_option
    
      ## Apply user defined function if specified
      user_funs <- if(!is.null(.user_funs)) {
        applyUserFuns(Est_temp, .user_funs = .user_funs, ...)
      }
      
      # user_funs <- if(!is.null(.user_funs)) {
      #   if(is.function(.user_funs)) {
      #     list("User_fun" = c(.user_funs(Est_temp, ...)))
      #   } else {
      #     x <- lapply(.user_funs, function(f) c(f(Est_temp, ...)))
      #     if(is.null(names(x))) {
      #       names(x) <- paste0("User_fun", 1:length(x))
      #     }
      #     x
      #   }
      # }
      
      ## Add output of the user functions to x1
      if(!is.null(.user_funs)) {
        x1 <- c(x1, user_funs)
      }
      
      ## Resampling from a bootstrap sample is required for the
      ## bootstraped t-interval CI (studentized CI), hence the second run
      # In the second run no sign change option is used. We can think about 
      # applying the same correction as in the first run
      if(.resample_method2 != "none") {
        
        Est_resamples2 <- resamplecSEMResults(
          .object               = Est_temp,
          .R                    = .R2,
          .handle_inadmissibles = .handle_inadmissibles2,
          .resample_method      = .resample_method2,
          .resample_method2     = "none",
          .user_funs            = .user_funs,
          .seed                 = .seed, 
          .sign_change_option   = "none",
          ...
        )
        x1 <- list("Estimates1" = x1, "Estimates2" = Est_resamples2)
      }
    } else if(status_code != 0 & .handle_inadmissibles != "ignore") {
      # Retrun NA if status is not okay and .handle_inadmissibles == "drop" or
      # "replace"
      x1 <- NA
    }
    ## Return
    x1
  }, future.seed = .seed)
  # })
          
  ## Process data --------------------------------------------------------------
  # Delete potential NA's
  out <- Filter(Negate(anyNA), Est_ls)
  
  ## Replace inadmissibles
  if(length(out) != .R && .resample_method == "bootstrap" && 
     .handle_inadmissibles == "replace") {
    
    n_attempt <- 0
    while (length(out) < .R) {
      n_attempt <- n_attempt + 1
      R_new <- .R - length(out)
      if(is.null(.seed)) {
        .seed <- sample(.Random.seed, 1)
      } else {
        .seed <- .seed + 1
      }
      
      Est_replace <- resamplecSEMResultsCore(
        .object               = .object,
        .R                    = R_new,
        .handle_inadmissibles = "drop",
        .handle_inadmissibles2= .handle_inadmissibles, 
        .resample_method      = .resample_method,
        .resample_method2     = .resample_method2,
        .R2                   = .R2,
        .user_funs            = .user_funs,
        .eval_plan            = .eval_plan,
        .seed                 = .seed,
        .sign_change_option   = .sign_change_option,
        ...
      )
      
      out <- c(out, Est_replace)
      
      if(n_attempt == 10*.R) {
        warning2(
          "The following issue was encountered in the `resamplecSEMResults()` function:\n",
          "Attempting to replace inadmissibles failed after ", 10*.R, " attempts.",
          " Returning all admissible results up to the last attempt."
          )
        break
      }
    }
  }
  
  out
}
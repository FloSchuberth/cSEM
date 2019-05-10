#' Data: satisfaction
#'
#' This dataset contains the variables from a customer satisfaction study of 
#' a Spanish credit institution on 250 customers. The data is identical to
#' the dataset provided by the \href{https://github.com/gastonstat/plspm}{plspm} package
#' but with the last column  (`gender`) removed. If you are looking for the original
#' dataset use the [satisfaction_gender] dataset.
#'
#' @docType data
#'
#' @description A data frame with 250 observations and 27 variables. 
#' Variables from 1 to 27 refer to six latent concepts: \code{IMAG}=Image, 
#' \code{EXPE}=Expectations, \code{QUAL}=Quality, \code{VAL}=Value, 
#' \code{SAT}=Satisfaction, and \code{LOY}=Loyalty.
#' \describe{
#'   \item{imag1-imag5}{Indicators attached to concept `IMAG` which is supposed to
#'                      capture aspects such as the institutions reputation, 
#'                      trustworthiness, seriousness, solidness, and caring 
#'                      about customer.}
#'   \item{expe1-expe5}{Indicators attached to concept `EXPE` which is supposed to
#'                      capture aspects concerning products and 
#'                      services provided, customer service, providing solutions,
#'                      and expectations for the overall quality.}
#'   \item{qual1-qual5}{Indicators attached to concept `QUAL` which is supposed to
#'                      capture aspects concerning reliability of products and services, 
#'                      the range of products and services, personal advice, 
#'                      and overall perceived quality.}
#'   \item{val1-val4}{Indicators attached to concept `VAL` which is supposed to
#'                      capture aspects related to beneficial services and 
#'                      products, valuable investments, quality relative to 
#'                      price, and price relative to quality.}
#'   \item{sat1-sat4}{Indicators attached to concept `SAT` which is supposed to
#'                      capture aspects concerning overall rating of satisfaction, 
#'                      fulfillment of expectations, satisfaction relative to 
#'                      other banks, and performance relative to customer's 
#'                      ideal bank.}
#'   \item{loy1-loy4}{Indicators attached to concept `LOY` which is supposed to
#'                    capture aspects concerning propensity to choose the 
#'                    same bank again, propensity to switch to other bank, 
#'                    intention to recommend the bank to friends, 
#'                    and the sense of loyalty.}
#' }
#' 
#' 
#' 
#' @source The \href{https://github.com/gastonstat/plspm}{plspm} package (version  0.4.9). 
#' Original source according to \pkg{plspm}:
#' "Laboratory of Information Analysis and Modeling (LIAM). 
#' Facultat d'Informatica de Barcelona, Universitat Politecnica de Catalunya".
"satisfaction"

#' Data: satisfaction including gender
#'
#' This data set contains the variables from a customer satisfaction study of 
#' a Spanish credit institution on 250 customers. The data is taken from the
#' \href{https://github.com/gastonstat/plspm}{plspm} package. For convenience, 
#' there is a version of the dataset with the last column (`gender`) removed: [satisfaction].
#'
#' @docType data
#'
#' @description  A data frame with 250 observations and 28 variables. 
#' Variables from 1 to 27 refer to six latent concepts: \code{IMAG}=Image, 
#' \code{EXPE}=Expectations, \code{QUAL}=Quality, \code{VAL}=Value, 
#' \code{SAT}=Satisfaction, and \code{LOY}=Loyalty.
#' \describe{
#'   \item{imag1-imag5}{Indicators attached to concept `IMAG` which is supposed to
#'                      capture aspects such as the institutions reputation, 
#'                      trustworthiness, seriousness, solidness, and caring 
#'                      about customer.}
#'   \item{expe1-expe5}{Indicators attached to concept `EXPE` which is supposed to
#'                      capture aspects concerning products and 
#'                      services provided, customer service, providing solutions,
#'                      and expectations for the overall quality.}
#'   \item{qual1-qual5}{Indicators attached to concept `QUAL` which is supposed to
#'                      capture aspects concerning reliability of products and services, 
#'                      the range of products and services, personal advice, 
#'                      and overall perceived quality.}
#'   \item{val1-val4}{Indicators attached to concept `VAL` which is supposed to
#'                      capture aspects related to beneficial services and 
#'                      products, valuable investments, quality relative to 
#'                      price, and price relative to quality.}
#'   \item{sat1-sat4}{Indicators attached to concept `SAT` which is supposed to
#'                      capture aspects concerning overall rating of satisfaction, 
#'                      fulfillment of expectations, satisfaction relative to 
#'                      other banks, and performance relative to customer's 
#'                      ideal bank.}
#'   \item{loy1-loy4}{Indicators attached to concept `LOY` which is supposed to
#'                    capture aspects concerning propensity to choose the 
#'                    same bank again, propensity to switch to other bank, 
#'                    intention to recommend the bank to friends, 
#'                    and the sense of loyalty.}
#'   \item{gender}{The sex of the respondent.}
#' }
#' @source The \href{https://github.com/gastonstat/plspm}{plspm} package (version  0.4.9). 
#' Original source according to \pkg{plspm}:
#' "Laboratory of Information Analysis and Modeling (LIAM). 
#' Facultat d'Informatica de Barcelona, Universitat Politecnica de Catalunya".
"satisfaction_gender"

#' Data: threecommonfactors
#'
#' A dataset containing 500 standardized observations on 9 indicator generated from a 
#' population model with three concepts modeled as common factors.
#' 
#' @docType data
#' 
#' @format A matrix with 400 rows and 9 variables:
#' \describe{
#'   \item{y11-y13}{Indicators attachted to the first common factor (`eta1`). 
#'                  Population loadings are: 0.7; 0.7; 0.7}
#'   \item{y21-y23}{Indicators attachted to the second common factor (`eta2`).
#'                  Population loadings are: 0.5; 0.7; 0.8}
#'   \item{y31-y33}{Indicators attachted to the third common factor (`eta3`).
#'                  Population loadings are: 0.8; 0.75; 0.7}
#' }
#'                  
#' The model is:
#' \deqn{`eta2` = gamma1 * `eta1` + zeta1}
#' \deqn{`eta3` = gamma2 * `eta1` + beta * `eta2` + zeta2}
#'
#' with population values `gamma1` = 0.6, `gamma2` = 0.4 and `beta` = 0.35.
#' @examples 
#' #============================================================================
#' # Correct model (the model used to generate the data)
#' #============================================================================
#' model_correct <- "
#' # Structural model
#' eta2 ~ eta1
#' eta3 ~ eta1 + eta2
#' 
#' # Measurement model
#' eta1 =~ y11 + y12 + y13
#' eta2 =~ y21 + y22 + y23
#' eta3 =~ y31 + y32 + y33 
#' "
#' 
#' a <- csem(threecommonfactors, model_correct)
#' 
#' ## The overall model fit is evidently almost perfect:
#' testOMF(a, .R = 50) # .R = 50 to speed up the example
"threecommonfactors"

#' Data: political democracy
#'
#' The Industrialization and Political Democracy dataset. This dataset is
#' used throughout Bollen's 1989 book (see pages 12, 17, 36 in chapter 2, pages
#' 228 and following in chapter 7, pages 321 and following in chapter 8; 
#' \insertCite{Bollen1989;textual}{cSEM}).
#' The dataset contains various measures of political democracy and
#' industrialization in developing countries.
#' 
#' @docType data
#' 
#' @format A data frame of 75 observations of 11 variables.
#'  \describe{
#'    \item{\code{y1}}{Expert ratings of the freedom of the press in 1960}
#'    \item{\code{y2}}{The freedom of political opposition in 1960}
#'    \item{\code{y3}}{The fairness of elections in 1960}
#'    \item{\code{y4}}{The effectiveness of the elected legislature in 1960}
#'    \item{\code{y5}}{Expert ratings of the freedom of the press in 1965}
#'    \item{\code{y6}}{The freedom of political opposition in 1965}
#'    \item{\code{y7}}{The fairness of elections in 1965}
#'    \item{\code{y8}}{The effectiveness of the elected legislature in 1965}
#'    \item{\code{x1}}{The gross national product (GNP) per capita in 1960}
#'    \item{\code{x2}}{The inanimate energy consumption per capita in 1960}
#'    \item{\code{x3}}{The percentage of the labor force in industry in 1960}
#' }
#'
#' @source The \href{http://lavaan.ugent.be/}{lavaan} package (version 0.6-3).
#' @references
#'   \insertAllCited{}            
#' @examples 
#' #============================================================================
#' # Example is taken from the lavaan website
#' #============================================================================
#' model <- "
#'   # measurement model
#'     ind60 =~ x1 + x2 + x3
#'     dem60 =~ y1 + y2 + y3 + y4
#'     dem65 =~ y5 + y6 + y7 + y8
#'     
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#'     
#'   # residual correlations
#'     y1 ~~ y5
#'     y2 ~~ y4 + y6
#'     y3 ~~ y7
#'     y4 ~~ y8
#'     y6 ~~ y8
#' "
#' 
#' aa <- csem(PoliticalDemocracy, model)
"PoliticalDemocracy"

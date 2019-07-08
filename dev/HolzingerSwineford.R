#' Holzinger and Swineford data
#'
#' The Holzinger and Swineford dataset \insertCite{Holzinger1939}{cSEM} 
#' consists of mental
#' ability test scores of seventh- and eighth-grade children from two
#' different schools (Pasteur and Grant-White). In the original dataset
#' (available in the \pkg{MBESS} package), there are scores for 26 tests.
#' The present dataset is taken from \pkg{lavaan} and consists of a smaller 
#' subset with 9 variables that is more widely used in the 
#' literature (e.g., \insertCite{Joereskog1969;textual}{cSEM}). 
#'
#' @docType data
#'
#' @format A data frame with 301 observations of 15 variables.
#' \describe{
#'   \item{\code{id}}{Identifier}
#'   \item{\code{sex}}{Gender}
#'   \item{\code{ageyr}}{Age, year part}
#'   \item{\code{agemo}}{Age, month part}
#'   \item{\code{school}}{School (Pasteur or Grant-White)}
#'   \item{\code{grade}}{Grade}
#'   \item{\code{x1}}{Visual perception}
#'   \item{\code{x2}}{Cubes}
#'   \item{\code{x3}}{Lozenges}
#'   \item{\code{x4}}{Paragraph comprehension}
#'   \item{\code{x5}}{Sentence completion}
#'   \item{\code{x6}}{Word meaning}
#'   \item{\code{x7}}{Speeded addition}
#'   \item{\code{x8}}{Speeded counting of dots}
#'   \item{\code{x9}}{Speeded discrimination straight and curved capitals}
#' }
#' 
#' @source The \href{http://lavaan.ugent.be/}{lavaan} package (version 0.6-3).
#' @references
#'   \insertAllCited{}
#'   
#' @examples 
#' #============================================================================
#' # Correct model (the model used to generate the data)
#' #============================================================================
#' model <- '
#' # Structural model
#'   textual ~ visual
#'   speed   ~ visual + textual
#' 
#' # Measurement model
#'   visual  =~ x1 + x2 + x3 
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9 
#' '
#' 
#' aa <- csem(HolzingerSwineford1939, model, .estimate_structural = FALSE)
"HolzingerSwineford1939"
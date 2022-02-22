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
#' @format A matrix with 500 rows and 9 variables:
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
#' testOMF(a, .R = 30) # .R = 30 to speed up the example
"threecommonfactors"

#' Data: Second order common factor of composites
#'
#' A dataset containing 500 standardized observations on 19 indicator generated from a 
#' population model with 6 concepts, three of which (`c1-c3`) are composites 
#' forming a second order common factor (`c4`). The remaining two (`eta1`, `eta2`)
#' are concepts modeled as common factors .
#' 
#' @docType data
#' 
#' @format A matrix with 500 rows and 19 variables:
#' \describe{
#'   \item{y11-y12}{Indicators attached  to `c1`. 
#'                  Population weights are: 0.8; 0.4.
#'                  Population loadings are: 0.925; 0.65}
#'   \item{y21-y24}{Indicators attached  to `c2`.
#'                  Population weights are: 0.5; 0.3; 0.2; 0.4.
#'                  Population loadings are: 0.804; 0.68; 0.554; 0.708}
#'   \item{y31-y38}{Indicators attached  to `c3`.
#'                  Population weights are: 0.3; 0.3; 0.1; 0.1; 0.2; 0.3; 0.4; 0.2.
#'                  Population loadings are: 0.496; 0.61; 0.535; 0.391; 0.391; 0.6; 0.5285; 0.53}
#'   \item{y41-y43}{Indicators attached  to `eta1`.
#'                  Population loadings are: 0.8; 0.7; 0.7}      
#'   \item{y51-y53}{Indicators attached  to `eta1`.
#'                  Population loadings are: 0.8; 0.8; 0.7}           
#' }
#'                  
#' The model is:
#' \deqn{`c4` = gamma1 * `eta1` + zeta1}
#' \deqn{`eta2` = gamma2 * `eta1` + beta * `c4` + zeta2}
#'
#' with population values `gamma1` = 0.6, `gamma2` = 0.4 and `beta` = 0.35.
#' The second order common factor is
#' \deqn{`c4` = lambdac1 * `c1` + lambdac2 * `c2` + lambdac3 * `c3` + epsilon}
"dgp_2ndorder_cf_of_c"

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
#' @source The \href{https://lavaan.ugent.be/}{lavaan} package (version 0.6-3).
#' @references
#'   \insertAllCited{}            
#' @examples 
#' #============================================================================
#' # Example is taken from the lavaan website
#' #============================================================================
#' # Note: example is modified. Across-block correlations are removed
#' model <- "
#' # Measurement model
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   
#' # Regressions / Path model
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#'   
#' # residual correlations
#'   y2 ~~ y4
#'   y6 ~~ y8
#' "
#' 
#' aa <- csem(PoliticalDemocracy, model)
"PoliticalDemocracy"


#' Data: Anime
#'
#' The data set for the example on \href{https://github.com/ISS-Analytics/pls-predict}{github.com/ISS-Analytics/pls-predict}
#' with irrelevant variables removed.
#'
#' @docType data
#'
#' @description A data frame with 183 observations and 13 variables. 
#' 
#' @source Original source: \href{https://github.com/ISS-Analytics/pls-predict}{github.com/ISS-Analytics/pls-predict}
"Anime"



#' Data: Russett
#'
#' The dataset was initially compiled by \insertCite{Russett1964;textual}{cSEM}, 
#' discussed and reprinted by \insertCite{Gifi1990;textual}{cSEM}, 
#' and partially transformed by \insertCite{Tenenhaus2011;textual}{cSEM}.
#' It is also used in \insertCite{Henseler2021;textual}{cSEM} for demonstration 
#' purposes.
#'
#' @format A data frame containing the following variables for 47 countries:
#'  \describe{
#'    \item{\code{gini}}{The Gini index of concentration}
#'    \item{\code{farm}}{The percentage of landholders who collectively occupy
#'    one-half of all the agricultural land (starting with the farmers
#'    with the smallest plots of land and working toward the largest)}
#'    \item{\code{rent}}{The percentage of the total number of farms that rent all
#'    their land. Transformation: ln (x + 1)}
#'    \item{\code{gnpr}}{The 1955 gross national product per capita in U.S. dollars.
#'    Transformation: ln (x)}
#'    \item{\code{labo}}{The percentage of the labor force employed in agriculture.
#'    Transformation: ln (x)}
#'    \item{\code{inst}}{Instability of personnel based on the term of office of the
#'    chief executive. Transformation: exp (x - 16.3)}
#'    \item{\code{ecks}}{The total number of politically motivated violent incidents,
#'    from plots to protracted guerrilla warfare. Transformation: ln (x + 1)}
#'    \item{\code{deat}}{The number of people killed as a result of internal group 
#'    violence per 1,000,000 people. Transformation: ln (x + 1)}
#'    \item{\code{stab}}{One if the country has a stable democracy, and zero otherwise}
#'    \item{\code{dict}}{One if the country experiences a dictatorship, and zero otherwise}
#' }
#' 
#' @docType data
#'
#' @description A data frame containing 10 variables with 47 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2020)
#' #============================================================================
#' model_Russett="
#' # Composite model
#' AgrIneq <~ gini + farm + rent
#' IndDev  <~ gnpr + labo
#' PolInst <~ inst + ecks + deat + stab + dict
#' 
#' # Structural model
#' PolInst ~ AgrIneq + IndDev
#' "
#' 
#' out <- csem(.data = Russett, .model = model_Russett,
#'             .PLS_weight_scheme_inner = 'factorial',
#'             .tolerance = 1e-06
#' )
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source From: \insertCite{Henseler2021;textual}{cSEM}
"Russett"


#' Data: ITFlex
#'
#' The dataset was studied by \insertCite{Benitez2018;textual}{cSEM} and is used in 
#' \insertCite{Henseler2021;textual}{cSEM} for demonstration purposes, see the
#' corresponding tutorial. 
#' All questionnaire items are measured on a 5-point scale.
#' 
#' @format A data frame containing the following variables:
#'  \describe{
#'    \item{\code{ITCOMP1}}{Software applications can be easily transported and
#'    used across multiple platforms.}
#'    \item{\code{ITCOMP2}}{Our firm provides multiple interfaces or entry points
#'    (e.g., web access) for external end users.}
#'    \item{\code{ITCOMP3}}{Our firm establishes corporate rules and standards for
#'    hardware and operating systems to ensure platform compatibility.}
#'    \item{\code{ITCOMP4}}{Data captured in one part of our organization are 
#'    immediately available to everyone in the firm.}
#'    \item{\code{ITCONN1}}{Our organization has electronic links and connections
#'    throughout the entire firm.}
#'    \item{\code{ITCONN2}}{Our firm is linked to business partners through
#'    electronic channels (e.g., websites, e-mail, wireless devices, electronic data interchange).}
#'    \item{\code{ITCONN3}}{All remote, branch, and mobile offices are connected to 
#'    the central office.}
#'    \item{\code{ITCONN4}}{There are very few identifiable communications
#'    bottlenecks within our firm.}
#'    \item{\code{MOD1}}{Our firm possesses a great speed in developing new
#'    business applications ormodifying existing applications.}
#'    \item{\code{MOD2}}{Our corporate database is able to communicate in
#'    several different protocols.}
#'    \item{\code{MOD3}}{Reusable software modules are widely used in new
#'    systems development.}
#'    \item{\code{MOD4}}{IT personnel use object-oriented and prepackaged
#'    modular tools to create software applications.}
#'    \item{\code{ITPSF1}}{Our IT personnel have the ability to work effectively in
#'    cross-functional teams.}
#'    \item{\code{ITPSF2}}{Our IT personnel are able to interpret business problems
#'    and develop appropriate technical solutions.}
#'    \item{\code{ITPSF3}}{Our IT personnel are self-directed and proactive.}
#'    \item{\code{ITPSF4}}{Our IT personnel are knowledgeable about the key
#'    success factors in our firm.}
#' }
#' 
#' @docType data
#'
#' @description A data frame containing 16 variables with 100 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2020)
#' #============================================================================
#' model_IT_Fex="
#' # Composite models
#' ITComp  <~ ITCOMP1 + ITCOMP2 + ITCOMP3 + ITCOMP4
#' Modul   <~ MOD1 + MOD2 + MOD3 + MOD4
#' ITConn  <~ ITCONN1 + ITCONN2 + ITCONN3 + ITCONN4
#' ITPers  <~ ITPSF1 + ITPSF2 + ITPSF3 + ITPSF4
#' 
#' # Saturated structural model
#' ITPers ~ ITComp + Modul + ITConn
#' Modul  ~ ITComp + ITConn 
#' ITConn ~ ITComp 
#' "
#' 
#'out <- csem(.data = ITFlex, .model = model_IT_Fex,
#'            .PLS_weight_scheme_inner = 'factorial',
#'            .tolerance = 1e-06,
#'            .PLS_ignore_structural_model = TRUE)
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source The data was collected through a survey by \insertCite{Benitez2018;textual}{cSEM}.
"ITFlex"



#' Data: LancelotMiltgenetal2016
#'
#' The data was analysed by \insertCite{Lancelot-Miltgen2016;textual}{cSEM} 
#' to study young consumers’ adoption intentions of a location tracker technology 
#' in the light of privacy concerns. It is also used in 
#' \insertCite{Henseler2021;textual}{cSEM} for demonstration purposes, see the
#' corresponding tutorial.
#' 
#' @docType data
#'
#' @description A data frame containing 10 variables with 1090 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2020)
#' #============================================================================
#' model_Med <- "
#' # Reflective measurement model
#' Trust =~ trust1 + trust2
#' PrCon =~ privcon1 + privcon2 + privcon3 + privcon4
#' Risk  =~ risk1 + risk2 + risk3
#' Int   =~ intent1 + intent2
#' 
#' # Structural model
#' Int   ~ Trust + PrCon + Risk
#' Risk  ~ Trust + PrCon
#' Trust ~ PrCon
#' "
#' 
#' out <- csem(.data = LancelotMiltgenetal2016, .model = model_Med,
#'             .PLS_weight_scheme_inner = 'factorial',
#'             .tolerance = 1e-06
#' )
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source This data has been collected through a cooperation with the European Commission 
#' Joint Research Center Institute for Prospective Technological Studies, contract 
#' “Young People and Emerging Digital Services: An Exploratory Survey on Motivations, 
#' Perceptions, and Acceptance of Risk” (EC JRC Contract IPTS No: 150876-2007 F1ED-FR).
"LancelotMiltgenetal2016"


#' Data: Yooetal2000
#'
#' The data is simulated and has the identical correlation matrix as the data
#' that was analysed by \insertCite{Yoo2000;textual}{cSEM} 
#' to examine how five elements of the marketing mix, namely price, store
#' image, distribution intensity, advertising spending, and price deals, are
#' related to the so-called dimensions of brand equity, i.e., perceived brand
#' quality, brand loyalty, and brand awareness/associations. It is also used in 
#' \insertCite{Henseler2017;textual}{cSEM} and \insertCite{Henseler2021;textual}{cSEM} 
#' for demonstration purposes, see the corresponding tutorial.
#' 
#' @docType data
#'
#' @description A data frame containing 34 variables with 569 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2021)
#' #============================================================================
#' model_HOC="
#' # Measurement models FOC
#' PR =~ PR1 + PR2 + PR3
#' IM =~ IM1 + IM2 + IM3
#' DI =~ DI1 + DI2 + DI3
#' AD =~ AD1 + AD2 + AD3
#' DL =~ DL1 + DL2 + DL3
#' AA =~ AA1 + AA2 + AA3 + AA4 + AA5 + AA6
#' LO =~ LO1 + LO3
#' QL =~ QL1 + QL2 + QL3 + QL4 + QL5 + QL6
#' 
#' # Composite model for SOC
#' BR <~ QL + LO + AA
#' 
#' # Structural model
#' BR~ PR + IM + DI + AD + DL 
#' "
#' 
#' out <- csem(.data = Yooetal2000, .model = model_HOC,
#'             .PLS_weight_scheme_inner = 'factorial',
#'             .tolerance = 1e-06)
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source Simulated data with the same correlation matrix as the data studied 
#' by \insertCite{Yoo2000;textual}{cSEM}.  
"Yooetal2000"


#' Data: Switching
#'
#' The data contains variables about the consumers’ intention to switch a
#' service provider. It is also used in \insertCite{Henseler2021;textual}{cSEM} 
#' for demonstration purposes, see the corresponding tutorial.
#' 
#' @docType data
#'
#' @description A data frame containing 26 variables with 767 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2021)
#' #============================================================================
#' model_Int <-"
#' # Measurement models
#' INV =~ INV1 + INV2 + INV3 +INV4
#' SAT =~ SAT1 + SAT2 + SAT3
#' INT =~ INT1 + INT2
#' 
#' # Structural model containing an interaction term.
#' INT ~ INV + SAT + INV.SAT
#' "
#' 
#' out <- csem(.data = Switching, .model = model_Int,
#'             .PLS_weight_scheme_inner = 'factorial',
#'             .tolerance = 1e-06)
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source The dataset is provided by Jörg Henseler.  
"Switching"



#' Data: BergamiBagozzi2000
#'
#' The dataset contains 22 variables and originates 
#' from a larger survey among South Korean employees conducted and
#' reported by \insertCite{Bergami2000;textual}{cSEM}. It is
#' also used in  \insertCite{Hwang2004;textual}{cSEM} and 
#' \insertCite{Henseler2021;textual}{cSEM} 
#' for demonstration purposes, see the corresponding tutorial.
#' 
#' @docType data
#'
#' @description A data frame containing 22 variables with 305 observations. 
#' 
#' @examples 
#' #============================================================================
#' # Example is taken from Henseler (2021)
#' #============================================================================
#' model_Bergami_Bagozzi_Henseler="
#' # Measurement models
#' OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
#' OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
#' AffLove =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
#' AffJoy  =~ orgcmt5 + orgcmt8
#' Gender  <~ gender
#' 
#' # Structural model 
#' OrgIden ~ OrgPres
#' AffLove ~ OrgPres + OrgIden + Gender 
#' AffJoy  ~ OrgPres + OrgIden + Gender 
#' "
#' 
#' out <- csem(.data = BergamiBagozzi2000, 
#'             .model = model_Bergami_Bagozzi_Henseler,
#'             .PLS_weight_scheme_inner = 'factorial',
#'             .tolerance = 1e-06
#' )
#' 
#' #============================================================================
#' # Example is taken from Hwang et al. (2004)
#' #============================================================================ 
#' 
#' model_Bergami_Bagozzi_Hwang="
#' # Measurement models
#' OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
#' OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
#' AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
#' AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8
#' 
#' # Structural model 
#' OrgIden ~ OrgPres 
#' AffLove ~ OrgIden
#' AffJoy  ~ OrgIden"
#'
#'out_Hwang <- csem(.data = BergamiBagozzi2000, 
#'                  .model = model_Bergami_Bagozzi_Hwang,
#'                  .approach_weights = "GSCA",
#'                  .disattenuate = FALSE,
#'                  .id = "gender",
#'                  .tolerance = 1e-06) 
#'
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source Survey among South Korean employees conducted and
#' reported by \insertCite{Bergami2000;textual}{cSEM}. 
"BergamiBagozzi2000"


#' Data: Summers
#'
#' The indicator correlation matrix for a modified version of \insertCite{Summers1965;textual}{cSEM} 
#' model. All constructs are modeled as composites. 
#' 
#' @docType data
#'
#' @description A (18 x 18) indicator correlation matrix. 
#' 
#' @examples 
#' 
#' require(cSEM)
#' 
#' model <- "
#' ETA1 ~ ETA2 + XI1 + XI2
#' ETA2 ~ ETA1 + XI3 +XI4
#' 
#' ETA1 ~~ ETA2
#' 
#' XI1  <~ x1 + x2 + x3
#' XI2  <~ x4 + x5 + x6
#' XI3  <~ x7 + x8 + x9
#' XI4  <~ x10 + x11 + x12
#' ETA1 <~ y1 + y2 + y3
#' ETA2 <~ y4 + y5 + y6
#' "
#' 
#' ## Generate data
#' summers_dat <- MASS::mvrnorm(n = 300, mu = rep(0, 18), 
#'                              Sigma = Sigma_Summers_composites, empirical = TRUE)
#' 
#' ## Estimate
#' res <- csem(.data = summers_dat, .model = model) # inconsistent
#' 
#' ## 
#' # 2SLS
#' res_2SLS <- csem(.data = summers_dat, .model = model, .approach_paths = "2SLS",
#'                  .instruments = list(ETA1 = c('XI1', 'XI2', 'XI3', 'XI4'),
#'                                      ETA2 = c('XI1', 'XI2', 'XI3', 'XI4'))
#')
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source Own calculation based on \insertCite{Dijkstra2015;textual}{cSEM}.
"Sigma_Summers_composites"

#' Data: SQ
#'
#' The data comes from a European manufacturer of durable consumer goods and was 
#' studied by \insertCite{Bliemel2004;textual}{cSEM} who focused on service quality.
#' It is also used in \insertCite{Henseler2021;textual}{cSEM} 
#' for demonstration purposes, see the corresponding tutorial.
#' 
#' @docType data
#'
#' @description A data frame containing 23 variables with 411 observations. The original
#' indicators were measured on a 6-point scale. In this version of the dataset,
#' the indidcators are scaled to be between 0 and 100. 
#' 
#' 
#' @references
#'   \insertAllCited{}
#'     
#' @source The dataset is provided by Jörg Henseler.  
"SQ"
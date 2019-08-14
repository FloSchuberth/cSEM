#' GSCA models concepts as composites, i.e. as an exact linear 
#' combination of indicators. The main idea of GSCA is to 
#' divide the overall structural model into three submodels: the structural model, 
#' the measurement model and the weighted relation model. These sub-models and their
#' corresponding matrices are provided by the `.csem_model` argument.
#' 
#' Weights as well as structural coefficients and loadings are estimated with the 
#' Alternating Least Squares (ALS) algorithm.
#' This algorithm consists of two main steps. In the first step, structural coefficients
#' and loadings are updated for fixed weights. Subsequently, in the second step, weights
#' are updated keeping the other parameters constant. In doing so, the sum squared residuals 
#' of all equations is minimized.
#' 
#' These two steps are iterated until convergence is reached or the prescribed 
#' maximum number of iteration steps, '.iter_max', was carried out. In their example
#' in \insertCite{Hwang2014;textual}{cSEM}, p. 75, the authors set '.iter_max' 
#' equal to 100, but this is an arbitrary choice. The user might increase the maximum number of iterations
#' (.iter_max) if the algorithm has not converged (indicated by a `Conv_status` equal 
#' to `FALSE`). 
#' 
#' Several convergence criterion are possible, all consisting in calculating some
#' difference of the weight matrices W of two subsequent steps. Checking for convergence
#' this way results in the same interpretation as calculating the differences
#' the global optimization criterion values as proposed by Hwang and Takane.
#' The `.tolerance` argument provides the value below which the difference has to
#' fall such that the algorithm converges. The default value of `1e-05` is the same
#' value which Hwang and Takane use in their example in \insertCite{Hwang2014;textual}{cSEM}
#' , p. 75. Nevertheless, the argument `.tolerance` can take any other value as well.
#' 
#' The ALS algorithm as outlined in \insertCite{Hwang2014;textual}{cSEM} essentially 
#' consists of regressions of endogenous constructs on 
#' exogenous constructs for the structural model coefficients and of regressions 
#' of all constructs on their related indicators for the loadings 
#' (steps 1a and 1b). If a construct is not a common factor/latent variable,  but a pure composite
#' the corresponding "construct" loadings do not need to be calculated because the 
#' constructs are exact linear combinations of their indicators. Instead, "composite" 
#' loadings are calculated and provided as the estimated loadings. These "composite" 
#' loadings are the "true" correlations between a construct and indicator. This correlation 
#' always exists, no matter which is the type of the construct under consideration.
#' Considering common factors/latent variables, estimated ("construct") loadings 
#' are the correlations between indicators and proxies, but not the "true" correlation
#' between indicator and the latent variable because of attenuation.
#' 
#' In the case that there is at least one construct which is not modeled as a 
#' common factor but as a pure composite, the user imperatively has to use GSCA 
#' and GSCA_M is no option. The reason is that the implementation of GSCA_M in 
#' the cSEM-package is based on \insertCite{Hwang2017;textual}{cSEM}. 
#' The GSCA_M-approach presented by the authors only works out if all constructs
#' are of common factor type. Otherwise, i.e., if there is at least one construct
#' of composite type estimation using GSCA_M fails because calculating weight estimates 
#' with GSCA_M leads to a product involving the measurement matrix. 
#' This matrix does not have full rank in the case of at least one pure composite construct.
#' It has a zero row for every construct which is a pure composite (i.e. all related loadings are zero) 
#' and, therefore, leads to a non-invertible matrix when multiplying it with its transposed.
#' 
#' Therefore estimation must be carried out with GSCA if there is at least one composite
#' construct. Otherwise, i.e., if there are only constructs which are of common factor type,
#' calling [csem()] will lead to an estimation via GSCA_M except in the case that
#' the user explicitly sets the argument `.disattenuate` to `FALSE`. estimation is always done
#' by 'standard' GSCA.
#' 
#' 
#' 
#' 
#' 
#' 
#' #' The GSCA_M procedure is an alternative approach to PLSc in structural equation 
#' modeling with composites where indicators (or/and constructs) are observed with 
#' errors. In this case parameter estimators, especially loadings, might be biased, 
#' when using GSCA or PLS-PM for estimation. 
#' However, GSCA_M provides a way to estimate the parameters consistently also in 
#' this situation. The term 'GSCA_M' stands for 'GSCA with measurement errors incorporated'. 
#' This approach was first presented in \insertCite{Hwang2017;textual}{cSEM}.
#' 
#' The basic idea of GSCA_M is to model indicators in the measurement model as a 
#' combination of common parts (arising from the constructs) and unique parts.
#' The purpose of adding a unique part to each indicator is to account for measurement 
#' errors in the indicators. In a next step, latent variables are expressed in 
#' the weighted relation model as a linear combination of indicators but with their
#' unique parts removed.
#' 
#' Together with the structural model, all three sub-models are combined to derive
#' one model equation. The corresponding matrices of the sub-models are provided 
#' by the `.csem_model` argument. This leads to a global optimization criterion 
#' which is minimized in order to obtain the optimal estimators for all parameters. 
#' This is done via an iterative algorithm. Having assigned initial values to all
#' parameters, several alternating steps are carried out until convergence. Each 
#' step consists essentially of regressions such that every set of parameters is 
#' updated by a least squares calculation keeping the other parameters constant.
#' This algorithm is explained more detailly in the Appendix of \insertCite{Hwang2017;textual}{cSEM}.
#' 
#' The steps are iterated until convergence is reached or the prescribed 
#' maximum number of iteration steps, '.iter_max', was carried out. In their example
#' in \insertCite{Hwang2014;textual}{cSEM}, p. 75, the authors set '.iter_max' 
#' equal to 100 for an estimation with GSCA, but this is an arbitrary choice. 
#' The user might try several values for '.iter_max' and eventually increase its 
#' value if convergence happened previously due to reaching the maximum number 
#' of steps (indicated by a `Conv_status` equal to `FALSE`). 
#' 
#' The convergence criterion for GSCA_M is some difference of subsequent parameter
#' estimators involving those of loadings and path coefficients. In this case, 
#' convergence means that this difference falls below a certain value provided by 
#' the `.tolerance` argument. The default value of `1e-05` is the same value which 
#' Hwang and Takane use in their example in \insertCite{Hwang2014;textual}{cSEM}, p. 75. 
#' Nevertheless, the argument `.tolerance` can take any other value as well.
#' 
#' Parameters are automatically estimated via GSCA_M when calling [csem()]. However, 
#' if there is at least one construct which is not a common factor, but a pure composite,
#' parameters have to be estimated with GSCA. The reason is that the implementation 
#' of GSCA_M in the cSEM-package is based on \insertCite{Hwang2017;textual}{cSEM}.
#' In the case of GSCA_M, the estimated loadings are the correlations between indicators 
#' and proxies, not the "true" correlations between indicators and underlying but unknown
#' constructs.
#' The GSCA_M-approach presented by the authors only works out if all constructs
#' are of common factor type. Otherwise, i.e., if there is at least one construct
#' of composite type estimation using GSCA_M fails because calculating weight estimates 
#' with GSCA_M leads to a product involving the measurement matrix. This matrix 
#' does not have full rank in the case of at least one pure composite construct.
#' It has a zero row for every construct which is a pure composite (i.e. all related loadings are zero) 
#' and, therefore, leads to a non-invertible matrix when multiplying it with its transposed.
#' Thus, in this case, the user imperatively has to use GSCA and GSCA_M is no option.
#' Otherwise, i.e., if there are only constructs which are of common factor type,
#' calling [csem()] will lead to an estimation via GSCA_M except in the case that
#' the user explicitly sets the argument `.disattenuate` to `FALSE`. estimation is always done
#' by 'standard' GSCA.
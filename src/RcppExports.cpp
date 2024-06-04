// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kroneckerC
arma::mat kroneckerC(arma::mat& X, arma::mat& Y, arma::uvec& idx);
RcppExport SEXP _cSEM_kroneckerC(SEXP XSEXP, SEXP YSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(kroneckerC(X, Y, idx));
    return rcpp_result_gen;
END_RCPP
}
// kroneckerCsp
arma::sp_mat kroneckerCsp(arma::mat& X, arma::mat& Y, arma::uvec& idx);
RcppExport SEXP _cSEM_kroneckerCsp(SEXP XSEXP, SEXP YSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(kroneckerCsp(X, Y, idx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cSEM_kroneckerC", (DL_FUNC) &_cSEM_kroneckerC, 3},
    {"_cSEM_kroneckerCsp", (DL_FUNC) &_cSEM_kroneckerCsp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cSEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

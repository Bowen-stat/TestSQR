// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// soft
arma::mat soft(arma::mat X, double lambda);
RcppExport SEXP _TestSQR_soft(SEXP XSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(soft(X, lambda));
    return rcpp_result_gen;
END_RCPP
}
// projlasso
Rcpp::List projlasso(arma::mat Z, arma::mat covzx, arma::vec lambda, double err_abs, int maxIter, double rho);
RcppExport SEXP _TestSQR_projlasso(SEXP ZSEXP, SEXP covzxSEXP, SEXP lambdaSEXP, SEXP err_absSEXP, SEXP maxIterSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covzx(covzxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type err_abs(err_absSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(projlasso(Z, covzx, lambda, err_abs, maxIter, rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TestSQR_soft", (DL_FUNC) &_TestSQR_soft, 2},
    {"_TestSQR_projlasso", (DL_FUNC) &_TestSQR_projlasso, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_TestSQR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

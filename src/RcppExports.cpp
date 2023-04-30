// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// convolve_cpp
NumericVector convolve_cpp(const NumericVector a, const NumericVector b);
RcppExport SEXP _adept_convolve_cpp(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_cpp(a, b));
    return rcpp_result_gen;
END_RCPP
}
// sliding_cor_sd_cpp
NumericVector sliding_cor_sd_cpp(const NumericVector shortvec, const NumericVector longvec, double sd_shortvec, const NumericVector sd_longvec_current);
RcppExport SEXP _adept_sliding_cor_sd_cpp(SEXP shortvecSEXP, SEXP longvecSEXP, SEXP sd_shortvecSEXP, SEXP sd_longvec_currentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type shortvec(shortvecSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type longvec(longvecSEXP);
    Rcpp::traits::input_parameter< double >::type sd_shortvec(sd_shortvecSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type sd_longvec_current(sd_longvec_currentSEXP);
    rcpp_result_gen = Rcpp::wrap(sliding_cor_sd_cpp(shortvec, longvec, sd_shortvec, sd_longvec_current));
    return rcpp_result_gen;
END_RCPP
}
// sliding_cor_store_sd_cpp
List sliding_cor_store_sd_cpp(const NumericVector shortvec, const NumericVector longvec, double sd_shortvec);
RcppExport SEXP _adept_sliding_cor_store_sd_cpp(SEXP shortvecSEXP, SEXP longvecSEXP, SEXP sd_shortvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type shortvec(shortvecSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type longvec(longvecSEXP);
    Rcpp::traits::input_parameter< double >::type sd_shortvec(sd_shortvecSEXP);
    rcpp_result_gen = Rcpp::wrap(sliding_cor_store_sd_cpp(shortvec, longvec, sd_shortvec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_adept_convolve_cpp", (DL_FUNC) &_adept_convolve_cpp, 2},
    {"_adept_sliding_cor_sd_cpp", (DL_FUNC) &_adept_sliding_cor_sd_cpp, 4},
    {"_adept_sliding_cor_store_sd_cpp", (DL_FUNC) &_adept_sliding_cor_store_sd_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_adept(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gibbsC
NumericMatrix gibbsC(int N, int thin, int n, int aa, int bb);
RcppExport SEXP _StatComp21064_gibbsC(SEXP NSEXP, SEXP thinSEXP, SEXP nSEXP, SEXP aaSEXP, SEXP bbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< int >::type bb(bbSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsC(N, thin, n, aa, bb));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp21064_gibbsC", (DL_FUNC) &_StatComp21064_gibbsC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp21064(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
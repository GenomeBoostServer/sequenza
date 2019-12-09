#include <RcppArmadillo.h>

using namespace Rcpp;


// slide_matrix
DataFrame slide_matrix(NumericVector x, IntegerVector position, int w = 100, bool smooth = true, int method = 1, bool verbose = true);
RcppExport SEXP _sequenza_slide_matrix(SEXP xSEXP, SEXP positionSEXP, SEXP wSEXP, SEXP smoothSEXP, SEXP methodSEXP, SEXP verboseSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type position(positionSEXP);
    traits::input_parameter< int >::type w(wSEXP);
    traits::input_parameter< bool >::type smooth(smoothSEXP);
    traits::input_parameter< int >::type method(methodSEXP);
    traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = wrap(slide_matrix(x, position, w, smooth, method, verbose));
    return rcpp_result_gen;
END_RCPP
}


// get_peaks
IntegerVector get_peaks(arma::vec x, IntegerVector position, int w = 100);
RcppExport SEXP _sequenza_get_peaks(SEXP xSEXP, SEXP positionSEXP, SEXP wSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< arma::vec >::type x(xSEXP);
    traits::input_parameter< IntegerVector >::type position(positionSEXP);
    traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = wrap(get_peaks(x, position, w));
    return rcpp_result_gen;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"_sequenza_slide_matrix", (DL_FUNC) &_sequenza_slide_matrix, 6},
    {"_sequenza_get_peaks", (DL_FUNC) &_sequenza_get_peaks, 3},
    {NULL, NULL, 0}
};


RcppExport void R_init_sequenza(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
};

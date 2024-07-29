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


// rmvnorm
arma::mat rmvnorm (int n, arma::vec mean, arma::mat sigma);
RcppExport SEXP _sequenza_rmvnorm(SEXP nSEXP, SEXP meanSEXP, SEXP sigmaSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< int >::type n(nSEXP);
    traits::input_parameter< arma::vec >::type mean(meanSEXP);
    traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = wrap(rmvnorm(n, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}

// rwishart
arma::mat rwishart(int df,  arma::mat const& S);
RcppExport SEXP _sequenza_rwishart(SEXP dfSEXP, SEXP SSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< int >::type df(dfSEXP);
    traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = wrap(rwishart(df, S));
    return rcpp_result_gen;
END_RCPP
}

// ginv
arma::mat ginv(arma::mat m, Rcpp::Nullable<double> tol);
RcppExport SEXP _sequenza_ginv(SEXP mSEXP, SEXP tolSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< arma::mat >::type m(mSEXP);
    traits::input_parameter< Rcpp::Nullable<double> >::type tol(tolSEXP);
    rcpp_result_gen = wrap(ginv(m, tol));
    return rcpp_result_gen;
END_RCPP
}

// get_clusters
IntegerVector get_clusters(const arma::mat &traces, int tot_clusters, double discard = 0.1);
RcppExport SEXP _sequenza_get_clusters(SEXP tracesSEXP, SEXP tot_clustersSEXP, SEXP discardSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< arma::mat >::type traces(tracesSEXP);
    traits::input_parameter< int >::type tot_clusters(tot_clustersSEXP);
    traits::input_parameter< double >::type discard(discardSEXP);
    rcpp_result_gen = wrap(get_clusters(traces, tot_clusters, discard));
    return rcpp_result_gen;
END_RCPP
}

// get_clusters_stats
arma::mat get_clusters_stats(const arma::mat &traces, int tot_clusters, double discard = 0.1);
RcppExport SEXP _sequenza_get_clusters_stats(SEXP tracesSEXP, SEXP tot_clustersSEXP, SEXP discardSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< arma::mat >::type traces(tracesSEXP);
    traits::input_parameter< int >::type tot_clusters(tot_clustersSEXP);
    traits::input_parameter< double >::type discard(discardSEXP);
    rcpp_result_gen = wrap(get_clusters_stats(traces, tot_clusters, discard));
    return rcpp_result_gen;
END_RCPP
}

// sample_z_loop
List sample_z_loop(int N, arma::mat A, double alpha, arma::vec z, arma::mat mu_ar, arma::cube S_ar,
     arma::vec mu_a0, arma::mat R_a0, double beta_a0, arma::mat W_a0, int nclust);
RcppExport SEXP _sequenza_sample_z_loop(SEXP NSEXP, SEXP ASEXP, SEXP alphaSEXP, SEXP zSEXP,
    SEXP mu_arSEXP, SEXP S_arSEXP, SEXP mu_a0SEXP, SEXP R_a0SEXP, SEXP beta_a0SEXP,
    SEXP W_a0SEXP, SEXP nclustSEXP){
BEGIN_RCPP
    RObject rcpp_result_gen;
    RNGScope rcpp_rngScope_gen;
    traits::input_parameter< int >::type N(NSEXP);
    traits::input_parameter< arma::mat >::type A(ASEXP);
    traits::input_parameter< double >::type alpha(alphaSEXP);
    traits::input_parameter< arma::vec >::type z(zSEXP);
    traits::input_parameter< arma::mat >::type mu_ar(mu_arSEXP);
    traits::input_parameter< arma::cube >::type S_ar(S_arSEXP);
    traits::input_parameter< arma::vec >::type mu_a0(mu_a0SEXP);
    traits::input_parameter< arma::mat >::type R_a0(R_a0SEXP);
    traits::input_parameter< double >::type beta_a0(beta_a0SEXP);
    traits::input_parameter< arma::mat >::type W_a0(W_a0SEXP);
    traits::input_parameter< int >::type nclust(nclustSEXP);
    rcpp_result_gen = wrap(sample_z_loop(N, A, alpha, z, mu_ar, S_ar, mu_a0, R_a0, beta_a0, W_a0, nclust));
    return rcpp_result_gen;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"_sequenza_slide_matrix", (DL_FUNC) &_sequenza_slide_matrix, 6},
    {"_sequenza_get_peaks", (DL_FUNC) &_sequenza_get_peaks, 3},
    {"_sequenza_rmvnorm", (DL_FUNC) &_sequenza_rmvnorm, 3},
    {"_sequenza_rwishart", (DL_FUNC) &_sequenza_rwishart, 2},
    {"_sequenza_ginv", (DL_FUNC) &_sequenza_ginv, 2},
    {"_sequenza_get_clusters", (DL_FUNC) &_sequenza_get_clusters, 3},
    {"_sequenza_get_clusters_stats", (DL_FUNC) &_sequenza_get_clusters_stats, 3},
    {"_sequenza_sample_z_loop", (DL_FUNC) &_sequenza_sample_z_loop, 11},
    {NULL, NULL, 0}
};


RcppExport void R_init_sequenza(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
};

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>
// Correctly setup the build environment 
// [[Rcpp::depends(RcppArmadillo)]]
	
using namespace arma;
using namespace Rcpp;

#ifndef utils_h
#define utils_h

// Add include guards for all header files
#include <RcppArmadillo.h>
#include <Rcpp.h>

// Add inline keywords for small functions
inline int length_unique(arma::vec const & x);

// Add constexpr for compile-time constants
constexpr int DEFAULT_WINDOW_SIZE = 100;
constexpr double DEFAULT_DISCARD = 0.1;

IntegerVector tabulate_cpp (arma::vec const & x, unsigned int const max);
IntegerVector tabulate2 (arma::vec const & x);
IntegerVector unique_sort (arma::vec const & x);
IntegerVector concat_iv(IntegerVector a, IntegerVector b);
IntegerVector sample_r(IntegerVector v, int n, arma::vec p);
arma::vec dmvnorm(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd);
arma::mat rmvnorm (int n, arma::vec mean, arma::mat sigma);
arma::mat rwishart(int df,  arma::mat const& S);
arma::mat ginv(arma::mat m, Rcpp::Nullable<double> tol);
arma::cube nrwishart(int n, int df, const arma::mat& S);
DataFrame slide_matrix(NumericVector x, IntegerVector position, int w = 100, 
                      bool smooth = true, int method = 1, bool verbose = true);

#endif

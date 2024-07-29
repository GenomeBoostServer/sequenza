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

int length_unique (arma::vec const & x);
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

#endif

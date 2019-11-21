#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>
// Correctly setup the build environment 
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

// Implementation of a function that given 2 vectors
// returns the D value of the Kolmogorov-Smirnov test 
double kstest(const NumericVector& a, const NumericVector&  b) {
	NumericVector x = na_omit(a);
	NumericVector y = na_omit(b);
	int n1 = x.size(), n2 = y.size();
    std::sort(x.begin(),x.end());
    std::sort(y.begin(),y.end());
    int j1 = 0, j2 = 0;
    double d = 0.0, fn1 = 0.0, fn2 = 0.0, d1, d2, dtemp;
    while (j1<n1 && j2<n2) {
        d1 = x[j1];
        d2 = y[j2];
        if (d1 <= d2) {
            fn1 = (j1+1.0)/n1;
            j1++;
        }
        if(d2 <= d1){
            fn2 = (j2+1.0)/n2;
            j2++;
        }
        dtemp = std::abs(fn2-fn1);
        if (dtemp>d) {
            d=dtemp;
        }
    }
    return(d);
}


double diff_test(const NumericVector& x, int start, int end, int method = 1) {
    int N = end - start;
    int su = (int)std::round(N / 2);
    int sd = N - su;
    NumericVector up(su), down(sd);
    for (int i = start, ii = 0; i < end; i++) {
        if (ii < su) {
            up[ii] = x[i];
        } else {
            down[ii-su] = x[i];
        }
        ii++;
    }

    if (method == 1)  {
        return(kstest(up, down));
    } else if (method == 2) {
        return(std::abs(Rcpp::mean(up - down)));
    } else if (method == 3) {
        return((1 + kstest(up, down)) * std::abs(Rcpp::mean(up - down)));
    } else {
        return(kstest(up, down));
    }
}


void  coord_quart(NumericVector res, const IntegerVector& x, const NumericVector& y, int w_half, int method = 1) {
    int N = res.size();
    Progress p(N, true);
    for (int i = 0; i < N; i++){
        if ( ! Progress::check_abort() ) {
            p.increment();
            int x1 = x[i] - w_half;
            int x2 = x[i] + w_half - 1;
            res[i] = diff_test(y, x1, x2, method);
        }
    }
}


DataFrame slide_matrix(NumericVector x, IntegerVector position, int w = 100, bool smooth = true, int method = 1) {
    int N = x.size();
    //x = x / mean(x, na.rm = TRUE)
    int w_half = std::round(w / 2);
    IntegerVector steps = Rcpp::Range(w_half, N - w_half);
    NumericVector ksres(steps.size());
    coord_quart(ksres, steps, x, w_half, method);
    IntegerVector pos_sub =  position[steps];
    if (smooth) {
        // sloppy smoothing using R smooth.spline
        // hopefully replaceable with a c++ native function.
        Environment pkg = Environment::namespace_env("stats");
        Function smooth_spline = pkg["smooth.spline"];
        Rcpp::List smooth_ksres = smooth_spline(pos_sub, ksres);
        return(DataFrame::create( Named("x") = smooth_ksres["x"],  Named("y") = smooth_ksres["y"]));
    } else {
        return(DataFrame::create( Named("x") = pos_sub,  Named("y") = ksres));      
    }
}

IntegerVector get_peaks(arma::vec x, IntegerVector position, int w = 100) {
    int w_half = std::round(w / 2);
    int N = x.n_elem;
    IntegerVector res(N - w_half);
    int ri = 0;
    for (int i = w_half; i < (N - w_half); i++){
        int x1 = i - w_half;
        int x2 = i + w_half - 1;
        arma::vec xr = x(arma::span(x1, x2));
        int max_i = arma::index_max(xr);
        if (max_i == w_half) {
            res(ri) = position(i);
            ri++;
        }
    }
    IntegerVector res2(ri);
    for (int i = 0; i < ri; i++){
        res2(i) = res(i);
    }
    return(res2);
}
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <atomic>
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
    size_t n1 = x.size(), n2 = y.size();
    if(n1 > 10000) {
        #pragma omp parallel
        {
            #pragma omp single
            std::sort(x.begin(), x.end());
        }
    } else {
        std::sort(x.begin(), x.end());
    }
    if(n2 > 10000) {
        #pragma omp parallel
        {
            #pragma omp single
            std::sort(y.begin(), y.end());
        }
    } else {
        std::sort(y.begin(), y.end());
    }
    size_t j1 = 0, j2 = 0;
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


inline double diff_test(const NumericVector& x, int start, int end, int method = 1) {
    const int N = end - start;
    const int su = (int)std::round(N / 2);
    const int sd = N - su;
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


void  coord_quart(NumericVector res, const IntegerVector& x, const NumericVector& y, int w_half, int method = 1, bool verbose = true) {
    const int N = res.size();
    Progress p(N, verbose);
    std::atomic<int> counter(0);
    int step = N / 100;
    #pragma omp parallel for if(N > 1000) schedule(dynamic)
    for (int i = 0; i < N; i++){
        if (verbose == true) {
            if ( ! Progress::check_abort() ) {
                if (counter++ % step == 0) {
                    #pragma omp critical
                    p.increment(step);
                }
            }
        }
        int x1 = x[i] - w_half;
        int x2 = x[i] + w_half - 1;
        res[i] = diff_test(y, x1, x2, method);
    }
}

// Modify the function signature to match exactly what's in RcppExports.cpp
DataFrame slide_matrix(NumericVector x, IntegerVector position, int w = 100, bool smooth = true, int method = 1, bool verbose = true) {
    // Input validation
    if(x.size() != position.size()) {
        stop("x and position must have the same length");
    }
    if(w <= 0) {
        stop("Window size must be positive");
    }

    const int N = x.size();
    const int w_half = std::round(w / 2);
    
    // Pre-calculate size to avoid reallocation
    const int result_size = N - 2 * w_half;
    IntegerVector steps = Rcpp::Range(w_half, N - w_half - 1);
    NumericVector ksres(result_size);
    
    coord_quart(ksres, steps, x, w_half, method, verbose);
    IntegerVector pos_sub = position[steps];

    if (smooth) {
        // Use R's smooth.spline with optimized parameters
        Environment pkg = Environment::namespace_env("stats");
        Function smooth_spline = pkg["smooth.spline"];
        Function predict = pkg["predict"];
        
        // Calculate optimal smoothing parameters based on data size
        double spar = std::min(1.0, std::max(0.5, log10(result_size) / 10.0));
        
        // Fit the smooth spline
        List spline_fit = smooth_spline(Named("x") = pos_sub,
                                      Named("y") = ksres,
                                      Named("spar") = spar,
                                      Named("keep.data") = false);
        
        // Predict values at original x positions to maintain dimension consistency
        List predicted = predict(spline_fit, pos_sub);
        NumericVector smoothed_y = as<NumericVector>(predicted["y"]);
        
        return DataFrame::create(
            _["x"] = pos_sub,
            _["y"] = smoothed_y,
            _["raw_y"] = ksres
        );
    } else {
        return DataFrame::create(
            _["x"] = pos_sub,
            _["y"] = ksres
        );
    }
}

IntegerVector get_peaks(arma::vec x, IntegerVector position, int w = 100) {
    int w_half = std::round(w / 2);
    int N = x.n_elem;
    if (N <= w_half ) {
        w_half = N / 2;
    }
    IntegerVector res(N - w_half);
    int ri = 0;
    
    for (int i = w_half; i < (N - w_half); i++) {
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
    for (int i = 0; i < ri; i++) {
        res2(i) = res(i);
    }
    return res2;
}
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

// Add this helper function for C++ smoothing
inline NumericVector smooth_vector(const NumericVector& y, const IntegerVector& x, int window = 5) {
    const int n = y.size();
    NumericVector smoothed(n);
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n; i++) {
        int start = std::max(0, i - window/2);
        int end = std::min(n-1, i + window/2);
        double sum = 0.0;
        int count = 0;
        
        for(int j = start; j <= end; j++) {
            if(!R_IsNA(y[j])) {
                sum += y[j];
                count++;
            }
        }
        smoothed[i] = count > 0 ? sum/count : NA_REAL;
    }
    return smoothed;
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
    
    // Use parallel processing for larger datasets
    #pragma omp parallel sections if(N > 10000)
    {
        #pragma omp section
        {
            coord_quart(ksres, steps, x, w_half, method, verbose);
        }
    }
    
    IntegerVector pos_sub = position[steps];

    if (smooth) {
        // Use our C++ smoothing implementation instead of R's smooth.spline
        NumericVector smoothed_y = smooth_vector(ksres, pos_sub);
        return DataFrame::create(
            _["x"] = pos_sub,
            _["y"] = smoothed_y,
            _["raw_y"] = ksres // Optional: include raw values
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
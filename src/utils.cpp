#include "utils.h"

using namespace Rcpp;
using namespace std;

bool isSymmetric(arma::mat mat) 
{
    arma::mat tr = arma::trans(mat);
    int N = mat.n_cols;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (mat(i, j) != tr(i, j)) {
                return false;
            }
        }
    }
    return true; 
}

void makeSymmetric(arma::mat &X) 
{
    arma::mat tr = arma::trans(X);
    int N = X.n_cols;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (X(i, j) != tr(i, j)) {
                if (i > j) {
                    X(i, j) = tr(i, j);
                }
            }
        }
    }
}

bool chol_fb (arma::mat &R, arma::mat X) {
    bool status = false;
    NumericMatrix R2;
    Environment pkg = Environment::namespace_env("base");
    Function cholR = pkg["chol"];
    try {
        R2 = cholR(X);
        status = true;
        R = Rcpp::as<arma::mat>(R2);
    } catch(std::exception &ex) {
        warning("chol decomposition fallback failed!");
    }
    return status;
}

IntegerVector sample_r (IntegerVector v, int n, arma::vec p) {
    //Environment pkg = Environment::namespace_env("base");
    //Function sampleR = pkg["sample"];
    //return(sampleR(v, n, true, p));
    // UGLY, FIX implement with RcppArmadillo::sample....
    NumericVector pv = NumericVector(p.begin(), p.end());
    IntegerVector ret = Rcpp::sample(v, n, true, pv) ;
    return ret ;
}


arma::mat chol_fallback (arma::mat X) {
    arma::mat R;
    bool staus = chol_fb(R, X);
    return R;
}

arma::mat chol2 (arma::mat X) {
    bool check = false;
    arma::mat R, L, U, P;
    // Workaround to the warning of not symmetric matrix
    // added in RcppArmadillo_0.9.400.2.0
    makeSymmetric(X);
    while(check == false) {
        check = arma::chol(R, X);
        //if(check == false) {
        //    check = arma::lu(L, U, P, X);
        //    R = sqrt(U % L.t());
        //}
        if(check == false) {
            check = chol_fb(R, X);
        }
        if(check == false) {
            Rcout << X << "\n";
            X += arma::eye(X.n_rows, X.n_rows) * 1e-7;
        }
    }
    return R;
}


IntegerVector tabulate_cpp (arma::vec const & x, int const m){
    IntegerVector counts(m);
    for (auto& now : x) {
        if (now > 0 && now <= m)
            counts[now - 1]++;
    }
    return counts;
}

IntegerVector tabulate2(arma::vec const & x){
    int m, nbins;
    nbins =  max(x);
    if (nbins > 1) {
        m = nbins;
    } else {
        m = 1;
    }
    return tabulate_cpp(x, m);
}

IntegerVector unique_sort (arma::vec const & x) {
    arma::vec xu = arma::sort(arma::unique(x));
    IntegerVector xuv = IntegerVector(xu.begin(), xu.end());
    return xuv;
}

int length_unique (arma::vec const & x) {
    int m;
    arma::vec xu = arma::unique(x);
    //int m = arma::max(x);
    m = xu.n_elem;
    return m;
}

IntegerVector concat_iv(IntegerVector a, IntegerVector b) {
    int a_s = a.size(), b_s = b.size();
    IntegerVector v(a_s + b_s);
    for(int i = 0; i < a_s + b_s; i++) {
        if (i < a_s) {
            v(i) = a[i];
        } else {
            int b_i = i - a_s;
            v(i) = b[b_i];
        }
    }
    return v;
}

const double log2pi = std::log(2.0 * M_PI);


arma::vec dmvnorm(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false) {
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(chol2(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
    }

    if (logd == false) {
        out = exp(out);
    }
    return(out);
}


// [[Rcpp::export]]
arma::mat rmvnorm (int n, arma::vec mean, arma::mat sigma){
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mean, 1, n).t() + Y * chol2(sigma);
}


// [[Rcpp::export]]
arma::mat rwishart(int df, const arma::mat& S){
  // Dimension of returned wishart
   int m = S.n_rows;

  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);

  // Fill the diagonal
  for( int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }

  // Fill the lower matrix with random guesses
  for( int j = 0; j < m; j++){
    for( int i = j+1; i < m; i++){
      Z(i,j) = R::rnorm(0,1);
    }
  }

  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * chol2(S);

  // Return random wishart
  return C.t()*C;
}

arma::cube nrwishart (int n, int df, const arma::mat& S) {
    int m = S.n_rows;
    arma::cube cb(m, m, n);
    for( int i = 0; i < n; i++){
        cb.slice(i) = rwishart(df, S);
    }
    return cb;
}


// [[Rcpp::export]]
arma::mat ginv(arma::mat m, Rcpp::Nullable<double> tol = R_NilValue) {
    arma::mat invM;
    bool check;
    if (tol.isNotNull()) {
        double i = Rcpp::as<double>(tol);
        check = arma::pinv(invM, m, i);
    } else {
        check = arma::pinv(invM, m);
    }
    if (check) {
        return invM;
    } else {
        invM = 1;
        return invM;
    }
}


IntegerVector get_clusters(const arma::mat &traces, int tot_clusters, double discard = 0.1) {
    int N = traces.n_rows, S = traces.n_cols;
    double d = N * discard;
    int keep_out = std::round(d);
    if (keep_out >= N) {
        keep_out = N - 1;
    }
    arma::mat z = traces.rows(keep_out, N - 1);
    IntegerVector clusters(S), sample_tabulate(tot_clusters);
    for (int i = 0; i < S; i++) {
        sample_tabulate = tabulate_cpp(z.col(i), tot_clusters);
        clusters[i] = arma::index_max(Rcpp::as<arma::vec>(sample_tabulate)) + 1;
    }
    return(clusters);
}

arma::mat get_clusters_stats(const arma::mat &traces, int tot_clusters, double discard = 0.1) {
    int N = traces.n_rows, S = traces.n_cols;
    double d = N * discard;
    int keep_out = std::round(d);
    arma::mat z = traces.rows(keep_out, N - 1), clusters(S, tot_clusters);
    IntegerVector sample_tabulate(tot_clusters);
    for (int i = 0; i < S; i++) {
        sample_tabulate = tabulate_cpp(z.col(i), tot_clusters);
        clusters.row(i) = Rcpp::as<arma::rowvec>(sample_tabulate) / (N - keep_out);
    }
    return(clusters);
}
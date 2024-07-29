#include "utils.h"


void _sample_z(int u, arma::mat A, double alpha, arma::vec &z, arma::mat &mu_ar, arma::cube &S_ar, arma::vec mu_a0, arma::mat R_a0, double beta_a0, arma::mat W_a0, int &nclust){

  // number of auxiliary tables to approximate the infinite number of empty tables
  int m = 5, K = 0;

  IntegerVector n_tot = tabulate2(z);
  // Reuse the computation from tabulate to find
  // the number and the labels of the non zeros
  // clusters. Make an IntegerVector w/o zeros
  // clusters (simlar to table())
  for (int i = 0; i < n_tot.size(); i++) {
    if (n_tot[i] > 0) {
      K += 1;
    }
  }
  IntegerVector n(K), clust_labels(K);
  for (int i = 0, i2 = 0; i < n_tot.size(); i++) {
    if (n_tot[i] > 0) {
      n[i2] = n_tot[i];
      clust_labels[i2] = i + 1;
      i2 += 1;
    }
  }

  int D = A.n_rows;

  arma::mat R_a0_inv = R_a0.i();
  arma::mat W_a0_inv = W_a0.i();
  int z_u_k, z_u = z(u);
  double logp;
  // Compute likelihood for every table
  arma::vec logprobs(K+m);

  for(int ki = 0; ki < K; ki++ ){
    int k = clust_labels[ki] - 1;
    z_u_k = 0;
    if (z_u == (k + 1)) {
       z_u_k = 1;
    }

    logp = log(n(ki) - z_u_k);
    logp += dmvnorm(trans(A.col(u)), trans(mu_ar.col(k)), ginv(S_ar.slice(k), R_NilValue), true).at(0,0);
    logprobs(ki) = logp;
  }

  // Create m auxiliary tables from the base distribution
  arma::mat aux_tables_mean = rmvnorm(m, mu_a0, R_a0_inv).t();
  int df;
  if (beta_a0 > D) {
    df = int(beta_a0 + 0.5);
  } else {
    df = D;
  }

  arma::cube aux_tables_covariance = nrwishart(m, df, W_a0 / beta_a0);

  // If last point in cluster, re-use its parameters for the auxiliary table so that
  // it has some probability of staying 
  if (n_tot(z_u - 1) == 1){
    aux_tables_mean.col(0) = mu_ar.col(z_u - 1);
    aux_tables_covariance.slice(0) = S_ar.slice(z_u - 1);
  }

  // Compute likelihood for auxilary table
  for(int k = 0; k < m; k++){
    logp = log(alpha/m);
    logp += dmvnorm(trans(A.col(u)), trans(aux_tables_mean.col(k)), aux_tables_covariance.slice(k), true).at(0,0);
    logprobs[K + k] = logp;
  }
  arma::vec probs = logprobs - max(logprobs);
  probs = exp(probs);
  probs = probs / sum(probs);
  
  // Choose table
  IntegerVector aux_labels =  Rcpp::Range(nclust+1, nclust+m);
  IntegerVector clusters = concat_iv(clust_labels, aux_labels);


  IntegerVector chosen2 = sample_r(clusters, 1, probs);
  int chosen = chosen2[0];

  // Re-label tables
  // if auxiliary table is chosen, assign label K+i-th to this new table
  if(chosen > nclust){
    mu_ar.resize(mu_ar.n_rows, nclust + 1);
    S_ar.resize(S_ar.n_rows, S_ar.n_cols, nclust + 1);
    mu_ar.col(nclust) = aux_tables_mean.col(chosen - nclust - 1);
    S_ar.slice(nclust) = aux_tables_covariance.slice(chosen - nclust - 1).i();
    nclust += 1;
    z(u) = nclust;
  }
  else{
    z(u) = chosen;
  }
}

// [[Rcpp::export]]
Rcpp::List sample_z(int u, arma::mat A, double alpha, arma::vec z, arma::mat mu_ar, arma::cube S_ar, arma::vec mu_a0, arma::mat R_a0, double beta_a0, arma::mat W_a0, int nclust){
  u = u - 1;
  _sample_z(u, A, alpha, z, mu_ar, S_ar, mu_a0, R_a0,  beta_a0, W_a0, nclust);

  return(Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("nclust") = nclust,
    Rcpp::Named("mu_ar") = mu_ar,
    Rcpp::Named("S_ar") = S_ar)
  );
}

// [[Rcpp::export]]
Rcpp::List sample_z_loop(int N, arma::mat A, double alpha, arma::vec z, arma::mat mu_ar, arma::cube S_ar, arma::vec mu_a0, arma::mat R_a0, double beta_a0, arma::mat W_a0, int nclust){

  for (int n = 0; n < N; n++){
    _sample_z(n, A, alpha, z, mu_ar, S_ar, mu_a0, R_a0,  beta_a0, W_a0, nclust);
  }
  return(Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("nclust") = nclust,
    Rcpp::Named("mu_ar") = mu_ar,
    Rcpp::Named("S_ar") = S_ar)
  );
}

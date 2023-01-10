// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "commonvars.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

void InitChains(arma::cube &X,
                arma::vec &alpha, 
                arma::mat &p, arma::mat &z,
                arma::cube &means, 
                arma::field<arma::cube> &covs) {
  double prod_b;
  NumericVector probs(num_comps);
  IntegerVector ints = seq_len(num_comps);
  IntegerVector z_class;
  
  // X
  X.slice(0) = X_mat;
  
  // Means
  arma::vec colMeans(dim);
  for (int i = 0; i < dim; i++) {
    colMeans(i) = mean(X_mat.col(i));
  }
  
  for (int i = 0; i < num_comps; i++) {
    means.slice(0).col(i) = colMeans; // set each column to the mean of data
    covs(0).slice(i) = cov(X_mat); // set each matrix to the covariance of the data
  }
  
  alpha(0) = 1;
  NumericVector b = rbeta(num_comps, 1, alpha(0));
  p(0,0) = b(0);
  for (int i = 1; i < num_comps; i++) {
    prod_b = 1;
    for (int j = 0; j < i; j++) {
      prod_b *= (1-b(j));
    }
    p(0, i) = b(i) * prod_b;
  }
  
  probs = NumericVector(p.row(0).begin(), p.row(0).end());
  z_class = sample(ints, n, true, probs);
  z.row(0) = as<arma::rowvec>(wrap(z_class)); 
}
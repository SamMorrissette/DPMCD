// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
#include "commonvars.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

Params DrawUnequalUnrestricted(arma::rowvec z) {
  Params parameters;
  arma::cube covariance(dim, dim, num_comps);
  arma::mat mean(dim, num_comps); 
  arma::uvec pos;
  arma::vec colMeans(dim);
  arma::mat x_j, centered_x, W_k(dim, dim), T(dim, dim);
  int n_k;
  
  for (int k = 0; k < num_comps; k++) {
    pos = find(z == k+1);
    n_k = pos.n_elem;
    if (n_k == 0) {
      covariance.slice(k) = iwishrnd(S, nu0);
      mean.col(k) = mvnrnd(mu0, kappa*covariance.slice(k), 1);
      continue;
    }
    
    x_j = X_mat.rows(pos);
    centered_x.set_size(n_k, dim);
    W_k.fill(0);
    
    for (int d = 0; d < dim; d++) {
      colMeans(d) = arma::mean(x_j.col(d));
    }
    
    for (int i = 0; i < n_k; i++) {
      for (int j = 0; j < dim; j++) {
        centered_x(i, j) = x_j(i, j) - colMeans(j);
      }
    }
    for (int i = 0; i < n_k; i++) {
      W_k += centered_x.row(i).t() * centered_x.row(i);
    }
    double d_nk = (double) n_k;
    T = (kappa*d_nk / (d_nk+kappa)) * (colMeans * colMeans.t());
    W_k += T + S;
    covariance.slice(k) = iwishrnd(W_k, nu0+n_k);
    
    arma::vec pst_mean = (d_nk  * colMeans + (kappa*mu0)) / (d_nk + kappa);
    arma::mat pst_var = covariance.slice(k) / (d_nk  + kappa);
    mean.col(k) = mvnrnd(pst_mean, pst_var, 1);
  }
  parameters.mean = mean;
  parameters.covariance = covariance;
  return parameters;
}


Params DrawUnequalSpherical(arma::rowvec z) {
  Params parameters;
  arma::cube covariance(dim, dim, num_comps);
  arma::mat mean(dim, num_comps); 
  arma::uvec pos;
  arma::vec colMeans(dim);
  arma::mat x_j, centered_x, temp_mat;
  int n_k;
  double lambda_k, W_k, T, pst_IG_alpha, pst_IG_beta;
  
  for (int k = 0; k < num_comps; k++) {
    pos = find(z == k+1);
    n_k = pos.n_elem;
    if (n_k == 0) {
      lambda_k = 1.0 / (R::rgamma(prior_IG_alpha, 1.0/prior_IG_beta));
      covariance.slice(k) = lambda_k * arma::mat(dim, dim, arma::fill::eye);
      mean.col(k) = mvnrnd(mu0, covariance.slice(k) / kappa, 1);
      continue;
    }
    
    x_j = X_mat.rows(pos);
    centered_x.set_size(n_k, dim);
    W_k = 0;
    
    for (int d = 0; d < dim; d++) {
      colMeans(d) = arma::mean(x_j.col(d));
    }
    
    for (int i = 0; i < n_k; i++) {
      for (int j = 0; j < dim; j++) {
        centered_x(i, j) = x_j(i, j) - colMeans(j);
      }
    }
    for (int i = 0; i < n_k; i++) {
      temp_mat = centered_x.row(i).t() * centered_x.row(i);
      W_k += sum(temp_mat.diag());
    }
    
    double d_nk = (double) n_k;
    T = as<double>(wrap((kappa * d_nk / (d_nk+kappa)) * (colMeans.t() * colMeans)));
    pst_IG_alpha = prior_IG_alpha + ((d_nk * dim) / 2);
    pst_IG_beta = prior_IG_beta + ((W_k + T) / 2);
    
    lambda_k = 1.0 / (R::rgamma(pst_IG_alpha, 1.0/pst_IG_beta));
    covariance.slice(k) = lambda_k * arma::mat(dim, dim, arma::fill::eye);
    
    arma::vec pst_mean = (d_nk  * colMeans + (kappa*mu0)) / (d_nk + kappa);
    arma::mat pst_var = covariance.slice(k) / (d_nk  + kappa);
    mean.col(k) = mvnrnd(pst_mean, pst_var, 1);
  }
  parameters.mean = mean;
  parameters.covariance = covariance;
  return parameters;
}




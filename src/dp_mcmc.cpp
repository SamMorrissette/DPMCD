// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <progress.hpp>
#include "helpers.h"
#include "commonvars.h"

arma::mat X_mat;
int n, dim, num_comps, num_iters;
double sigmasq_prior_shape;
double sigmasq_prior_scale;
const double prior_shape = 1, prior_rate = 1;
arma::vec mu0;
arma::mat S;
double nu0, kappa;
double prior_IG_alpha;
double prior_IG_beta;


using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
List DP_MCMC(arma::mat obs_dist, 
             arma::mat init_X, double init_sigmasq, 
             const int K, const int iters,
             int modelIndex) {
  
  // Global variables used in other files
  X_mat = init_X;
  n = init_X.n_rows;
  dim = init_X.n_cols;
  num_comps = K;
  num_iters = iters;
  
  // Helper variables
  NumericVector b;
  Params newParams;
  double s_res = CalcSSR(distRcpp(wrap(init_X)), wrap(obs_dist));
  
  // Priors for NIW (UU model)
  mu0 = arma::vec(dim, arma::fill::zeros);
  nu0 = dim+2;
  //S = arma::mat(dim, dim, arma::fill::eye);
  //S = ((nu0 - dim - 1) * arma::cov(init_X)) / 4;
  kappa = 1;
  S = arma::cov(init_X);
  
  // Priors for IG (US model)
  arma::mat cov_X_init = cov(X_mat);
  prior_IG_alpha = (dim+2) / 2;
  prior_IG_beta = (sum(cov_X_init.diag()) / dim) / 2; //(dim+2) / 2; 
  
  // Priors for sigma squared (Inverse-Gamma prior)
  sigmasq_prior_shape = 5;
  sigmasq_prior_scale = (sigmasq_prior_shape-1) * (s_res / (n*(n-1)/2));
  
  // Declare Chains
  arma::cube X(n, dim, iters);
  arma::vec sigmasq(iters);
  arma::vec alpha(iters);
  arma::mat p(iters, K);
  arma::mat z(iters, n);
  arma::cube class_probs(n, num_comps, iters);
  arma::cube means(dim, num_comps, iters);
  arma::field<arma::cube> covs(iters);
  covs.fill(arma::zeros<arma::cube>(dim, dim, num_comps));
  
  // Label-switching variables
  arma::mat Q(n, num_comps);
  arma::uvec perms;
  arma::mat perm_mat(iters, num_comps, arma::fill::zeros);
  arma::cube old_probs(n, num_comps, iters, arma::fill::zeros);
  
  // Initialize the chains
  sigmasq(0) = init_sigmasq;
  InitChains(X, alpha, p, z, means, covs);
  
  Progress prog(iters, true);
  
  // Begin Gibbs
  for (int t = 1; t < iters; t++) {
    if (Progress::check_abort() )
      return -1.0;
    prog.increment();
    
    // Step 1 (Update X)
    X.slice(t) = UpdateX(obs_dist, X.slice(t-1), sigmasq(t-1), z.row(t-1), means.slice(t-1), covs(t-1));
    
    // Step 1.5 (Transform X)
    X.slice(t) = TransformX(X.slice(t), init_X);
    
    X_mat = X.slice(t); // Update global X matrix
    
    // Step 2 (Update sigma_sq measurement error)
    sigmasq(t) = UpdateSigmasq(obs_dist, X.slice(t), s_res, sigmasq(t-1));
    
    // Step 3 (drawing classes for each data point)
    class_probs.slice(t) = CalculateClassProbs(p.row(t-1), means.slice(t-1), covs(t-1));
    old_probs.slice(t) = class_probs.slice(t);
    
    //Step 4 (Label switching)
    if (t == 100) {
      Q = InitializeQ(class_probs);
    }  else if (t > 100) {
      perms = UndoLabelSwitching(Q, class_probs.slice(t));
      perm_mat.row(t) = arma::conv_to<arma::rowvec>::from(perms);
      for (int k = 0; k < num_comps; k++) {
        class_probs.slice(t).col(k) = old_probs.slice(t).col(perms(k)-1);
      }
      Q = UpdateQ(Q, class_probs.slice(t), t);
    }
    
    // Step 5 (update allocations) 
    z.row(t) = UpdateClasses(class_probs.slice(t));

    // Step 6 (update stick-breaking weights) 
    b = StickBreaking(z.row(t), alpha(t-1));
    p.row(t) = UpdateWeights(b);
    
    // Step 7 (Update model parameters)
    newParams = UpdateTheta(z.row(t), modelIndex);
    means.slice(t) = newParams.mean;
    covs(t) = newParams.covariance;

    // Step 8 (Update alpha)
    alpha(t) = UpdateAlpha(b);
    //alpha(t) = 1;
  }
  
  return List::create(
    Named("X") = X,
    Named("sigmasq") = sigmasq,
    Named("alpha") = alpha,
    Named("means") = means,
    Named("covs") = covs,
    Named("z") = z,
    Named("class_probs") = class_probs, 
    Named("p") = p,
    Named("perms") = perm_mat);
}

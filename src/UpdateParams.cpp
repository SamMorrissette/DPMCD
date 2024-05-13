// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
#include "commonvars.h"
#include "helpers.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


arma::mat CalculateClassProbs(arma::rowvec p,
                              arma::mat mean,
                              arma::cube cov) {
  IntegerVector ints = seq_len(num_comps);
  arma::rowvec probs(num_comps);
  double numerator, denom;
  arma::rowvec out_vec(n);
  arma::mat prob_mat(n, num_comps);
  arma::vec result(1);
  
  for (int i = 0; i < n; i++) {
    denom = 0;
    for (int k = 0; k < num_comps; k++) {
      result = dmvnorm(X_mat.row(i), mean.col(k), cov.slice(k));
      denom = denom + p(k) * result(0);
      }
    
    for (int k = 0; k < num_comps; k++) {
      result = dmvnorm(X_mat.row(i), mean.col(k), cov.slice(k));
      numerator = p(k) * result(0);
      probs(k) = numerator/denom;
      }
    prob_mat.row(i) = probs;
    }
  return prob_mat;
}

arma::rowvec UpdateClasses(arma::mat prob_mat) {
  IntegerVector ints = seq_len(num_comps);
  NumericVector probs(num_comps);
  arma::rowvec out_vec(n);
  for (int i = 0; i < n; i++) {
    probs = as<NumericVector>(wrap(prob_mat.row(i)));
    out_vec(i) = as<int>(wrap(sample(ints, 1, true, probs)));
  }
  return out_vec;
}

NumericVector StickBreaking(arma::rowvec z, double alpha) {
  NumericVector out_vec(num_comps);
  //Rprintf("%f \n", alpha);
  int n_k;
  double total = 0;
  double gamma1, gamma2;
  arma::uvec pos;
  for (int k = 0; k < (num_comps-1); k++) {
    pos = find(z == k+1);
    n_k = pos.n_elem;
    total += n_k;
    
 
    
    out_vec(k) = R::rbeta(1+n_k, alpha + (n-total));
    
    
    
    
    //gamma1 = as<double>(wrap(arma::randg(arma::distr_param(1+n_k, 1))));
    //gamma2 = as<double>(wrap(arma::randg(arma::distr_param(alpha+(n-total), 1.0))));
    //out_vec(k) = gamma1 / (gamma1+gamma2);
    //Rprintf("%f out_vec \n", out_vec(k));
  }
  out_vec(num_comps-1) = 1;

  return out_vec;
}

arma::rowvec UpdateWeights(NumericVector b) {
  arma::rowvec out_vec(num_comps);
  double prod_b;
  out_vec(0) = b(0);
  for (int k = 1; k < num_comps; k++) {
    prod_b = 1;
    for (int j = 0; j < k; j++) {
      prod_b *= (1-b(j));
    }
    out_vec(k) = b(k) * prod_b;
  }
  return out_vec;
}

Params UpdateTheta(arma::rowvec z, int modelIndex) {
  Params params_out;
  switch (modelIndex) {
    case 1: // Unequal Unrestricted
      {
        params_out = DrawUnequalUnrestricted(z);
        break;
      }
    case 2: // Unequal Spherical
      {
        params_out = DrawUnequalSpherical(z);
        break;
      }
    case 3: // Unequal Diagonal;
      {
        params_out = DrawUnequalDiagonal(z);
        break;
      }
    case 4: // Equal Unrestricted;
      { 
        params_out = DrawEqualUnrestricted(z);
        break;
      }
    case 5: // Equal Spherical;
      {
        params_out = DrawEqualSpherical(z);
        break;
      }
    case 6: // Equal Diagonal;
      {
        params_out = DrawEqualDiagonal(z);
        break;
      }
  }
  return params_out;
}

double UpdateAlpha(NumericVector b) {
  double alpha, sum_logs, pst_rate;
  alpha = 1;
  // sum_logs = 0;
  // 
  // for (int i = 0; i < (num_comps-1); i++) {
  //   sum_logs += log(1-b(i));
  // }
  // //Rprintf("%f \n", sum_logs);
  // 
  // pst_rate = prior_rate - sum_logs;
  // //when most of the observations are allocated to a single cluster, then the stick weight, b(i),
  // //will be close to 1 for that cluster. 
  // 
  // //when b(i) is close to 1, then sum_logs is large and negative which makes pst_rate very high and alpha very small
  // 
  // //when alpha is very small, then the stick breaking weight for a component might be 1 (i.e. b(i) = 1) and then when 
  // //you calculate log(1-1) = -Inf
  // 
  // //
  // try{
  //   alpha = as<double>(wrap(arma::randg(1, arma::distr_param(prior_shape + num_comps - 1, 1/pst_rate))));
  // } catch(...) {
  //   // for (int i = 0; i < (num_comps-1); i++) {
  //   //   Rprintf("%f \n", b(i));
  //   // }
  //   Rprintf("sum_logs %f \n", sum_logs);
  //   Rprintf("pst_rate %f \n", pst_rate);
  //   warning("Error!");
  //   stop("error!");
  //   alpha = 1;
  // }
  return alpha;
}







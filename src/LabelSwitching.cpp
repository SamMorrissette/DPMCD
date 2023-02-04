// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "commonvars.h"
#include "RcppHungarian.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppHungarian)]]

arma::rowvec NormalizeRow(arma::rowvec x) {
  const double UPPER_THRESH = 0.999999;
  const double LOWER_THRESH = 0.000001;
  double row_sum;
  arma::rowvec norm_x(num_comps);
  
  for (int i = 0; i < num_comps; i++) {
    if (x(i) > UPPER_THRESH) {
      x(i) = UPPER_THRESH;
    } else if (x(i) < LOWER_THRESH) {
      x(i) = LOWER_THRESH;
    }
  }
  row_sum = sum(x);
  for (int i = 0; i < num_comps; i++) {
    norm_x(i) = x(i)/row_sum;
  }

  return norm_x;
}

arma::mat InitializeQ(arma::cube class_probs) {
  arma::mat Q(n, num_comps);
  int switch_iters = 100;
  
  for (int m = 0; m < switch_iters; m++) {
    for (int i = 0; i < n; i++) {
      class_probs.slice(m).row(i) = NormalizeRow(class_probs.slice(m).row(i));
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < num_comps; j++) {
      for (int m = 0; m < switch_iters; m++) {
        Q(i, j) += class_probs.slice(m)(i, j);
      }
    }
  }
  Q = Q / (switch_iters-0);
  return Q;
}

// This function was taken from the Rcpp hungarian package! 
IntegerVector HungarianSolver(NumericMatrix costMatrix) {
  int nr = costMatrix.nrow();
  int nc = costMatrix.ncol();
  
  vector<double> c(nc);
  vector<vector<double>> cm(nr, c);
  for (int i=0; i < nr; i++){
    for (int j=0; j < nc; j++){
      c[j] = costMatrix(i,j);
    }
    cm[i] = c;
  }
  
  HungarianAlgorithm HungAlgo;
  vector<int> assignment;
  double cost = HungAlgo.Solve(cm, assignment);
  IntegerVector assign(nr);
  for (int i=0; i < nr; i++){
    assign(i) = assignment[i]+1;
  }
  return assign;
}

arma::uvec UndoLabelSwitching(arma::mat Q, arma::mat prob_mat) {
  arma::mat cost_mat(num_comps, num_comps, arma::fill::zeros);
  NumericMatrix hung_input;
  IntegerVector hung_output;
  arma::uvec perm_vec(num_comps);
  
  // Normalize first:
  for (int i = 0; i < n; i++) {
    prob_mat.row(i) = NormalizeRow(prob_mat.row(i));
  }
  
  // Build cost matrix
  for (int j = 0; j < num_comps; j++) {
    for (int l = 0; l < num_comps; l++) {
      for (int i = 0; i < n; i++) {
        cost_mat(j, l) += prob_mat(i, l) * (log(prob_mat(i, l)) - log(Q(i, j)));
      }
    }
  }
  
  // Solve using RcppHungarian
  hung_input = as<NumericMatrix>(wrap(cost_mat));  
  hung_output = HungarianSolver(hung_input);
  perm_vec = as<arma::uvec>(wrap(hung_output));
  
  return perm_vec;
}

arma::mat UpdateQ(arma::mat Q, arma::mat perm_prob_mat, int t) {
  arma::mat new_Q(num_comps, num_comps, arma::fill::zeros);
  t = (double) (t - 100);
  new_Q = ((t*Q) + perm_prob_mat) / (t+1);
  return new_Q;
}
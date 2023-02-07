
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
//#define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// [[Rcpp::export]]
double CalcLPML(List dpobj) {
  int B, n;
  arma::uword z_i;
  double LPML = 0, total;
  arma::vec alpha = dpobj["alpha"], density;
  arma::mat bmds_X = dpobj["BMDS_X"], Z = dpobj["z"];
  arma::cube X = dpobj["X"], means = dpobj["means"];
  arma::field<arma::cube> covs = dpobj["covs"];
  B = alpha.n_elem;
  n = bmds_X.n_rows;
  
  for (int i = 0; i < n; i++) {
    total = 0;
    for (int s = 0; s < B; s++) {
      z_i = Z(s, i) - 1;
      density = dmvnorm(X.slice(s).row(i), means.slice(s).col(z_i), covs(s).slice(z_i));
      total = total + (1 / density(0));
    }
    LPML = LPML + log(1.0 / ((1.0 / (double) B) * total));
  }
  return LPML;
}

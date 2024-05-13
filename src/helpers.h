#ifndef __UTILITIES__
#define __UTILITIES__

typedef struct Params Params;

arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false);

double CalcSSR(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2);

Rcpp::NumericMatrix distRcpp(Rcpp::NumericMatrix X);

void InitChains(arma::cube &X,
                arma::vec &alpha, 
                arma::mat &p, arma::mat &z, 
                arma::cube &means, 
                arma::field<arma::cube> &covs);

arma::mat UpdateX(arma::mat obs_distances, arma::mat current_X, 
                  double sigmasq, arma::rowvec z,
                  arma::mat means, arma::cube covs);

arma::mat TransformX(arma::mat X, arma::mat X_star);

double UpdateSigmasq(arma::mat obs_distances, arma::mat current_X, double s_res, double sigmasq);

arma::mat CalculateClassProbs(arma::rowvec p,
                                 arma::mat mean,
                                 arma::cube cov);

arma::rowvec UpdateClasses(arma::mat prob_mat);

Rcpp::NumericVector StickBreaking(arma::rowvec z, double alpha);
arma::rowvec UpdateWeights(Rcpp::NumericVector b);

double UpdateAlpha(Rcpp::NumericVector b);

// Updating GMM parameters
Params UpdateTheta(arma::rowvec z, int modelIndex);
Params DrawUnequalUnrestricted(arma::rowvec z);
Params DrawUnequalSpherical(arma::rowvec z);
Params DrawUnequalDiagonal(arma::rowvec z);
Params DrawEqualUnrestricted(arma::rowvec z);
Params DrawEqualSpherical(arma::rowvec z);
Params DrawEqualDiagonal(arma::rowvec z);

// Label switching
arma::mat InitializeQ(arma::cube class_probs);
arma::uvec UndoLabelSwitching(arma::mat Q, arma::mat prob_mat);
arma::mat UpdateQ(arma::mat Q, arma::mat prob_mat, int t);


#endif // __UTILITIES__
extern arma::mat X_mat;
extern int n, dim, num_comps, num_iters;
extern double sigmasq_prior_shape, sigmasq_prior_scale;
extern const double prior_shape, prior_rate;
extern arma::vec mu0;
extern arma::mat S;
extern double nu0, kappa;
extern double prior_IG_alpha, prior_IG_beta; // Priors for spherical models
typedef struct Params {
  arma::mat mean;
  arma::cube covariance;
} Params;

RunDPM <- function(distances, init_X, init_sigmasq, K, burn, iters, modelIndices) {
  
  # Run MCMC
  for (i in 1:length(modelIndices)) {
    dpobj <- DP_MCMC(distances, init_X, init_sigmasq, K, iters, modelIndices[i])
  }
  
  # Get rid of burn-in
  ind <- (burn+1):iters
  out <- list(BMDS_X = init_X,
              X = dpobj$X[,,ind],
              sigmasq = dpobj$sigmasq[ind,],
              alpha = dpobj$alpha[ind,],
              means = dpobj$means[,,ind],
              covs = dpobj$covs[ind,],
              z = dpobj$z[ind,],
              class_probs = dpobj$class_probs[,,ind],
              p = dpobj$p[ind,],
              perms = dpobj$perms)
  return(out)
}

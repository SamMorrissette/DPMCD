createOutput <- function(dpobj, init_X, burn, iters) {
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
              perms = dpobj$perms[ind,])
  return(out)
}

RunDPM <- function(distances, init_X, init_sigmasq, K, burn, iters, modelIndices, parallel, cores, seed) {
  if (parallel == FALSE) {
    # Run MCMC one model at a time
    allModels <- vector("list", length = length(modelIndices))
    for (i in 1:length(modelIndices)) {
      dpobj <- DP_MCMC(distances, init_X, init_sigmasq, K, iters, modelIndices[i])
      allModels[[i]] <- createOutput(dpobj, init_X, burn, iters)
    }
  } else if (parallel == TRUE) {
    if (!foreach::getDoParRegistered()) {
      doParallel::registerDoParallel(cores=cores)
    }
    allModels <- foreach::foreach(j=1:length(modelIndices), .packages ="DPMCD") %dopar% {
      set.seed(seed*2)
      dpobj <- DP_MCMC(distances, init_X, init_sigmasq, K, iters, modelIndices[j])
      output <- createOutput(dpobj, init_X, burn, iters)
    }
  }
  
  return(allModels)
}

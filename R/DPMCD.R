DPMCD <- function(distances, max_p, parallel, K = 25, 
                  burn = 0, iters = 5000, 
                  modelNames = c("UU", "US", "EU")) {
  
  match.arg(modelNames, several.ok = TRUE)
  
  allModels <- c("UU", "US", "EU")
  modelIndices <- match(modelNames, allModels)
  
  # Initial run of BMDS for initialization and dimension estimation
  BMDS_out <- RunBMDS(distances, max_p, parallel=parallel) 
  plot(BMDS_out$mdsics)
  p <- 3 #which.min(BMDS_out$mdsics)
  X_est <- BMDS_out$X[[p]]
  sigmasq_est <- BMDS_out$sigma_sq[[p]]
  
  # Run Dirichlet Process Mixture Model
  DPM_out <- RunDPM(distances, X_est, sigmasq_est, K, burn, iters, modelIndices)
  DPM_out
}
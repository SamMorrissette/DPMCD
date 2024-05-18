DPMCD <- function(distances, max_p, K = 25, 
                  burn = 1000, iters = 5000, 
                  modelNames = c("UU", "US", "UD", "EU", "ES", "ED"),
                  parallel = FALSE, cores) {
  
  match.arg(modelNames, several.ok = TRUE)
  
  if ((isTRUE(parallel) & missing(cores))) {
    stop("Please provide the number of cores for parallelization (must be an integer greater than 0)")
  } else if (isFALSE(parallel) & !missing(cores)) {
    stop("Number of cores should not be specified when parallel is set to FALSE")
  }  
  
  if (missing(cores)) {
    cores = 0
  }

  
  allModels <- c("UU", "US", "UD", "EU", "ES", "ED")
  modelIndices <- match(modelNames, allModels)
  
  # Initial run of BMDS for initialization and dimension estimation
  BMDS_out <- RunBMDS(distances, max_p, 
                      parallel = parallel, cores = cores)
  p <- which.min(BMDS_out$mdsics)
  X_est <- BMDS_out$X[[p]]
  sigmasq_est <- BMDS_out$sigma_sq[[p]]
  print(paste("Dimension:", p))
  print("Done running BMDS, starting DPM")
  
  # Run Dirichlet Process Mixture Model
  DPM_out <- RunDPM(distances, X_est, sigmasq_est, K, burn, iters, modelIndices, 
                    parallel = parallel, cores = cores)
  
  # Add model names to each element of the list
  for (i in 1:length(modelNames)) {
    DPM_out[[i]]$modelName <- modelNames[i]
  }
  
  DPM_out
}
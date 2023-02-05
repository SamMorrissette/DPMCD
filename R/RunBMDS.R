#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @export

RunBMDS <- function(distances, max_p, parallel = FALSE, cores) {
  bmds_burn = 1000
  bmds_iter = 5000
  if(parallel == TRUE) {
    doParallel::registerDoParallel(cores=cores)
  }
  n <- nrow(distances)
  
  X <- vector("list", length = max_p)
  sigma_sq <- vector("list", length = max_p)
  
  if (parallel == FALSE) {
    for (i in 1:max_p) {
      # Doesn't include delta matrix (reduce memory required)
      temp_bmds <- bmdsMCMC(DIST = distances, p = i, nwarm = bmds_burn, niter = bmds_iter)
      X[[i]] <- temp_bmds$x_bmds
      sigma_sq[[i]] <- temp_bmds$e_sigma
      print(i)
    }
  } else if (parallel == TRUE) {
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    out_list <- foreach::foreach(j=1:max_p, .packages ="DPMCD", .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
      output <- bmdsMCMC(distances, j, nwarm = bmds_burn, niter = bmds_iter)
      list(output$x_bmds, output$e_sigma)
    }
    X <- out_list[[1]]
    sigma_sq <- out_list[[2]]
  }
  
  if(parallel == TRUE) {
    doParallel::stopImplicitCluster()
  }
  
  mdsics <- MDSIC(distances, X)
  
  return(list(X = X,
              sigma_sq = sigma_sq,
              mdsics = mdsics))
}

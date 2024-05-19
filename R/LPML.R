LPML <- function(dp_object) {
  num_models <- length(dp_object)
  LPML_vec <- rep(NA, num_models)
  for (i in 1:num_models) {
    # This is to ensure proper conversion to armadillo cubes
    dim = ncol(dp_object[[i]]$BMDS_X)
    if (dim == 1) {
      dp_object[[i]]$X <- simplify2array(lapply(seq_len(ncol(dp_object[[i]]$X)), 
                                                function(t) dp_object[[i]]$X[,t,drop=FALSE]))
      dp_object[[i]]$means <- simplify2array(lapply(seq_len(ncol(dp_object[[i]]$means)), 
                                                    function(t) t(dp_object[[i]]$means[,t])))
      }


    LPML_vec[i] <- CalcLPML(dp_object[[i]])
  }
  return(LPML_vec)
}
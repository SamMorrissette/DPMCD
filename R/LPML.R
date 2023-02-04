LPML <- function(dp_object) {
  num_models <- length(dp_object)
  LPML_vec <- rep(NA, num_models)
  for (i in 1:num_models) {
    LPML_vec[[i]] <- CalcLPML(dp_object[[i]])
  }
  return(LPML_vec)
}
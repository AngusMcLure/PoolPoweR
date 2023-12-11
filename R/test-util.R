# To address issue [#27](https://github.com/AngusMcLure/PoolPoweR/issues/27).
# Ensures that relative tolerance is used between small values when testing for
# equality.
numeric_rel_diff <- function(act, exp, tolerance) {
  rd <- abs(act-exp)/abs(exp)
  if(rd < tolerance) TRUE
  else {
    warning(glue::glue("Relative difference: {rd} is >= tolerance: {tolerance}"))
    FALSE
  }
}

matrix_rel_diff <- function(act, exp, tolerance) {
  cm <- matrix(c(1,1,1,1), nrow = 2)
  comp <- all.equal(act/exp, cm)
  #comp <- mean(abs(act-exp)/mean(exp))
  if(is.logical(comp)) return(TRUE) # All values are perfectly equal
  else {
    mean_rel_diff <- as.numeric(strsplit(comp, " ")[[1]][4])
    if(mean_rel_diff < tolerance) return(TRUE)
    else {
      warning(glue::glue("mean_rel_diff: {mean_rel_diff} is >= tolerance: {tolerance}"))
      return(FALSE)
    }
  }
}

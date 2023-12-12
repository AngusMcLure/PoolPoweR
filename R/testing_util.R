# To address issue [#27](https://github.com/AngusMcLure/PoolPoweR/issues/27).
# Ensures that relative tolerance is used between small values when testing for
# equality.
numeric_rel_diff <- function(act, exp, tolerance, var_name=FALSE) {
  rd <- abs(act-exp)/abs(exp)
  if(rd < tolerance) return(TRUE)
  else {
    if(is.character(var_name)) message(glue::glue("In {var_name}:"))
    message(glue::glue("Relative difference: {rd} is >= tolerance: {tolerance}"))
    return(FALSE)
  }
}

matrix_rel_diff <- function(act, exp, tolerance) {
  mean_rd <- mean(abs(act-exp)/mean(abs(exp)))
  if(mean_rd < tolerance) return(TRUE)
  else {
    message(glue::glue("Mean relative difference: {mean_rd} is >= tolerance: {tolerance}"))
    return(FALSE)
  }
}

list_rel_diff <- function(act, exp, tolerance) {
  var_names <- names(act)
  comp <- mapply(numeric_rel_diff, act, exp, var_names,
                 MoreArgs = list(tolerance = tolerance))
  if(any(!comp)) return(FALSE)
  return(TRUE)
}



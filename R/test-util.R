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
  mean_rd <- mean(abs(act-exp)/mean(exp))
  if(mean_rd < tolerance) TRUE
  else {
    warning(glue::glue("Mean relative difference: {mean_rd} is >= tolerance: {tolerance}"))
    FALSE
  }
}

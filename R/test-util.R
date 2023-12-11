relative_difference <- function(actual, expected, tolerance) {
  # https://github.com/AngusMcLure/PoolPoweR/issues/27
  # Always uses relative tolerance when testing for equality
  # all.equal() switches to absolute tolerance when values are small
  cm <- matrix(c(1,1,1,1), nrow = 2)
  comp <- all.equal(actual/expected, cm)
  if(class(comp) == "logical") return(TRUE) # All values are perfectly equal
  else {
    mean_rel_diff <- as.numeric(strsplit(comp, " ")[[1]][4])
    if(mean_rel_diff < tolerance) return(TRUE)
    else {
      warning(glue::glue("mean_rel_diff: {mean_rel_diff} is >= tolerance: {tolerance}"))
      return(FALSE)
    }
  }
}

check_geq2 <- function(val, min) {
  name <- deparse(substitute(val)) # get name of variable
  if (!is.numeric(val)) {
    stop(glue::glue("{name} must be numeric, not {class(val)}."))
  }
  if (val < min) {
    stop(glue::glue("{name} must be >= {min}."))
  }
}

check_in_range2 <- function(val) {
  name <- deparse(substitute(val))
  msg <- message(glue::glue("{name} must be a numeric value between 0 and 1, inclusive."))
  if(!is.numeric(val)) {
    stop(glue::glue("{name} must be numeric, not {class(val)}."))
  }
  if(val < 0) {
    msg
    stop(glue::glue("{val} is < 0"))
  } else {
    msg
    if(val > 1) stop(glue::glue("{val} is > 1"))
  }
}
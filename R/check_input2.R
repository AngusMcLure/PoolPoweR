#' Check input helpers
#'
#' Helper functions to validate input parameters.
#' 
#' @param val numeric/NA Input parameter
#' @param min numeric for check_geq2
#' @param allow_na boolean whether NAs are acceptable inputs or not. i.e. for 
#' clustering cases correlation/cost_cluster are NA
check_geq2 <- function(val, min, allow_na = FALSE) {
  name <- deparse(substitute(val)) # get name of variable
  if (allow_na) {
  	if (!is.numeric(val) & !is.na(val)) {
  	  stop(glue::glue("{name} must be numeric or NA, not {class(val)}."))
  	}
  } else { # allow_na = FALSE
  	if (!is.numeric(val)) {
  	  stop(glue::glue("{name} must be numeric, not {class(val)}."))
  	}
  }
  if (is.numeric(val) && val < min) {
    stop(glue::glue("{name} must be >= {min}."))
  }
}

#' @rdname check_geq2
check_in_range2 <- function(val, allow_na = FALSE) {
  name <- deparse(substitute(val)) # get name of variable
  if (allow_na) {
  	if (!is.numeric(val) & !is.na(val)) {
  	  stop(glue::glue("{name} must be numeric or NA, not {class(val)}."))
    }
  } else {
  	if (!is.numeric(val)) {
  	  stop(glue::glue("{name} must be numeric, not {class(val)}."))
  	}
  }
  if (is.numeric(val) && (val < 0 || val > 1)) {
    message(glue::glue("{name} must be a numeric value between 0 and 1, inclusive."))
    stop(glue::glue("{name} = {val}"))
  }
}
